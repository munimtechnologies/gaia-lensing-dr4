#!/usr/bin/env python3
"""Gaia DR4 release-day job of gaia-lensing-dr4: weigh the lens stars of the astrometric
microlensing events Klueter, Bastian & Wambsganss (2020) predicted inside Gaia's DR4 window.

MV2, Munim, Inc. (mv2@munimtech.com). Apache-2.0. Made with AI agents working in MV2's name.
Frozen with PREREGISTRATION.md on 2026-11-24, before DR4 (2026-12-02).

Steps per event (events.py): DR2 astrometry of lens and source -> DR4 counterparts by cone
search at J2017.5 -> pair check (ID_MISMATCH if the DR4 tracks never come within 3 arcsec, or not
within 1 yr of the paper's T_CA) -> the source's DR4 epoch astrometry by DataLink -> per-transit
normal points (gaiaorb.epochs, unchanged) -> lensfit.fit_event at r_cut 500 mas (primary) and
1000 mas (robustness) -> verdict, mass or 95% upper limit, comparison with the paper's
forecast. Gate G0 (geometry) and the validation events (Proxima Cen, LAWD 37) decide the
overall label (lensfit.overall).
Usage: python lens_job.py --data-dir /data --out-dir /out [--release "Gaia DR4"] [--only 12,37]
"""
import argparse
import csv
import json
import math
import os
import time

import events as EV
import lensfit as L

T0 = time.time()
TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
DR4_TABLES = ("gaiadr4.gaia_source", "gaiadr4.all_source_astrometry")
LENS_RADIUS_MAS = 2000.0
SOURCE_RADIUS_MAS = 1000.0
PAIR_MAX_SEP_MAS = 3000.0
PAIR_MAX_DT_YR = 1.0


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def tap(query, tries=4):
    import io_gaia
    last = None
    for k in range(tries):
        try:
            return io_gaia.tap_csv(query, TAP)
        except Exception as e:  # archive busy on release day: back off and retry
            last = e
            log("TAP retry %d: %s" % (k + 1, str(e)[:200]))
            time.sleep(20 * (k + 1))
    raise RuntimeError("TAP failed: %s" % last)


def propagate(ra, dec, pmra, pmdec, from_ep, to_ep):
    """Linear propagation (flat sky; fine for a few years and arcsec-level matching)."""
    dt = to_ep - from_ep
    ddec = (pmdec or 0.0) * dt / L.MAS_PER_DEG
    dra = (pmra or 0.0) * dt / L.MAS_PER_DEG / max(1e-6, math.cos(math.radians(dec)))
    return ra + dra, dec + ddec


def dr2_lookup(ids):
    """{dr2_id: dict(ra, dec, pmra, pmdec, parallax, g, ref_epoch)} from DR2, else DR3 via the
    dr2_neighbourhood table."""
    idl = ",".join(str(i) for i in sorted(ids))
    out = {}
    try:
        rows = tap("SELECT source_id, ra, dec, parallax, pmra, pmdec, phot_g_mean_mag, ref_epoch "
                   "FROM gaiadr2.gaia_source WHERE source_id IN (%s)" % idl)
        for r in rows:
            out[int(r["source_id"])] = dict(ra=num(r["ra"]), dec=num(r["dec"]), pmra=num(r["pmra"]),
                                            pmdec=num(r["pmdec"]), parallax=num(r["parallax"]),
                                            g=num(r["phot_g_mean_mag"]), ref_epoch=num(r["ref_epoch"]) or 2015.5,
                                            via="gaiadr2")
    except Exception as e:
        log("DR2 lookup failed: %s" % e)
    missing = [i for i in ids if i not in out]
    if missing:
        try:
            rows = tap("SELECT n.dr2_source_id, g.ra, g.dec, g.parallax, g.pmra, g.pmdec, g.phot_g_mean_mag "
                       "FROM gaiadr3.dr2_neighbourhood AS n JOIN gaiadr3.gaia_source AS g "
                       "ON g.source_id = n.dr3_source_id WHERE n.dr2_source_id IN (%s)"
                       % ",".join(str(i) for i in missing))
            for r in rows:
                i = int(r["dr2_source_id"])
                if i not in out:
                    out[i] = dict(ra=num(r["ra"]), dec=num(r["dec"]), pmra=num(r["pmra"]),
                                  pmdec=num(r["pmdec"]), parallax=num(r["parallax"]),
                                  g=num(r["phot_g_mean_mag"]), ref_epoch=2016.0, via="gaiadr3")
        except Exception as e:
            log("DR3 fallback failed: %s" % e)
    return out


def dr4_match(pred_ra, pred_dec, radius_mas, g_ref=None, tables=DR4_TABLES):
    """Nearest DR4 entry to the predicted J2017.5 position; prefers |dG| < 1 when both G known."""
    for tbl in tables:
        q = ("SELECT * FROM %s WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', %.9f, %.9f, %.9f))"
             % (tbl, pred_ra, pred_dec, radius_mas / L.MAS_PER_DEG))
        try:
            rows = tap(q)
        except Exception as e:
            log("cone on %s failed: %s" % (tbl, e))
            continue
        cands = []
        for r in rows:
            ra, dec = num(r.get("ra")), num(r.get("dec"))
            if ra is None or dec is None:
                continue
            ep = num(r.get("ref_epoch")) or L.DR4_REF_JYEAR
            ra5, dec5 = propagate(ra, dec, num(r.get("pmra")), num(r.get("pmdec")), ep, L.DR4_REF_JYEAR)
            off = L.tangent_offset_mas(ra5, dec5, pred_ra, pred_dec)
            g = num(r.get("phot_g_mean_mag"))
            bad_g = g_ref is not None and g is not None and abs(g - g_ref) >= 1.0
            cands.append((bad_g, math.hypot(*off), dict(
                source_id=int(r["source_id"]), ra=ra5, dec=dec5, parallax=num(r.get("parallax")),
                parallax_error=num(r.get("parallax_error")), pmra=num(r.get("pmra")),
                pmdec=num(r.get("pmdec")), g=g, table=tbl, match_mas=math.hypot(*off))))
        if cands:
            cands.sort(key=lambda c: (c[0], c[1]))
            return cands[0][2]
    return None


def fetch_points(source_id, data_dir, release):
    import io_gaia
    from gaiaorb import epochs
    d = os.path.join(data_dir, "epochs")
    for url in io_gaia.datalink_urls(source_id, release):
        for k in range(3):
            try:
                blob = io_gaia.fetch(url, timeout=300)
                break
            except Exception as e:
                log("DataLink retry %d %s: %s" % (k + 1, source_id, str(e)[:160]))
                blob = None
                time.sleep(10 * (k + 1))
        if not blob:
            continue
        paths = io_gaia.save_epochs(blob, d, "src_%d" % source_id, log)
        if not paths:
            continue
        rows = []
        for p in paths:
            rows.extend(io_gaia.read_table_rows(p))
        by = epochs.rows_by_source(rows)
        mine = by.get(source_id) or (rows if len(by) == 1 else [])
        if mine:
            pts, counts = epochs.normal_points(mine)
            return pts, counts
    return None, None


def analyse_event(e, dr2, data_dir, release):
    rec = dict(num=e["num"], lens=e["lens"], lens_dr2=str(e["lens_dr2"]), source_dr2=str(e["source_dr2"]),
               t_ca_paper=e["t_ca"], m_in=e["m_in"], sm5=e["sm5"], sm10=e["sm10"])
    lo, so = dr2.get(e["lens_dr2"]), dr2.get(e["source_dr2"])
    if lo is None or lo["ra"] is None:
        rec["verdict"] = "LENS_NOT_FOUND_DR2"
        return rec, None
    if so is None or so["ra"] is None:
        rec["verdict"] = "SOURCE_NOT_FOUND_DR2"
        return rec, None
    lra, ldec = propagate(lo["ra"], lo["dec"], lo["pmra"], lo["pmdec"], lo["ref_epoch"], L.DR4_REF_JYEAR)
    sra, sdec = propagate(so["ra"], so["dec"], so["pmra"], so["pmdec"], so["ref_epoch"], L.DR4_REF_JYEAR)
    lens4 = dr4_match(lra, ldec, LENS_RADIUS_MAS, lo["g"])
    if lens4 is None or lens4["parallax"] is None or lens4["pmra"] is None:
        rec["verdict"] = "LENS_NOT_FOUND_DR4"
        return rec, None
    src4 = dr4_match(sra, sdec, SOURCE_RADIUS_MAS, so["g"])
    if src4 is None:
        rec["verdict"] = "SOURCE_NOT_FOUND_DR4"
        return rec, None
    rec.update(lens_dr4=str(lens4["source_id"]), source_dr4=str(src4["source_id"]),
               lens_g=lens4["g"], source_g=src4["g"], lens_plx=lens4["parallax"],
               lens_plx_err=lens4["parallax_error"], source_table=src4["table"],
               source_catalog_plx=src4["parallax"], source_match_mas=src4["match_mas"])
    geom = dict(alpha=src4["ra"], delta=src4["dec"],
                lens_offset=L.tangent_offset_mas(lens4["ra"], lens4["dec"], src4["ra"], src4["dec"]),
                lens_pm=(lens4["pmra"], lens4["pmdec"]), lens_plx=lens4["parallax"])
    src0 = (0.0, 0.0, src4["parallax"] or 0.0, src4["pmra"] or 0.0, src4["pmdec"] or 0.0)
    tmin, smin = L.closest_approach(geom, src0, -3.5, 3.5)
    rec.update(t_min_dr4=L.DR4_REF_JYEAR + tmin, sep_min_dr4_mas=smin)
    if smin > PAIR_MAX_SEP_MAS or abs(L.DR4_REF_JYEAR + tmin - e["t_ca"]) > PAIR_MAX_DT_YR:
        rec["verdict"] = "ID_MISMATCH"
        return rec, None
    pts, counts = fetch_points(src4["source_id"], data_dir, release)
    if not pts:
        rec["verdict"] = "NO_EPOCHS"
        return rec, None
    rec["epoch_counts"] = counts
    resid = L.geometry_residuals(geom["alpha"], geom["delta"], pts)
    rec["geometry_median_au"] = L.median(resid)
    fr = L.fit_event(pts, geom, src0)
    fr_alt = L.fit_event(pts, geom, src0, r_cut=L.R_CUT_ALT_MAS)
    rec["fit"] = {k: v for k, v in fr.items()}
    rec["fit_alt"] = {k: v for k, v in fr_alt.items()}
    rec["verdict"] = L.event_verdict(fr)
    rec["verdict_alt"] = L.event_verdict(fr_alt)
    if fr.get("status") == "OK":
        plx_err = lens4["parallax_error"] or 0.0
        pred = L.predicted_thetaE2(e["m_in"], lens4["parallax"], fr["plx_S"])
        rec.update(thetaE2=fr["thetaE2"], thetaE2_err=fr["thetaE2_err"], snr=fr["snr"],
                   thetaE2_predicted=pred, snr_predicted=pred / fr["thetaE2_err"] if fr["thetaE2_err"] > 0 else None)
        M, sM = L.mass_from_fit(fr, lens4["parallax"], plx_err)
        rec.update(mass=M, mass_err=sM)
        if rec["verdict"] != "DEFLECTION_DETECTED":
            rec["mass_upper95"] = L.upper_limit_mass(fr, lens4["parallax"])
        if fr_alt.get("status") == "OK" and M is not None:
            Ma, sMa = L.mass_from_fit(fr_alt, lens4["parallax"], plx_err)
            rec.update(mass_alt=Ma, mass_alt_err=sMa)
            if Ma is not None and sM is not None and sMa is not None:
                # pre-registered: ROBUST if the two masses agree within the larger 1-sigma error
                rec["robust"] = abs(M - Ma) <= max(sM, sMa)
        return rec, (pts, geom, fr)
    return rec, None


def write_points(path, pts, geom, fr):
    src = [fr["params"][k][0] for k in L.PARAMS[:5]]
    D = L.separation_track(pts, geom, src)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_yr_from_J2017.5", "w_mas", "sig_mas", "sin_psi", "cos_psi", "parallax_factor_al",
                    "sep_mas", "lens_column_per_mas", "used_primary"])
        for p, d in zip(pts, D):
            r = math.hypot(*d)
            g = (d[0] * p["s"] + d[1] * p["c"]) / (r * r) if r > 0 else 0.0
            w.writerow(["%.6f" % p["t"], "%.4f" % p["w"], "%.4f" % p["sig"], "%.6f" % p["s"], "%.6f" % p["c"],
                        "%.6f" % p["f"], "%.1f" % r, "%.8f" % g, int(r >= L.R_CUT_MAS)])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--release", default="Gaia DR4")
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    os.makedirs(a.data_dir, exist_ok=True)
    os.makedirs(os.path.join(a.out_dir, "points"), exist_ok=True)
    evs = EV.as_dicts()
    only = {int(x) for x in a.only.split(",") if x.strip()}
    if only:
        evs = [e for e in evs if e["num"] in only]
    ids = {e["lens_dr2"] for e in evs} | {e["source_dr2"] for e in evs}
    dr2 = dr2_lookup(ids)
    log("DR2/DR3 astrometry for %d of %d objects" % (len(dr2), len(ids)))
    recs = []
    geo = []
    for e in evs:
        try:
            rec, extra = analyse_event(e, dr2, a.data_dir, a.release)
        except Exception as ex:  # one bad event must not stop the rest
            rec, extra = dict(num=e["num"], lens=e["lens"], verdict="ERROR", error=str(ex)[:500]), None
        if extra is not None:
            write_points(os.path.join(a.out_dir, "points", "event_%02d.csv" % e["num"]), *extra)
        if rec.get("geometry_median_au") is not None:
            geo.append(rec["geometry_median_au"])
        log("event %d %s: %s %s" % (e["num"], e["lens"], rec.get("verdict"),
                                   "" if rec.get("mass") is None else "M = %.3f +- %.3f" % (rec["mass"], rec["mass_err"])))
        recs.append(rec)
    g0 = bool(geo) and L.median(geo) < L.GEOMETRY_TOL
    val = []
    for r in recs:
        if r["num"] in EV.VALIDATION and r.get("verdict") == "DEFLECTION_DETECTED" and r.get("mass") is not None:
            mref, sref = EV.VALIDATION[r["num"]]
            r["validation"] = L.validation_check(r["mass"], r["mass_err"], mref, sref)
            val.append(r["validation"])
    head = [r.get("verdict") for r in recs if r["num"] in EV.HEADLINE]
    outcome = L.overall(g0, val, head)
    res = dict(project="gaia-lensing-dr4", lensfit_version=L.VERSION, release=a.release,
               gate_G0=dict(passed=g0, median_au=L.median(geo) if geo else None, n_events=len(geo)),
               validation=val, outcome=outcome, events=recs, runtime_s=round(time.time() - T0, 1))
    with open(os.path.join(a.out_dir, "results.json"), "w") as f:
        json.dump(res, f, indent=1, default=str)
    cols = ["num", "lens", "verdict", "verdict_alt", "lens_dr4", "source_dr4", "source_g", "lens_plx",
            "sep_min_dr4_mas", "t_min_dr4", "thetaE2", "thetaE2_err", "snr", "thetaE2_predicted",
            "snr_predicted", "mass", "mass_err", "mass_upper95", "mass_alt", "mass_alt_err", "m_in", "sm5",
            "validation"]
    with open(os.path.join(a.out_dir, "events.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in recs:
            w.writerow([r.get(c, "") for c in cols])
    log("OUTCOME %s" % outcome)
    return res


if __name__ == "__main__":
    main()
