#!/usr/bin/env python3
"""Pre-DR4 rehearsal of gaia-lensing-dr4 on ESA's Gaia DR4 pre-release epoch astrometry
(12 sources, published 2026-06-29). Run BEFORE any DR4 event data exist.

MV2, Munim, Inc. (mv2@munimtech.com). Apache-2.0. Made with AI agents working in MV2's name.

R1 geometry: for every pre-release source, find the sky position whose model AL parallax
   factors (lensfit.observer_xyz + parallax_offsets) best match the published
   parallax_factor_al. PASS if the median rms over sources is below lensfit.GEOMETRY_TOL
   (0.03 au) and no source exceeds 0.05. This checks the ephemeris, the L2 factor and every sign
   convention the release-day fit relies on.
R2 null: fit each source with a phantom lens (offset (1500, 300) mas, relative motion
   (-600, -100) mas/yr, parallax 200 mas; no real lens is there). PASS if at least 9 of the 12
   sources give |S/N| < 3 (the pre-release includes three orbital binaries, which may not).
If R1 fails, lensfit's geometry may be fixed before DR4 (a dated amendment); nothing else
changes. Usage: python rehearsal.py --data-dir /data --out-dir /out
"""
import argparse
import json
import os
import time

import lensfit as L

T0 = time.time()
PHANTOM = dict(lens_offset=(1500.0, 300.0), lens_pm=(-600.0, -100.0), lens_plx=200.0)


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def grade(rms_list, snr_list):
    rms_ok = [r for r in rms_list if r is not None]
    r1 = bool(rms_ok) and L.median(rms_ok) < L.GEOMETRY_TOL and max(rms_ok) <= 0.05
    n_quiet = sum(1 for s in snr_list if s is not None and abs(s) < 3.0)
    r2 = n_quiet >= 9
    return dict(R1_geometry="PASS" if r1 else "FAIL", R2_null="PASS" if r2 else "FAIL",
                n_sources=len(rms_list), n_quiet=n_quiet,
                median_rms=L.median(rms_ok) if rms_ok else None,
                max_rms=max(rms_ok) if rms_ok else None,
                outcome="REHEARSAL_PASS" if (r1 and r2) else "REHEARSAL_FAIL")


def main(argv=None):
    import io_gaia
    from gaiaorb import epochs
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    a = ap.parse_args(argv)
    os.makedirs(a.data_dir, exist_ok=True)
    os.makedirs(a.out_dir, exist_ok=True)
    zp = io_gaia.download(io_gaia.PRERELEASE_URLS, os.path.join(a.data_dir, "dr4-prerelease.zip"), log=log)
    paths = io_gaia.extract(zp, os.path.join(a.data_dir, "prerelease"))
    tables = io_gaia.usable_tables([p for p in paths if p.lower().rsplit(".", 1)[-1] in io_gaia.TABLE_EXTS], log)
    log("usable epoch tables: %d" % len(tables))
    per_source = {}
    for p in tables:
        rows = io_gaia.read_table_rows(p)
        for sid, rs in epochs.rows_by_source(rows).items():
            key = sid if sid is not None else os.path.basename(p)
            per_source.setdefault(key, []).extend(rs)
    results = []
    for sid in sorted(per_source, key=str):
        pts, counts = epochs.normal_points(per_source[sid])
        rec = dict(source=str(sid), n_points=len(pts), counts=counts)
        if len(pts) < L.MIN_POINTS:
            rec.update(status="TOO_FEW", rms=None, snr=None)
            results.append(rec)
            continue
        alpha, delta, rms = L.locate_from_parallax_factors(pts)
        geom = dict(PHANTOM, alpha=alpha, delta=delta)
        fr = L.fit_event(pts, geom, (0.0, 0.0, 0.0, 0.0, 0.0))
        rec.update(status=fr["status"], alpha_fit=alpha, delta_fit=delta, rms=rms,
                   snr=fr.get("snr"), thetaE2=fr.get("thetaE2"), thetaE2_err=fr.get("thetaE2_err"),
                   n_used=fr.get("n_points"), jitter_mas=fr.get("jitter_mas"),
                   params=fr.get("params"))
        log("source %s: %d points, position (%.3f, %.3f), rms %.4f au, phantom S/N %s"
            % (sid, len(pts), alpha, delta, rms, rec["snr"]))
        results.append(rec)
    g = grade([r["rms"] for r in results], [r["snr"] for r in results])
    out = dict(project="gaia-lensing-dr4 rehearsal", lensfit_version=L.VERSION, grade=g, sources=results,
               runtime_s=round(time.time() - T0, 1))
    with open(os.path.join(a.out_dir, "rehearsal.json"), "w") as f:
        json.dump(out, f, indent=1)
    log("OUTCOME %s %s" % (g["outcome"], json.dumps(g)))
    return out


if __name__ == "__main__":
    main()
