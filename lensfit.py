"""Astrometric microlensing of a background star by a known nearby star, fitted to Gaia DR4
along-scan (AL) epoch astrometry. Standard library only.

Part of gaia-lensing-dr4 (MV2, Munim, Inc., mv2@munimtech.com). Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it yet.

Model (frozen 2026-11-24, before Gaia DR4):
  The source's AL abscissa per transit (Gaia DR4 / gaiasupdate convention, as in gaiaorb.model):
    w = da* s + dd c + plx_S (f - w fr k) + mua* t s + mud t c + thetaE2 * g
  with s = sin(psi), c = cos(psi), f the published AL parallax factor, t in Julian years from
  J2017.5, and the lensing column
    g = (D_a* s + D_d c) / |D|^2,
  D(t) = unlensed source position minus lens position (mas, tangent plane). For a resolved
  source (|D| >> thetaE) the bright lensed image is pushed AWAY from the lens by
  thetaE^2 / |D| (exact: thetaE (sqrt(u^2 + 4) - u) / 2 with u = |D| / thetaE; the relative
  difference is below 1/u^2). The model is linear in (da*, dd, plx_S, mua*, mud, thetaE2);
  D(t) depends on the source parameters, so the fit is iterated (geometry from the previous
  solution, 3 passes).
  Mass: thetaE^2 = KAPPA * M * (plx_L - plx_S), KAPPA = 4G/(c^2 au) = 8.144 mas / Msun.
Geometry: Gaia's barycentric position is taken as 1.01 x the Earth's heliocentric position
(Gaia orbits L2), with the Sun from the Astronomical Almanac low-precision formulae
(about 0.01 deg). Error budget: below 0.02 au in the parallax offsets, i.e. below 2% of
plx_rel in |D|, which changes g by under 4% at the separations used.
"""
import math

VERSION = "0.1.0"
KAPPA = 8.144              # mas per solar mass
MAS_PER_DEG = 3.6e6
MAS_TO_RAD = math.pi / (180.0 * 3600.0 * 1000.0)
DR4_REF_JYEAR = 2017.5
JD_J2000 = 2451545.0
L2_FACTOR = 1.01
PARAMS = ("ra_offset_mas", "dec_offset_mas", "parallax_mas", "pmra_masyr", "pmdec_masyr",
          "thetaE2_mas2")

# Pre-registered thresholds (PREREGISTRATION.md section 4)
R_CUT_MAS = 500.0          # primary: drop transits with lens-source separation below this
R_CUT_ALT_MAS = 1000.0     # robustness re-fit
MIN_POINTS = 12            # fewer usable transits after the cut -> TOO_FEW
SNR_DETECT = 3.0
GEOMETRY_TOL = 0.03        # au; median |f_model - f_published| for gate G0
CLIP = 5.0
N_ITER = 3


# ----------------------------------------------------------------------------- geometry
def jd_from_t(t):
    return JD_J2000 + (DR4_REF_JYEAR + t - 2000.0) * 365.25


def observer_xyz(t):
    """Gaia's approximate barycentric position (au, ICRS equatorial) at t (yr from J2017.5)."""
    n = jd_from_t(t) - JD_J2000
    L = math.radians((280.460 + 0.9856474 * n) % 360.0)
    g = math.radians((357.528 + 0.9856003 * n) % 360.0)
    lam = L + math.radians(1.915) * math.sin(g) + math.radians(0.020) * math.sin(2.0 * g)
    R = 1.00014 - 0.01671 * math.cos(g) - 0.00014 * math.cos(2.0 * g)
    eps = math.radians(23.439 - 0.0000004 * n)
    xs = R * math.cos(lam)
    ys = R * math.cos(eps) * math.sin(lam)
    zs = R * math.sin(eps) * math.sin(lam)
    return (-L2_FACTOR * xs, -L2_FACTOR * ys, -L2_FACTOR * zs)


def parallax_offsets(alpha_deg, delta_deg, xyz):
    """Apparent displacement per unit parallax (east a*, north d) for an observer at xyz (au)."""
    a, d = math.radians(alpha_deg), math.radians(delta_deg)
    X, Y, Z = xyz
    pa = X * math.sin(a) - Y * math.cos(a)
    pd = X * math.cos(a) * math.sin(d) + Y * math.sin(a) * math.sin(d) - Z * math.cos(d)
    return pa, pd


def model_parallax_factor(alpha_deg, delta_deg, p, xyz=None):
    xyz = observer_xyz(p["t"]) if xyz is None else xyz
    pa, pd = parallax_offsets(alpha_deg, delta_deg, xyz)
    return pa * p["s"] + pd * p["c"]


def geometry_residuals(alpha_deg, delta_deg, points):
    """|f_model - f_published| per point (gate G0)."""
    return [abs(model_parallax_factor(alpha_deg, delta_deg, p) - p["f"]) for p in points]


def median(xs):
    v = sorted(xs)
    n = len(v)
    if n == 0:
        return None
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def locate_from_parallax_factors(points, coarse_deg=2.0, fine_deg=0.1, span_deg=3.0):
    """Sky position whose model parallax factors best match the published ones (rehearsal
    check of the ephemeris and sign conventions when no catalogue position is at hand).
    Returns (alpha, delta, rms)."""
    xyz = [observer_xyz(p["t"]) for p in points]

    def rms(a, d):
        s = 0.0
        for p, o in zip(points, xyz):
            r = model_parallax_factor(a, d, p, o) - p["f"]
            s += r * r
        return math.sqrt(s / len(points))

    best = (1e9, 0.0, 0.0)
    d = -90.0 + coarse_deg / 2.0
    while d < 90.0:
        a = 0.0
        while a < 360.0:
            v = rms(a, d)
            if v < best[0]:
                best = (v, a, d)
            a += coarse_deg
        d += coarse_deg
    _, a0, d0 = best
    n = int(round(span_deg / fine_deg))
    for i in range(-n, n + 1):
        dd = max(-90.0, min(90.0, d0 + i * fine_deg))
        for j in range(-n, n + 1):
            aa = (a0 + j * fine_deg) % 360.0
            v = rms(aa, dd)
            if v < best[0]:
                best = (v, aa, dd)
    return best[1], best[2], best[0]


def tangent_offset_mas(ra_deg, dec_deg, ra0_deg, dec0_deg):
    """(a*, d) offset in mas of (ra, dec) from (ra0, dec0); flat-sky, fine below 1 arcmin."""
    da = (ra_deg - ra0_deg + 180.0) % 360.0 - 180.0
    return (da * math.cos(math.radians(dec0_deg)) * MAS_PER_DEG, (dec_deg - dec0_deg) * MAS_PER_DEG)


def separation_track(points, geom, src):
    """D = source - lens (mas) at each point. geom: alpha, delta (deg; source reference
    position), lens_offset (lens minus source reference position at J2017.5, mas),
    lens_pm (mas/yr), lens_plx (mas). src: (da*, dd, plx, mua*, mud) of the source."""
    out = []
    la, ld = geom["lens_offset"]
    lma, lmd = geom["lens_pm"]
    lp = geom["lens_plx"]
    for p in points:
        pa, pd = parallax_offsets(geom["alpha"], geom["delta"], observer_xyz(p["t"]))
        sa = src[0] + src[3] * p["t"] + src[2] * pa
        sd = src[1] + src[4] * p["t"] + src[2] * pd
        xa = la + lma * p["t"] + lp * pa
        xd = ld + lmd * p["t"] + lp * pd
        out.append((sa - xa, sd - xd))
    return out


def closest_approach(geom, src, t0=-3.0, t1=3.0, step_days=1.0):
    """(t_min, sep_min_mas) of the unlensed relative track over [t0, t1] yr from J2017.5."""
    best = (None, 1e300)
    n = int((t1 - t0) * 365.25 / step_days) + 1
    for i in range(n):
        t = t0 + i * step_days / 365.25
        D = separation_track([dict(t=t, s=0.0, c=1.0)], geom, src)[0]
        r = math.hypot(D[0], D[1])
        if r < best[1]:
            best = (t, r)
    return best


# ----------------------------------------------------------------------------- linear algebra
def solve_wls(cols, y, sig):
    """Weighted least squares. Returns dict(beta, cov, chi2, resid)."""
    m = len(cols)
    n = len(y)
    A = [[0.0] * m for _ in range(m)]
    b = [0.0] * m
    for i in range(n):
        w = 1.0 / (sig[i] * sig[i])
        row = [cols[a][i] for a in range(m)]
        for a in range(m):
            ra = row[a] * w
            b[a] += ra * y[i]
            for c in range(a, m):
                A[a][c] += ra * row[c]
    for a in range(m):
        for c in range(a):
            A[a][c] = A[c][a]
    cov = invert(A)
    beta = [sum(cov[a][c] * b[c] for c in range(m)) for a in range(m)]
    resid = [y[i] - sum(beta[a] * cols[a][i] for a in range(m)) for i in range(n)]
    chi2 = sum((resid[i] / sig[i]) ** 2 for i in range(n))
    return dict(beta=beta, cov=cov, chi2=chi2, resid=resid)


def invert(A):
    """Gauss-Jordan inverse with partial pivoting (small symmetric matrices)."""
    n = len(A)
    M = [list(map(float, A[i])) + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-300:
            raise ValueError("singular normal matrix")
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [v / pv for v in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0.0:
                fct = M[r][col]
                M[r] = [vr - fct * vc for vr, vc in zip(M[r], M[col])]
    return [row[n:] for row in M]


# ----------------------------------------------------------------------------- fit
def design(points, D):
    k = MAS_TO_RAD
    y = [p["w"] for p in points]
    cols = [[p["s"] for p in points],
            [p["c"] for p in points],
            [p["f"] - p["w"] * p.get("fr", 0.0) * k for p in points],
            [p["t"] * p["s"] for p in points],
            [p["t"] * p["c"] for p in points],
            [(d[0] * p["s"] + d[1] * p["c"]) / (d[0] * d[0] + d[1] * d[1]) for p, d in zip(points, D)]]
    return cols, y


def _sig(base, eps):
    return [math.sqrt(s * s + eps * eps) for s in base]


def fit_jitter(cols, y, sig0, n_par, iters=60):
    """Smallest eps >= 0 with chi2 / (N - n_par) <= 1 (bisection)."""
    target = len(y) - n_par
    if target <= 0:
        return 0.0

    def chi2(eps):
        return solve_wls(cols, y, _sig(sig0, eps))["chi2"]

    if chi2(0.0) <= target:
        return 0.0
    hi = max(1e-4, median(sig0))
    k = 0
    while chi2(hi) > target:
        hi *= 2.0
        k += 1
        if k > 80:
            raise ValueError("jitter search diverged")
    lo = 0.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if chi2(mid) > target:
            lo = mid
        else:
            hi = mid
    return hi


def fit_event(points, geom, src0, r_cut=R_CUT_MAS, n_iter=N_ITER, clip=CLIP, with_lens=True):
    """Fit the source's 5 astrometric parameters and thetaE^2 to its normal points.
    src0: starting (da*, dd, plx, mua*, mud) (catalogue values; offsets 0).
    Returns a dict; 'status' is OK or TOO_FEW."""
    src = list(src0)
    pts_all = list(points)
    fit = None
    used = []
    for _ in range(n_iter):
        D_all = separation_track(pts_all, geom, src)
        keep = [i for i, d in enumerate(D_all) if math.hypot(d[0], d[1]) >= r_cut]
        pts = [pts_all[i] for i in keep]
        D = [D_all[i] for i in keep]
        if len(pts) < MIN_POINTS:
            return dict(status="TOO_FEW", n_points=len(pts), n_total=len(pts_all))
        n_par = 6 if with_lens else 5
        n_clipped = 0
        for rnd in range(2):
            cols, y = design(pts, D)
            cols = cols[:n_par]
            base = [p["sig"] for p in pts]
            eps = fit_jitter(cols, y, base, n_par)
            sig = _sig(base, eps)
            fit = solve_wls(cols, y, sig)
            if rnd == 1:
                break
            good = [i for i in range(len(pts)) if abs(fit["resid"][i]) <= clip * sig[i]]
            if len(good) == len(pts) or len(good) < MIN_POINTS:
                break
            n_clipped = len(pts) - len(good)
            pts = [pts[i] for i in good]
            D = [D[i] for i in good]
        used = pts
        src = fit["beta"][:5]
    seps = [math.hypot(d[0], d[1]) for d in D]
    out = dict(status="OK", n_points=len(used), n_total=len(pts_all), n_clipped=n_clipped,
               jitter_mas=eps, chi2=fit["chi2"], r_cut_mas=r_cut,
               sep_min_used_mas=min(seps), sep_max_used_mas=max(seps),
               t_span=[min(p["t"] for p in used), max(p["t"] for p in used)],
               params={name: [fit["beta"][a], math.sqrt(fit["cov"][a][a])]
                       for a, name in enumerate(PARAMS[:len(fit["beta"])])})
    if with_lens:
        th, sth = fit["beta"][5], math.sqrt(fit["cov"][5][5])
        out.update(thetaE2=th, thetaE2_err=sth, snr=th / sth if sth > 0 else 0.0,
                   cov_thetaE2_plx=fit["cov"][5][2], plx_S=fit["beta"][2],
                   plx_S_err=math.sqrt(fit["cov"][2][2]))
    return out


def mass_from_fit(fr, plx_L, plx_L_err):
    """Lens mass (Msun) and 1-sigma error from thetaE^2 and the relative parallax."""
    th, sth = fr["thetaE2"], fr["thetaE2_err"]
    ps, sps, cov = fr["plx_S"], fr["plx_S_err"], fr["cov_thetaE2_plx"]
    rel = plx_L - ps
    if rel <= 0:
        return None, None
    M = th / (KAPPA * rel)
    dth = 1.0 / (KAPPA * rel)
    dps = th / (KAPPA * rel * rel)
    dpl = -dps
    var = dth * dth * sth * sth + dps * dps * sps * sps + 2.0 * dth * dps * cov + dpl * dpl * plx_L_err ** 2
    return M, math.sqrt(max(0.0, var))


def predicted_thetaE2(m_in, plx_L, plx_S):
    return KAPPA * m_in * (plx_L - plx_S)


def event_verdict(fr):
    """Pre-registered per-event verdict from a fit_event result."""
    if fr.get("status") != "OK":
        return fr.get("status", "FAILED")
    if fr["snr"] >= SNR_DETECT:
        return "DEFLECTION_DETECTED"
    if fr["snr"] <= -SNR_DETECT:
        return "NEGATIVE_DEFLECTION"
    return "NOT_DETECTED"


def upper_limit_mass(fr, plx_L):
    """95% one-sided upper limit on the mass (Msun) when not detected (thetaE^2 floored at 0)."""
    rel = plx_L - fr["plx_S"]
    if rel <= 0:
        return None
    return (max(0.0, fr["thetaE2"]) + 1.645 * fr["thetaE2_err"]) / (KAPPA * rel)


def validation_check(M, sM, m_ref, s_ref):
    """PASS if |M - m_ref| <= 3 sqrt(sM^2 + s_ref^2)."""
    return "PASS" if abs(M - m_ref) <= 3.0 * math.sqrt(sM * sM + s_ref * s_ref) else "FAIL"


def overall(g0_ok, validation_results, headline_verdicts):
    """Overall outcome. validation_results: list of PASS/FAIL for validation events measured
    with S/N >= 3. headline_verdicts: list of per-event verdicts of the headline set."""
    if not g0_ok:
        return "GEOMETRY_FAIL"
    if any(v == "FAIL" for v in validation_results):
        val = "SYSTEMATICS"
    elif validation_results:
        val = "VALIDATED"
    else:
        val = "UNVALIDATED"
    n_det = sum(1 for v in headline_verdicts if v == "DEFLECTION_DETECTED")
    return "%s_%d_HEADLINE_DETECTED" % (val, n_det)
