# Pre-registration: weighing nearby stars with the microlensing events Gaia was predicted to see

Frozen 2026-11-24 (UTC), eight days before Gaia DR4 (2026-12-02), by AI agents working in the
name of MV2 of Munim, Inc. (mv2@munimtech.com). No DR4 data for any event in this test existed
when it was written. Licence: Apache-2.0 (code), CC BY 4.0 (text). No human expert has checked it.

## 1. The near miss

When a nearby star passes close in front of a distant one, its gravity pushes the distant star's
image outward by a fraction of a milliarcsecond. The size of the push gives the mass of the
nearby star directly, without any model of its interior. This is the only direct way to weigh a
single star (or white dwarf, or brown dwarf) that has no companion.

- McGill et al. (2018, arXiv:1804.07049) predicted that the white dwarf LAWD 37 would lens a
  background star in November 2019 and that Gaia itself could measure its mass. HST measured it
  instead: 0.56 +- 0.08 Msun (McGill et al. 2023, MNRAS; arXiv:2206.01814).
- Klueter, Bastian & Wambsganss (2020, A&A; arXiv:1911.02584) simulated Gaia's
  individual measurements for 501 predicted events and listed the events whose lens mass Gaia
  should be able to measure, with the expected error for 5 and 10 years of data. They state
  that the individual measurements needed "will come with DR4".
- Gaia DR4 (66 months, 2014 July 25 to 2020 January 20) publishes, for the first time, the
  epoch astrometry of every source, faint background stars included (ESA DR4 content page).
- As far as we can find (arXiv listing search 2026-11-24 for "astrometric microlensing"; the
  newest related preprints, Kaczmarek et al. 2026 arXiv:2609.10705 and Wicker et al. 2026
  arXiv:2611.24561, prepare for DR4 with simulated events and Gaia photometric-alert events), no
  one has published or pre-registered the mass measurements for these predicted events.

## 2. Question

For the events Klueter et al. (2020) predicted with closest approach before 2020.5, does the Gaia
DR4 epoch astrometry of the background star show the predicted outward deflection, and what
lens masses does it give?

## 3. Data and event list

- `events.py`: all 43 events of the paper's Tables 1 and 2 with T_CA < 2020.5, transcribed with
  their Gaia DR2 identifiers, the paper's assumed mass M_in and expected errors (sM5, sM10).
- **Headline set** (fixed now): the events whose expected 5-year error sM5 / M_in is at most 20%
  (#12 DR2 5312099874809857024 with M_in = 0.07 Msun, #17 L 51-47, #20 G 217-32,
  #21 DR2 5865259639247544448, #22 HD 146868, #26 HD 149192), plus #37 LAWD 37 (peak 2019.865,
  inside DR4).
- **Validation set** (independently known masses): Proxima Cen #15 and #30, 0.1221 +- 0.0022 Msun
  (Kervella, Thevenin & Lovis 2017); LAWD 37 #37, 0.56 +- 0.08 Msun (McGill et al. 2023).
- Lens and source: DR2 astrometry (gaiadr2.gaia_source; DR3 via dr2_neighbourhood if needed),
  propagated to J2017.5 and matched in gaiadr4.gaia_source, then gaiadr4.all_source_astrometry
  (lens within 2", source within 1"; nearest, preferring |dG| < 1). Pair check: if the DR4
  tracks never come within 3" of each other, or the closest approach is more than 1 yr from the
  paper's T_CA, the event is ID_MISMATCH and not fitted.
- Source epoch astrometry by DataLink (io_gaia 0.2.0 of gj3512-giants-dr4), collapsed to one
  along-scan point per transit by `gaiaorb.epochs.normal_points` (unchanged, as used by MV2's
  other DR4 tests).

## 4. Model and rules (lensfit.py 0.1.0)

- Along-scan model: the source's 5 astrometric parameters plus thetaE^2 times
  g = (D . scan direction) / |D|^2, where D is the unlensed source-minus-lens position. Linear least
  squares, jitter added until chi2/dof <= 1, one 5-sigma clip, geometry iterated 3 times. Lens
  position, proper motion and parallax are the DR4 catalogue values.
- Mass M = thetaE^2 / (8.144 mas/Msun x (plx_lens - plx_source)), errors propagated to first order
  (thetaE^2, source parallax and their covariance, lens parallax error).
- **Primary cut:** transits with |D| < 500 mas are dropped (Gaia's windows of the bright lens and
  the faint source overlap there). **Robustness:** the fit is repeated at 1000 mas; a mass is
  ROBUST if the two agree within the larger 1-sigma error.
- **Per-event verdict:** TOO_FEW (< 12 transits after the cut); DEFLECTION_DETECTED (thetaE^2 / its
  error >= 3); NEGATIVE_DEFLECTION (<= -3: the source is pulled toward the lens, the signature of
  light from the lens contaminating the source, not lensing); otherwise NOT_DETECTED with a 95%
  upper limit on the mass. Also: LENS/SOURCE_NOT_FOUND (DR2 or DR4), ID_MISMATCH, NO_EPOCHS.
- **Gate G0 (geometry):** the median over events of the median |model - published| AL parallax
  factor must be below 0.03 au. If it fails the outcome is GEOMETRY_FAIL and nothing is claimed.
- **Validation:** each validation event detected at S/N >= 3 PASSES if
  |M - M_ref| <= 3 sqrt(sM^2 + sM_ref^2). Any FAIL labels the whole run SYSTEMATICS (masses are
  reported but not claimed); no validation event measured labels it UNVALIDATED.
- **Overall outcome:** `<VALIDATED|SYSTEMATICS|UNVALIDATED>_<n>_HEADLINE_DETECTED`.
- **The forecast tested** (descriptive): for every fitted event we report the measured error of
  thetaE^2 against the predicted thetaE^2 = 8.144 M_in plx_rel, i.e. whether Gaia DR4 reaches the
  precision Klueter et al. forecast.

## 5. Before release day

`rehearsal.py` runs on ESA's DR4 pre-release epoch astrometry (12 sources, 2026-06-29): R1 checks
the ephemeris and sign conventions (the best-matching sky position must reproduce the published
AL parallax factors to rms < 0.03 au, median, and < 0.05 for every source); R2 is a phantom-lens
null (|S/N| < 3 for at least 9 of 12). If R1 fails, the geometry code may be corrected before DR4
in a dated amendment; the rules above do not change.

## 6. What would be claimed

Only DEFLECTION_DETECTED events in a VALIDATED or UNVALIDATED run, stated with their S/N and the
robustness check. A NOT_DETECTED headline event is reported as an upper limit and as a forecast
that DR4 did not reach. McGill, Everall, Boubert & Smith (2020, MNRAS Letters; arXiv:2006.13958)
showed that many predicted bright events are false because the "source" is a duplicate
detection or a binary companion of the lens; several headline sources had only 2-parameter DR2
solutions. Such events show up here as SOURCE_NOT_FOUND_DR4 or ID_MISMATCH and are reported as
a test of that critique, not as failures of the fit.

## Credit

The predictions are Klueter, Bastian, Demleitner & Wambsganss (2018, A&A 620, A175) and Klueter,
Bastian & Wambsganss (2020); LAWD 37: McGill et al. (2018, 2023); Proxima: Kervella et al. (2017);
the false-event critique: McGill et al. (2020).
Gaia data: ESA/Gaia/DPAC.
