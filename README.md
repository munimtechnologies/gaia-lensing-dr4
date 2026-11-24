# gaia-lensing-dr4

**Can Gaia weigh the nearby stars whose gravity it was predicted to see bending starlight?**
A pre-registered test for Gaia DR4 (released 2026-12-02). By MV2 of Munim, Inc.
(mv2@munimtech.com). Licence: Apache-2.0 (code); CC BY 4.0 (text and result files). The code
and text were made with AI agents working in his name, and no human expert has checked them yet.

## The idea

When a nearby star passes a fraction of an arcsecond from a distant one, its gravity pushes the
distant star's image outward by up to a few milliarcseconds. The size of that push gives the
nearby star's mass directly. In 2018-2020 Klueter, Bastian, Demleitner & Wambsganss predicted dozens of
such events during Gaia's mission and forecast which lens masses Gaia's own measurements could
determine; McGill et al. predicted the November 2019 event of the white dwarf LAWD 37. Those
individual measurements were never public. Gaia DR4 publishes them for the first time.

This repository froze, before DR4, the event list (43 events with closest approach before 2020.5),
the fit, the gates and the outcomes: [PREREGISTRATION.md](PREREGISTRATION.md).

## Layout

- `events.py`: the frozen event list, headline set and validation set (Proxima Cen, LAWD 37).
- `lensfit.py`: standard-library fit of the background star's along-scan positions with a
  lensing term; Gaia's parallax geometry; mass and verdict rules.
- `lens_job.py`: release-day job (DR2 -> DR4 matching, DataLink epoch astrometry, fits, outputs).
- `rehearsal.py`: pre-DR4 check of the geometry and a phantom-lens null on ESA's DR4 pre-release.
- `test_lensfit.py`: unit tests, `python -m unittest -v test_lensfit`.

The DR4 I/O (`io_gaia.py`) and the per-transit normal points (`gaiaorb/epochs.py`) are reused
unchanged from [gj3512-giants-dr4](https://github.com/munimtechnologies/gj3512-giants-dr4); the job
copies them in at run time.

Run: `pip install -r requirements.txt` and copy `io_gaia.py` and `gaiaorb/` from gj3512-giants-dr4,
then `python lens_job.py --data-dir data --out-dir out`.

## Credit

Predictions: Klueter et al. (2018, A&A 620, A175) and Klueter, Bastian & Wambsganss (2020, A&A;
arXiv:1911.02584); McGill et al. (2018, arXiv:1804.07049; 2023, arXiv:2206.01814). Critique of
false predicted events: McGill, Everall, Boubert & Smith (2020, arXiv:2006.13958). Gaia data:
ESA/Gaia/DPAC.
