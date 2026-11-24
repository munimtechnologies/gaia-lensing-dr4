"""Frozen event list (2026-11-24, before Gaia DR4).

Part of gaia-lensing-dr4 (MV2, Munim, Inc.). Apache-2.0. Made with AI agents working in MV2's name.

Source: Klueter, Bastian & Wambsganss 2020, A&A, "Expectations on the mass
determination using astrometric microlensing by Gaia" (arXiv:1911.02584), Table 1 (closest
approach before 2019.5; with the 5-year-mission error sM5) and Table 2 (after 2019.5), all
events with closest approach T_CA before 2020.5 (DR4 covers 2014.57-2020.05; later events
contribute their leading tail). Columns: event number in the paper, lens name, Gaia DR2
source_id of lens and background source, T_CA (Julian year), the mass the paper assumed
(M_in, Msun), its expected 1-sigma mass error for the 10-year (sM10) and 5-year (sM5)
missions (None where the paper gives none or >= 100%), and whether DR2 had only a
2-parameter solution for the source (asterisk in the paper).
Transcribed from the paper's HTML; the release-day job checks every pair (the DR2 tracks must
come within 3 arcsec of each other within 1 year of T_CA), so a transcription slip shows up
as ID_MISMATCH instead of a wrong answer.
"""

EVENTS = [
    # (num, lens, lens_dr2, source_dr2, t_ca, m_in, sm10, sm5, src_2p)
    (1, "HD 22399", 488099359330834432, 488099363630877824, 2013.711, 1.3, 0.60, None, True),
    (2, "HD 177758", 4198685678421509376, 4198685678400752128, 2013.812, 1.1, 0.41, None, True),
    (3, "L 820-19", 5736464668224470400, 5736464668223622784, 2014.419, 0.28, 0.058, 0.11, True),
    (4, "", 2081388160068434048, 2081388160059813120, 2014.526, 0.82, 0.39, None, True),
    (5, "", 478978296199510912, 478978296204261248, 2014.692, 0.65, 0.24, 0.49, False),
    (6, "G 123-61B", 1543076475514008064, 1543076471216523008, 2014.763, 0.26, 0.099, 0.18, True),
    (7, "L 601-78", 5600272625752039296, 5600272629670698880, 2014.783, 0.21, 0.041, 0.068, True),
    (8, "UCAC3 27-74415", 6368299918479525632, 6368299918477801728, 2015.284, 0.36, 0.047, 0.11, True),
    (9, "BD+00 5017", 2646280705713202816, 2646280710008284416, 2015.471, 0.58, 0.18, 0.44, True),
    (10, "G 123-61A", 1543076475509704192, 1543076471216523008, 2016.311, 0.32, 0.079, 0.14, True),
    (11, "EC 19249-7343", 6415630939116638464, 6415630939119055872, 2016.650, 0.26, 0.099, None, True),
    (12, "", 5312099874809857024, 5312099870497937152, 2016.731, 0.07, 0.0088, 0.0098, False),
    (13, "PM J08503-5848", 5302618648583292800, 5302618648591015808, 2017.204, 0.65, 0.17, 0.29, False),
    (14, "", 5334619419176460928, 5334619414818244992, 2017.258, 0.65, 0.23, 0.43, False),
    (15, "Proxima Cen", 5853498713160606720, 5853498708818460032, 2017.392, 0.12, 0.034, 0.082, False),
    (16, "Innes' star", 5339892367683264384, 5339892367683265408, 2017.693, 0.33, 0.16, 0.22, False),
    (17, "L 51-47", 4687511776265158400, 4687511780573305984, 2018.069, 0.28, 0.035, 0.046, False),
    (18, "", 4970215770740383616, 4970215770743066240, 2018.098, 0.65, 0.076, 0.21, False),
    (19, "BD-06 855", 3202470247468181632, 3202470247468181760, 2018.106, 0.8, 0.38, 0.72, True),
    (20, "G 217-32", 429297924157113856, 429297859741477888, 2018.134, 0.23, 0.018, 0.040, True),
    (21, "", 5865259639247544448, 5865259639247544064, 2018.142, 0.8, 0.053, 0.12, True),
    (22, "HD 146868", 1625058605098521600, 1625058605097111168, 2018.183, 0.92, 0.060, 0.14, True),
    (23, "Ross 733", 4516199240734836608, 4516199313402714368, 2018.282, 0.43, 0.067, 0.13, False),
    (24, "LP 859-51", 6213824650812054528, 6213824650808938880, 2018.359, 0.43, 0.19, None, True),
    (25, "L 230-188", 4780100658292046592, 4780100653995447552, 2018.450, 0.15, 0.073, 0.11, False),
    (26, "HD 149192", 5930568598406530048, 5930568568425533440, 2018.718, 0.67, 0.043, 0.12, True),
    (27, "HD 85228", 5309386791195469824, 5309386795502307968, 2018.751, 0.82, 0.043, 0.65, True),
    (28, "HD 155918", 5801950515627094400, 5801950515623081728, 2018.773, 1.0, 0.13, 0.35, True),
    (29, "HD 77006", 1015799283499485440, 1015799283498355584, 2018.796, 1.0, 0.42, None, True),
    (30, "Proxima Cen", 5853498713160606720, 5853498713181091840, 2018.819, 0.12, 0.020, 0.052, False),
    (31, "HD 2404", 2315857227976341504, 2315857227975556736, 2019.045, 0.84, 0.19, None, True),
    (32, "LP 350-66", 2790883634570755968, 2790883634570196608, 2019.299, 0.32, 0.15, None, True),
    (33, "HD 110833", 1568219729458240128, 1568219729456499584, 2019.360, 0.78, 0.056, None, True),
    (34, "L 702-43", 4116840541184279296, 4116840536790875904, 2019.673, 0.22, 0.085, None, False),
    (35, "G 251-35", 1136512191212093440, 1136512191210614272, 2019.856, 0.64, 0.26, None, True),
    (36, "G 245-47A", 546488928621555328, 546488928619704320, 2019.858, 0.27, 0.075, None, False),
    (37, "LAWD 37", 5332606522595645952, 5332606350796955904, 2019.865, 0.65, 0.13, None, False),
    (38, "LSPM J2129+4720", 1978296747258230912, 1978296747268742784, 2019.960, 0.34, 0.13, None, True),
    (39, "", 5788178166117584640, 5788178170416392192, 2020.126, 0.54, 0.12, None, True),
    (40, "G 16-29", 4451575895403432064, 4451575895400387968, 2020.214, 0.51, 0.19, None, True),
    (41, "HD 66553", 654826970401335296, 654826970399770368, 2020.309, 0.87, 0.11, None, True),
    (42, "HD 120065", 1251328585567327744, 1251328589861875840, 2020.344, 1.2, 0.54, None, True),
    (43, "HD 78663", 3842095911266162432, 3842095915561269888, 2020.356, 0.8, 0.28, None, True),
]

FIELDS = ("num", "lens", "lens_dr2", "source_dr2", "t_ca", "m_in", "sm10", "sm5", "src_2p")

# Headline set: events whose 5-year expected relative error sM5 / M_in is at most 20%, plus
# LAWD 37 (#37), whose November 2019 peak falls inside DR4 and whose 10-year error is 19%.
HEADLINE = (12, 17, 20, 21, 22, 26, 37)

# Validation set: lenses with independently known masses (Msun, 1-sigma).
#  Proxima Cen: 0.1221 +- 0.0022 (Kervella, Thevenin & Lovis 2017, A&A, from the
#    Mann et al. 2015 mass-luminosity relation).
#  LAWD 37: 0.56 +- 0.08 (McGill et al. 2023, MNRAS; HST astrometric microlensing of
#    its November 2019 event, arXiv:2206.01814).
VALIDATION = {15: (0.1221, 0.0022), 30: (0.1221, 0.0022), 37: (0.56, 0.08)}


def as_dicts():
    return [dict(zip(FIELDS, e)) for e in EVENTS]
