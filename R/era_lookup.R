# era_lookup.R
# Treatment-era fallback when an explicit protocol arm is missing.
# Eras chosen around major shifts in pediatric/AYA B-ALL therapy:
#   - Early MRD-stratification era (2000-2005)
#   - Dex/HD-MTX optimization + MRD-driven intensification (2006-2011)
#   - TKI in frontline Ph+, expanded MRD-driven arms (2012-2017)
#   - Frontline immunotherapy (blinatumomab/inotuzumab), Ph-like TKI trials,
#     CAR-T in salvage (2018-2025)
#
# Reference protocols by era (NOT exhaustive). Useful when populating
# protocol_arm from chart review; era is auto-derived from dx_year otherwise.

PROTOCOLS_BY_ERA <- list(
  era_2000_2005 = c(
    # Pediatric B-ALL (frontline)
    "CCG-1991", "CCG-1961", "CCG-1952",                # SR / HR / infant (pre-COG merge)
    "POG-9904", "POG-9905", "POG-9906",
    "COG-AALL0331",                                    # SR B-ALL (2005-)
    "COG-AALL0232",                                    # HR B-ALL (2004-)
    "COG-AALL0031",                                    # Ph+ (imatinib)
    "SJCRH-TotalXV",                                   # 2000-2007
    "ALL-BFM-2000", "AIEOP-BFM-ALL-2000",
    "NOPHO-ALL2000",
    "UKALL-2003",                                      # 2003-2011
    "DCOG-ALL10",                                      # 2004-2013
    "INTERFANT-99",
    # AYA / adult
    "CALGB-19802", "CALGB-10102",
    "GMALL-07/2003", "GRAALL-2003",
    "UKALL-XII / ECOG-E2993",
    "HyperCVAD"
  ),
  era_2006_2011 = c(
    "COG-AALL0232", "COG-AALL0331",
    "COG-AALL0434",                                    # T-ALL, nelarabine
    "COG-AALL0622",                                    # Ph+, dasatinib
    "COG-AALL0932",                                    # SR B-ALL (2010-)
    "SJCRH-TotalXVI",                                  # 2007-2017
    "AIEOP-BFM-ALL-2009",                              # 2010-2017
    "NOPHO-ALL2008",                                   # 2008-2018
    "UKALL-2003 (late)",
    "DCOG-ALL10/11",
    "INTERFANT-06",
    "EsPhALL",                                         # Ph+, imatinib
    "CALGB-10403",                                     # AYA pediatric-inspired (2007-2012)
    "GRAALL-2005", "GMALL-08/2013-pilot",
    "DFCI-05-001", "DFCI-11-001"
  ),
  era_2012_2017 = c(
    "COG-AALL1131",                                    # HR B-ALL, clofarabine, MRD-driven
    "COG-AALL1231",                                    # T-ALL, bortezomib
    "COG-AALL1331",                                    # Relapsed B-ALL, blinatumomab
    "COG-AALL1521",                                    # Ph-like + TKI
    "COG-AALL1631",                                    # Ph+, imatinib + chemo
    "SJCRH-TotalXVI",
    "AIEOP-BFM-ALL-2009 / 2017 transition",
    "NOPHO-ALL2008 (late)",
    "UKALL-2011",                                      # 2012-2019
    "DCOG-ALL11",
    "DFCI-11-001",
    "GRAALL-2014",
    "GMALL-08/2013"
  ),
  era_2018_2025 = c(
    "COG-AALL1731",                                    # SR B-ALL + blinatumomab
    "COG-AALL1732",                                    # HR B-ALL + inotuzumab
    "COG-AALL1631 (ongoing)",
    "COG-AALL2031", "COG-AALL2121",                    # newer ongoing
    "SJCRH-TotalXVII",                                 # 2017-
    "AIEOP-BFM-ALL-2017",
    "ALLTogether1",                                    # pan-European, 2020-
    "INTERFANT-21",
    "UKALL-2019" ,
    "EsPhALL/COG AALL1631 hybrid",
    # AYA / adult
    "ECOG-E1910 (blinatumomab in MRD-neg adult B-ALL)",
    "GIMEMA-LAL2317", "GRAALL-2014 / Quest",
    "GMALL-08/2013 + inotuzumab", "HyperCVAD + blinatumomab/inotuzumab",
    "MDACC-Mini-HCVD + inotuzumab (older adults)"
  )
)

# Map dx_year -> era label. Vectorized.
derive_treatment_era <- function(dx_year) {
  cut(dx_year,
      breaks = c(-Inf, 2005, 2011, 2017, Inf),
      labels = c("era_2000_2005", "era_2006_2011",
                 "era_2012_2017", "era_2018_2025"),
      right  = TRUE)
}
