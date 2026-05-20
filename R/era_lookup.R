# era_lookup.R
# Treatment-era fallback when an explicit protocol arm is missing.
# Cohort restricted to dx_year >= 2010 (see config/features.yaml).
# Eras chosen around major shifts in pediatric/AYA/adult B-ALL therapy
# within the 2010-2025 window:
#   - era_2010_2014 : MRD-driven escalation/de-escalation entrenched;
#                     dex + HD-MTX standardized; frontline imatinib in Ph+
#   - era_2015_2019 : blinatumomab in salvage (AALL1331), Ph-like targeted
#                     trials (AALL1521), bortezomib in T-ALL (AALL1231)
#   - era_2020_2025 : frontline immunotherapy (AALL1731 SR + blina,
#                     AALL1732 HR + inotuzumab), ALLTogether1, CAR-T
#                     in salvage routine, E1910 MRD-neg adult B-ALL

PROTOCOLS_BY_ERA <- list(
  era_2010_2014 = c(
    # Pediatric B-ALL frontline
    "COG-AALL0232 (late)", "COG-AALL0331 (late)",
    "COG-AALL0932",                                    # SR B-ALL (2010-2018)
    "COG-AALL1131",                                    # HR B-ALL (2012-2019)
    "COG-AALL0434",                                    # T-ALL, nelarabine
    "COG-AALL0622",                                    # Ph+, dasatinib
    "COG-AALL1231",                                    # T-ALL, bortezomib (2014-)
    "SJCRH-TotalXVI",                                  # 2007-2017
    "AIEOP-BFM-ALL-2009",                              # 2010-2017
    "NOPHO-ALL2008",                                   # 2008-2018
    "UKALL-2011",                                      # 2012-2019
    "DCOG-ALL11",                                      # 2012-
    "EsPhALL",                                         # Ph+ + imatinib
    "INTERFANT-06",                                    # 2006-2016
    # AYA / adult
    "CALGB-10403",                                     # AYA pediatric-inspired (-2012)
    "DFCI-11-001",
    "GRAALL-2005 / 2014 transition",
    "GMALL-08/2013",
    "HyperCVAD (+/- rituximab)"
  ),
  era_2015_2019 = c(
    "COG-AALL1131",                                    # ongoing
    "COG-AALL1231",                                    # ongoing
    "COG-AALL1331",                                    # relapsed B-ALL, blinatumomab
    "COG-AALL1521",                                    # Ph-like + TKI
    "COG-AALL1631",                                    # Ph+ + imatinib
    "SJCRH-TotalXVI / Total-XVII transition",
    "AIEOP-BFM-ALL-2009 / 2017 transition",
    "UKALL-2011 (late)",
    "DCOG-ALL11",
    "ALLTogether1 pilot",
    # AYA / adult
    "GRAALL-2014",
    "GMALL-08/2013",
    "ECOG-E1910 (start)",                              # MRD-neg adult B-ALL + blina
    "HyperCVAD + blinatumomab/inotuzumab (early)"
  ),
  era_2020_2025 = c(
    "COG-AALL1731",                                    # SR + blinatumomab (2018-)
    "COG-AALL1732",                                    # HR + inotuzumab (2019-)
    "COG-AALL1631 (ongoing)",
    "COG-AALL2031", "COG-AALL2121",                    # newer ongoing
    "SJCRH-TotalXVII",                                 # 2017-
    "AIEOP-BFM-ALL-2017",
    "ALLTogether1",                                    # pan-European, 2020-
    "INTERFANT-21",
    "UKALL-2019",
    "EsPhALL / COG AALL1631 hybrid",
    # AYA / adult
    "ECOG-E1910",
    "GIMEMA-LAL2317",
    "MDACC Mini-HCVD + inotuzumab",
    "HyperCVAD + blinatumomab/inotuzumab",
    "CAR-T salvage (tisagenlecleucel / brexucabtagene)"
  )
)

# Map dx_year -> era label. Returns NA for years outside the configured
# cohort window so 03_qc_and_merge.R can drop those rows.
derive_treatment_era <- function(dx_year,
                                 min_year = 2010,
                                 max_year = 2025) {
  era <- cut(dx_year,
             breaks = c(-Inf, 2014, 2019, Inf),
             labels = c("era_2010_2014", "era_2015_2019", "era_2020_2025"),
             right  = TRUE)
  era[dx_year < min_year | dx_year > max_year] <- NA
  era
}
