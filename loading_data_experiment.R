###############################################################################
# Loading data experiment ----
###############################################################################
loading_data_experiment <- function() {
  
  data <- data.frame(
    # No. screened
    data_screened_cage1           = c(100, 148, 1253, 1066, 1036, 1033, 1427, 1000),
    data_screened_cage2           = c(100, 262, 1538, 1221, 1052, 1178, 1333, 1053),
    data_screened_cage3           = c(100, 1406, 1061, 1118, 1434, 1080, 1027, 1032),
    
    # phenotype ----
    data_phenotype_V9_cage1       = c(50, 81, 258, 524, 598, 362, 545, 421),
    data_phenotype_V9_cage2       = c(50, 167, 402, 570, 533, 505, 667, 472),
    data_phenotype_V9_cage3       = c(50, 208, 225, 458, 359, 432, 440, 390),
    
    # phenotype ----
    data_phenotype_LDonly_cage1   = c(0, 15, 147, 211, 129, 135, 257, 126),
    data_phenotype_LDonly_cage2   = c(0, 27, 163, 183, 139, 167, 209, 137),
    data_phenotype_LDonly_cage3   = c(0, 231, 138, 145, 323, 203, 189, 211),
    
    data_phenotype_GFPLonly_cage1 = c(50, 47, 481, 322, 408, 388, 513, 408),
    data_phenotype_GFPLonly_cage2 = c(50, 83, 632, 429, 405, 480, 603, 390),
    data_phenotype_GFPLonly_cage3 = c(50, 525, 401, 492, 352, 330, 337, 339),
    
    data_phenotype_LDGFPL_cage1   = c(50, 86, 625, 533, 499, 510, 657, 466),
    data_phenotype_LDGFPL_cage2   = c(50, 152, 743, 609, 508, 531, 521, 526),
    data_phenotype_LDGFPL_cage3   = c(50, 650, 522, 481, 759, 547, 501, 482),
    
    # No. genotyped ----
    data_genotyped_cage1          = c(50, 50, 50, 50, 50, 50, 50, 50),
    data_genotyped_cage2          = c(50, 50, 50, 50, 50, 50, 50, 50),
    data_genotyped_cage3          = c(50, 50, 50, 50, 50, 50, 50, 50),
    
    # allele frequency ----
    data_allelefreq_L224_cage1    = 50 * c(75, 23.21, 8.93, 12.18, 11.03, 10.94, 9.73, 3.08) / 100,
    data_allelefreq_L224_cage2    = 50 * c(75, 15.9, 28.94, 15.59, 7.13, 6.07, 4.59, 7.04) / 100,
    data_allelefreq_L224_cage3    = 50 * c(75, 23.3, 17.84, 16.28, 7.28, 11.41, 8.96, 3.68) / 100,
    
    data_allelefreq_Q224_cage1    = 50 * c(25, 74.76, 88.09, 84.28, 87.31, 88.38, 89.11, 94.85) / 100,
    data_allelefreq_Q224_cage2    = 50 * c(25, 74.82, 66.4, 76.73, 90.39, 93.33, 94.36, 90.85) / 100,
    data_allelefreq_Q224_cage3    = 50 * c(25, 71.8, 77.53, 76.49, 90.19, 84.37, 89.36, 94.58) / 100,
    
    data_allelefreq_NHEJ_cage1    = 50 * c(0, 2.03, 2.98, 3.54, 1.67, 0.68, 1.16, 2.07) / 100,
    data_allelefreq_NHEJ_cage2    = 50 * c(0, 9.28, 4.66, 7.68, 2.48, 0.6, 1.05, 2.11) / 100,
    data_allelefreq_NHEJ_cage3    = 50 * c(0, 4.9, 4.64, 7.23, 2.53, 4.22, 1.68, 1.74) / 100
  )
  
  return(data)
}
