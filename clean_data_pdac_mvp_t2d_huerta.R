# cleaning pdac_GCST90475570_mvp.tsv which is Verma MVP data
# with 1880 cases and 500k controls
# https://www.ebi.ac.uk/gwas/studies/GCST90475570

# Clean Huerta-Chagoya Nat Gen 2024 T2D which is a meta-analysis of 51k cases and 370k controls with
# 12.2% cases of non-European ancestry. It is UKB imputed with topmed, MGH BIoank, GERA cohort,
# and AllofUs.

# The raw datasets are available at the GWAS catalog.

library(data.table)
library(tidyverse)

# 19.7 million variants
setwd("/rsrch3/home/biostatistics/rsun3/pdac_pleio/data")
pdac <- fread("pdac_GCST90475570_mvp.tsv", select=c("chromosome", "base_pair_location",
                                               "effect_allele", "other_allele", "odds_ratio", "p_value", 
                                               "rsid", "effect_allele_frequency")) %>%
  mutate(name = paste0(chromosome, "_",  base_pair_location)) %>%
  mutate(Zpdac = qnorm(1 - p_value / 2)) %>%
  mutate(Zpdac= ifelse(odds_ratio > 1, Zpdac, -1 * Zpdac)) %>%
  mutate(pval_pdac = p_value) %>%
  select(-p_value)

# grch38
# pdac[7673538, ] 

# 146M million variants
t2d <- fread("diabetes_GCST90444202.tsv", select=c("chromosome", "base_pair_location", "beta", "standard_error", "p_value",
                                                   "effect_allele", "other_allele")) %>%
  mutate(name = paste0(chromosome, "_", base_pair_location)) 
colnames(t2d)[6:7] <- c("EA", "NEA")


# get the common positions
# 18.5M variants
commonPDACdiab <- merge(pdac, t2d, by="name") %>%
  mutate(match = ifelse(EA == effect_allele & NEA == other_allele, 1, 0)) %>%
  mutate(flip = ifelse(EA == other_allele & NEA == effect_allele, 1, 0)) %>%
  mutate(Zdiab = beta / standard_error) %>%
  mutate(Zdiab = ifelse(flip == 1, -1 * Zdiab, Zdiab)) %>%
  filter(match == 1 | flip == 1)

# save
write.table(commonPDACdiab, "final_pdac_mvp_diabetes_huerta.txt", append=F, quote=F, row.names=F, col.names=T, sep='\t')



