
##############################
## TwoSampleMR in Asian popl.
#############################


# remove environment variables
rm(list = ls())

# load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(TwoSampleMR)


#####################################################
# MR with metabolite as exposure and CAD as outcome
#####################################################

#read cad summary stats
cad <- fread(".../BBJCAD_2020_sumstats.txt",
                  header = TRUE, stringsAsFactors = FALSE, data.table = FALSE) 

#read metabolite summary stats
meta <- fread(".../helios_metabolite_gwas.txt",
             header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

#select instruments (variants) associated with metabolite (exposure)
meta %>% 
  filter(p < 5e-08) %>% 
  mutate(phenotype = "metabolite") -> meta
  
#format exposure data
exposure_data <- tryCatch({format_data(meta, 
                                       type = "exposure", 
                                       header = TRUE, 
                                       snp_col = "SNP", 
                                       beta_col = "b", 
                                       se_col = "se", 
                                       pval_col = "p",
                                       eaf_col = "Freq",
                                       effect_allele_col = "A1", 
                                       other_allele_col = "A2", 
                                       samplesize_col = "N", 
                                       phenotype_col = "phenotype")}, 
                          error = function(error) { 
                            message(error) 
                            return(NA)
                          }
                         )
 
#retrieve instruments from outcome cad data 
cad %>% 
  filter(SNP %in% exposure_data$SNP) -> cad

#remove instruments if independently associated with outcome cad data 
cad %>% filter(p> 5e-08) -> cad
  
#format outcome data
outcome_data <- tryCatch({format_data(cad,
                                      type = "outcome", 
                                      header = TRUE, 
                                      samplesize_col = "N",
                                      snp_col = "SNP", 
                                      beta_col = "b", 
                                      se_col = "se", 
                                      pval_col = "p",
                                      eaf_col = "Freq",
                                      effect_allele_col = "A1", 
                                      other_allele_col = "A2")}, 
                         error = function(error) {
                           message(error) 
                           return(NA)
                         }
                        )

#harmonize instruments
dat <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)

#clump instruments
dat <- clump_data(dat, clump_r2 = 0.05, clump_p1 = 5e-08, clump_kb = 1000, pop = "EAS")

#perform MR
set.seed(1234)    
mr_res <- mr(dat)
mr_res


#sensitivity tests
steiger_filtering(dat) %>%  
  mutate(f.stat = ((samplesize.exposure-1-1)/1) * (rsq.exposure/(1-rsq.exposure)))

Isq(dat$beta.exposure, dat$se.exposure)

mr_pleiotropy_test(dat)

run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)

mr_egger_regression(b_exp = dat$beta.exposure, b_out = dat$beta.outcome, se_exp = dat$se.exposure, se_out = dat$se.outcome)

mr_scatter_plot(mr_results = mr_res, dat = dat)

plotres <- mr_leaveoneout(dat = dat)
mr_leaveoneout_plot(plotres)

#####################################################
# MR with metabolite as exposure and CAD as outcome
#####################################################

#read cad summary stats
cad <- fread(".../BBJCAD_2020_sumstats.txt",
             header = TRUE, stringsAsFactors = FALSE, data.table = FALSE) 

#read metabolite summary stats
meta <- fread(".../helios_metabolite_gwas.txt",
              header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

#select instruments (variants) associated with cad (exposure)
cad %>% 
  filter(p < 5e-08) %>% 
  mutate(phenotype = "CAD") -> cad

#format exposure data
exposure_data <- tryCatch({format_data(cad, 
                                       type = "exposure", 
                                       header = TRUE, 
                                       snp_col = "SNP", 
                                       beta_col = "b", 
                                       se_col = "se", 
                                       pval_col = "p",
                                       eaf_col = "Freq",
                                       effect_allele_col = "A1", 
                                       other_allele_col = "A2", 
                                       samplesize_col = "N", 
                                       phenotype_col = "phenotype")}, 
                          error = function(error) { 
                            message(error) 
                            return(NA)
                          }
)


#retrieve instruments from outcome metabolite data 
meta %>% 
  filter(SNP %in% exposure_data$SNP) -> meta

#remove instruments if independently associated with outcome metabolite data 
meta %>% filter(p> 5e-08) -> meta

#format outcome data
outcome_data <- tryCatch({format_data(meta,
                                      type = "outcome", 
                                      header = TRUE, 
                                      samplesize_col = "N",
                                      snp_col = "SNP", 
                                      beta_col = "b", 
                                      se_col = "se", 
                                      pval_col = "p",
                                      eaf_col = "Freq",
                                      effect_allele_col = "A1", 
                                      other_allele_col = "A2")}, 
                         error = function(error) {
                           message(error) 
                           return(NA)
                         }
)


#harmonize instruments
dat <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)

#clump instruments
dat <- clump_data(dat, clump_r2 = 0.05, clump_p1 = 5e-08, clump_kb = 1000, pop = "EAS")

#perform MR
mr_res <- mr(dat)
mr_res

#sensitivity tests
steiger_filtering(dat) %>%  
  mutate(f.stat = ((samplesize.exposure-1-1)/1) * (rsq.exposure/(1-rsq.exposure))) %>% 
  filter(f.stat < 30)
