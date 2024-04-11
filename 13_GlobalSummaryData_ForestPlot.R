
##############################
## Forest Plot for Global data
#############################


# remove environment variables
rm(list = ls())

# load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(forestplot)

##########################
## Metabolites Plot ##
##########################
base_data <- tibble::tibble(study = c("HELIOS", "CLSA", "METSIM"),
                            samplesize = c("1,876", "8,299", "6,136"), 
                            trait = c("3BH5C", "3BH5C", "3BH5C"), 
                            ea = c("T", "T", "T"), 
                            freq = c("0.30", "0.12", "0.12"), 
                            beta = c(-0.48, -0.41, -0.51), 
                            se = c(0.04, 0.02, 0.03), 
                            p = c("8.93E-41*", "6.75E-72*", "5.60E-76*"))


base_data %>% 
  mutate(lower = beta-1.96*se, 
         upper = beta+1.96*se) -> base_data

forestplot(base_data, labeltext = c(study, samplesize, trait, ea, freq, p), 
           align = "l", 
           graph.pos = 7, 
           colgap = unit(1, "cm"), 
           lwd.ci = 2,
           vertices = FALSE,
           title = "rs10488763 (chr11:110373636_T/A)",
           mean = beta, 
           lower = lower, 
           upper = upper, 
           xlab = "beta per allele copy (95% CI)", 
           boxsize = 0.15,
           txt_gp = fpTxtGp(label = gpar(cex = 1.5), xlab = gpar(cex = 1.5), ticks = gpar(cex = 1.5)), 
           xticks = c(-0.75, -0.50, -0.25, 0, 0.25), 
           col = fpColors(box = "darkblue", lines = "royalblue", zero = "darkgray")) %>%
  fp_add_header(study = "GWAS\nCohort", samplesize = "Sample\nSize", trait = "Trait", ea = "Effect\nAllele", freq = "Allele\nFrequency", p = "P-value") 


##########################
## CAD Plot ##
##########################
base_data <- tibble::tibble(study = c("BBJ", "UKBB", "CARDIoGRAM+C4D"),
                            samplesize = c("168,228", "296,525", "184,305"), 
                            trait = c("CAD", "CAD", "CAD"), 
                            ea = c("T", "T", "T"), 
                            freq = c("0.41", "0.13", "0.15"), 
                            beta = c(0.08, 0.04, 0.05), 
                            se = c(0.01, 0.01, 0.01), 
                            p = c("2.03E-15*", "6.20E-04", "1.02E-04"))

base_data %>% 
  mutate(OR = exp(beta), 
         lower = exp(beta-1.96*se), 
         upper = exp(beta+1.96*se)) -> base_data

forestplot(base_data, labeltext = c(study, samplesize, trait, ea, freq, p), 
           align = "l", 
           graph.pos = 7, 
           colgap = unit(1, "cm"), 
           lwd.ci = 2,
           vertices = FALSE,
           title = "rs10488763 (chr11:110373636_T/A)",
           mean = OR, 
           lower = lower, 
           upper = upper, 
           xlab = "OR (95% CI)", 
           zero = 1, 
           boxsize = 0.15,
           txt_gp = fpTxtGp(label = gpar(cex = 1.5), xlab = gpar(cex = 1.5), ticks = gpar(cex = 1.5)), 
           xticks = c(0.9, 1, 1.1, 1.2),
           col = fpColors(box = "darkblue", lines = "royalblue", zero = "darkgray")) %>%
  fp_add_header(study = "GWAS\nCohort", samplesize = "Sample\nSize", trait = "Trait", ea = "Effect\nAllele", freq = "Allele\nFrequency", p = "P-value") 

