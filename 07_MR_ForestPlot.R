
##############################
## Forest Plot for MR results
#############################


# remove environment variables
rm(list = ls())

# load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(forestplot)
library(metafor)


base_data <- tibble::tibble(method = c("", "GSMR", "Inverse Variance Weighted", "", "", "Inverse Variance Weighted", "Weighted-median", "Weighted-mode", "MR-Egger"),
                            snp = c(NA, 13, 2, NA, NA, 12, 12, 12, 12), 
                            bxy = c(NA, -0.116, -0.151, NA, NA, -0.024, -0.023, -0.024, -0.029), 
                            se = c(NA, 0.014, 0.035, NA, NA, 0.006, 0.006, 0.006, 0.009), 
                            p = c("", "1.5E-16", "1.9E-05", "", "", "3.0E-05", "8.6E-05", "2.1E-03", "1.2E-02"))

base_data %>% 
  mutate(OR = exp(bxy), 
         lower = exp(bxy-1.96*se), 
         upper = exp(bxy+1.96*se)) -> base_data

forestplot(base_data, labeltext = c(method, snp, p), 
           align = "l", 
           graph.pos = 4, 
           colgap = unit(0.3, "cm"), 
           lwd.ci = 2,
           vertices = FALSE,
           mean = OR, 
           lower = lower, 
           upper = upper, 
           zero = 1, 
           xlab = "OR per 1SD increase in 3BH5C (95% CI)", 
           boxsize = 0.20,
           txt_gp = fpTxtGp(label = gpar(cex = 1.5), xlab = gpar(cex = 1.5), ticks = gpar(cex = 1.5)), 
           xticks = c(0.7, 0.8, 0.9, 1, 1.1),
           col = fpColors(box="darkblue",line="royalblue", zero="darkgray")) %>% 
  fp_add_header(method = "Method",snp = "N SNPs", p = "P-value") %>% 
  fp_set_zebra_style("#f2f2f2")

## test for heterogeneity of IVW estimates from discovery (asian) and replication (european)
rma(yi = c(-0.186, -0.024), sei = c(0.039, 0.006), method = "FE") #fixed effect
rma(yi = c(-0.186, -0.024), sei = c(0.039, 0.006)) # random effect
