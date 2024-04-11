##########################
## Colocalization analysis
##########################

# remove environment variables
rm(list = ls())

# load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(coloc)
library(stringr)

#metabolite summary stats
metab <- fread(".../meta_100006370.res_scale.loco.mlma", 
              sep = "\t",  header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

#pqtl
pqtl <- fread(".../coloc/combined_RS2_chunk1969_chr11_FDX1_P10109_OID30078_v1_Cardiometabolic_II.regenie",
                   header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

pqtl %>% 
  rename(Chr = CHROM, 
         bp = GENPOS) -> pqtl 

probe = 110371275 #hg38

window_bp = c( (probe - (500 *1000)) : (probe + (500 *1000)) )  # window around peak +/- 1mb

#for 1mb region around peak, merge data 
pqtl %>% 
  filter(bp %in% window_bp) %>% 
  inner_join(metab, by = c("Chr", "bp")) -> pqtl_metab

#check if alleles are present
pqtl_metab %>% 
  mutate(A1.present = ifelse(A1 == ALLELE1 | A1 == ALLELE0, 1, 0), 
         A2.present = ifelse(A2 == ALLELE1 | A2 == ALLELE0, 1, 0)) %>% 
  filter(A1.present == 1 & A2.present == 1) -> pqtl_metab 

#all beta values to be aligned to A1 allele
table(pqtl_metab$ALLELE1 == pqtl_metab$A1)
table(pqtl_metab$ALLELE0 == pqtl_metab$A1)

pqtl_metab %>%
  mutate(correct_beta = ifelse(ALLELE1 == A1, BETA, -BETA), 
         correct_freq = ifelse(ALLELE1 == A1, A1FREQ, 1-A1FREQ)) -> pqtl_metab

rm(pqtl)
rm(metab)

#extract and format pqtl dataset for coloc
pqtl_metab %>% 
  mutate(varbeta = SE^2) %>% 
  dplyr::select(correct_beta, varbeta, SNP, bp, correct_freq) %>% 
  dplyr::rename(beta = correct_beta,
                snp = SNP, 
                position = bp, 
                MAF = correct_freq) %>% 
  distinct() -> d1_pqtl

d1_pqtl <- as.list(d1_pqtl)
d1_pqtl$type = "quant"
d1_pqtl$N = 49748

check_dataset(d1_pqtl)

#extract and format metabolite dataset for coloc
pqtl_metab %>% 
  mutate(varbeta = se^2) %>%  
  dplyr::select(b, varbeta, SNP, bp, Freq) %>% 
  dplyr::rename(beta = b, 
                snp = SNP, 
                position = bp, 
                MAF = Freq) %>% 
  distinct() -> d2_metab

d2_metab <- as.list(d2_metab)
d2_metab$type = "quant"
d2_metab$N = 1876

check_dataset(d2_metab)

#colocalization analysis
my.res <- coloc.abf(dataset1=d1_pqtl, dataset2=d2_metab)
my.res

