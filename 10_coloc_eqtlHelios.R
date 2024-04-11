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

#eqtl
eqtl <- fread(".../coloc/sig_cis_HELIOS_compiled_final.txt",
                   header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

eqtl %>% 
  filter(str_starts(SNP, "11:")) -> eqtl

eqtl %>% 
  mutate(snp_id = SNP) %>% 
  separate(SNP, c("Chr", "bp", NA, NA)) -> eqtl #get alleles later

probe = 110371275 #hg38

window_bp = c( (probe - (500 *1000)) : (probe + (500 *1000)) )  # window around peak +/- 1mb

#for 1mb region around peak, merge data 
eqtl %>% 
  mutate(Chr = as.numeric(Chr), 
         bp = as.numeric(bp)) %>% 
  filter(bp %in% window_bp) %>% 
  inner_join(metab, by = c("Chr", "bp")) -> eqtl_metab

#get alleles and merge
eqtl_allele <- fread(".../coloc/helios_eqtl/sg10k_maf0.001_hwe1e-03_Rsq0.30.bim",
                     header = FALSE, stringsAsFactors = FALSE, data.table = FALSE)

colnames(eqtl_allele) <- c("Chr", "SNP", "V3", "bp", "eqtl_A1", "eqtl_A2")

eqtl_metab %>% 
  left_join(eqtl_allele, by = "SNP") %>% 
  drop_na() -> eqtl_metab

#check if alleles are present
eqtl_metab %>% 
  mutate(A1.present = ifelse(A1 == eqtl_A1 | A1 == eqtl_A2, 1, 0), 
         A2.present = ifelse(A2 == eqtl_A1 | A2 == eqtl_A2, 1, 0)) %>% 
  filter(A1.present == 1 & A2.present == 1) -> eqtl_metab 

#all beta values to be aligned to A1 allele
table(eqtl_metab$eqtl_A1 == eqtl_metab$A1)
table(eqtl_metab$eqtl_A2 == eqtl_metab$A1)

eqtl_metab %>%
  mutate(correct_beta = ifelse(eqtl_A1 == A1, beta, -beta)) -> eqtl_metab

rm(eqtl)
rm(metab)

#get all genes in the region
genelist <- unique(eqtl_metab$gene)

#extract and format eqtl dataset for coloc
eqtl_metab %>% 
  filter(gene == genelist[str_detect(genelist, "ENSG00000137714")]) %>% #change gene Ensembl ID as needed
  mutate(varbeta = (beta/`t-stat`)^2) %>% 
  dplyr::select(correct_beta, varbeta, SNP, bp.x, Freq) %>% 
  dplyr::rename(beta = correct_beta,
                snp = SNP, 
                position = bp.x, 
                MAF = Freq) %>% 
  distinct() -> d1_eqtl

d1_eqtl <- as.list(d1_eqtl)
d1_eqtl$type = "quant"
d1_eqtl$N = 1228

check_dataset(d1_eqtl)

#extract and format metabolite dataset for coloc
eqtl_metab %>% 
  filter(gene == genelist[str_detect(genelist, "ENSG00000137714")]) %>% #change gene Ensembl ID as needed
  mutate(varbeta = se^2) %>%  
  dplyr::select(b, varbeta, SNP, bp.x, Freq) %>% 
  dplyr::rename(beta = b, 
                snp = SNP, 
                position = bp.x, 
                MAF = Freq) %>% 
  distinct() -> d2_metab

d2_metab <- as.list(d2_metab)
d2_metab$type = "quant"
d2_metab$N = 1876

check_dataset(d2_metab)

#colocalization analysis
my.res <- coloc.abf(dataset1=d1_eqtl, dataset2=d2_metab)
my.res

