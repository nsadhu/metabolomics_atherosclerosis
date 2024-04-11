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
eqtl <- fread(".../coloc/GTEx_Analysis_v8_eQTL_all_associations_Whole_Blood.allpairs_chr11.txt",
                   header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

eqtl %>% 
  separate(variant_id, c("Chr", "bp", "ref_A2", "alt_A1", NA)) -> eqtl
  #The effect sizes of eQTLs are defined as the effect of the alternative allele (ALT) relative to the reference (REF) allele in the human genome reference. 
  #In other words, the eQTL effect allele is the ALT allele, not the minor allele.

probe = 110371275 #hg38

window_bp = c( (probe - (500 *1000)) : (probe + (500 *1000)) )  # window around peak +/- 1mb

#for 1mb region around peak, merge data 
eqtl %>% 
  mutate(Chr = as.numeric(str_sub(Chr, -2)), 
         bp = as.numeric(bp)) %>% 
  filter(bp %in% window_bp) %>% 
  inner_join(metab, by = c("Chr", "bp")) -> eqtl_metab

#check if alleles are present
eqtl_metab %>% 
  mutate(A1.present = ifelse(A1 == alt_A1 | A1 == ref_A2, 1, 0), 
         A2.present = ifelse(A2 == alt_A1 | A2 == ref_A2, 1, 0)) %>% 
  filter(A1.present == 1 & A2.present == 1) -> eqtl_metab 
 
#all beta values to be aligned to A1 allele
table(eqtl_metab$alt_A1 == eqtl_metab$A1)
table(eqtl_metab$ref_A2 == eqtl_metab$A1)

eqtl_metab %>%
  mutate(correct_slope = ifelse(alt_A1 == A1, slope, -slope), 
         correct_maf = ifelse(alt_A1 == A1, maf, 1-maf)) -> eqtl_metab

rm(eqtl)
rm(metab)

#get all genes in the region
genelist <- unique(eqtl_metab$gene_id)

#extract and format eqtl dataset for coloc
eqtl_metab %>% 
  filter(gene_id == genelist[str_detect(genelist, "ENSG00000137714")]) %>% #change gene Ensembl ID as needed
  mutate(varslope = (slope_se)^2) %>% 
  dplyr::select(correct_slope, varslope, SNP, bp, maf) %>% 
  dplyr::rename(beta = correct_slope,
                varbeta = varslope,
                snp = SNP, 
                position = bp, 
                MAF = maf) %>% 
  distinct() -> d1_eqtl

d1_eqtl <- as.list(d1_eqtl)
d1_eqtl$type = "quant"
d1_eqtl$N = 670

check_dataset(d1_eqtl)

#extract and format metabolite dataset for coloc
eqtl_metab %>% 
  filter(gene_id == genelist[str_detect(genelist, "ENSG00000137714")]) %>% #change gene Ensembl ID as needed
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
my.res <- coloc.abf(dataset1=d1_eqtl, dataset2=d2_metab)
my.res

