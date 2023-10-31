
##############################
## Linear regression Model 1
#############################


# remove environment variables
rm(list = ls())

# load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)

working = "C:/Users/Nilanjana/OneDrive - Nanyang Technological University/Metabolomics/NS_Helios_Metabolon/"

# read metabolon data
df <- fread(paste0(working, "Batch-normImpDataAll_Jan2023_missing25_zerovar_1073metabolites_235repVis2Excl_9pc1OutliersExcl_log_std.csv"), 
            sep = ",",  header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)
df[1:5, 1:5]

# read pheno data
pheno <- fread(paste0(working, "HELIOS_Metabolon_demo_carotidFinal11k_logtransform_boxBatch_riskfactors_8192_Nov2022.txt"), 
                 sep = "\t",  header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

head(pheno)

# subset and factorize phenotypes 
##ETHNICITY - remove OTHERS
pheno %>%
  filter(Ethnicity != "O") -> df_pheno

df_pheno$Ethnicity = as.factor(df_pheno$Ethnicity)
df_pheno$Sex = as.factor(df_pheno$Sex)
df_pheno$SSID = as.factor(df_pheno$SSID)

# drop na from outcome to analyze
df_pheno %>% 
  drop_na(ln_mean_cIMT, Age, Sex, Ethnicity) -> df_pheno

# subset and sort samples in same order
df %>% 
  filter(PARENT_SAMPLE_NAME %in% df_pheno$PARENT_SAMPLE_NAME) %>% 
  arrange(PARENT_SAMPLE_NAME) -> df

table(df$PARENT_SAMPLE_NAME == df_pheno$PARENT_SAMPLE_NAME)

df = df[, -1] # all columns except sample names
df[1:5, 1:5]

# remove metabolites with only one value i.e. column variance = 0
df %>% 
  dplyr::select_if(~var(., na.rm = TRUE) != 0) -> df

# save metabolite names
CHEM_ID <- colnames(df)
colnames(df) = NULL
df <- as.matrix(df)

# run linear regression

lm.formula <- as.formula('df_pheno$ln_mean_cIMT ~ df[, i] +
                         df_pheno$Age + df_pheno$Sex + df_pheno$Ethnicity + df_pheno$SSID')

res = data.frame(matrix(nrow = ncol(df), ncol = 4))
colnames(res) =c('Est','SE', 'z_value','P')

system.time(
  for(i in 1:ncol(df)) {
    tryCatch({reg.out <- summary(glm(lm.formula, family = gaussian(), na.action = na.omit))}, error = function(error) {return(NA)})
    if(!exists("reg.out")) {
      res[i,] = rep(NA,4)
    } else {
      res[i,] = tryCatch({reg.out$coefficients[2,]},error = function(error) {return(rep(NA,4))})
      rm(reg.out)
    }
    if(i%%100 == 0) {
      print(i)
    }
  }
)

final_out <- cbind(CHEM_ID, res, stringsAsFactors = FALSE)
head(final_out)
final_out %>% 
  arrange(P)

# get chemical annotation file
df_chem <- fread("C:/Users/Nilanjana/OneDrive - Nanyang Technological University/Metabolomics/HELIOS data-Metabolome/HELIOS-Metabolon_281222/Chemical Annotation All.csv", 
                 header = TRUE, stringsAsFactors = FALSE, fill = TRUE, data.table = FALSE)
head(df_chem)

library(stringr)
df_chem %>%
  mutate_if(is.character, str_trim) -> df_chem

df_chem$CHEM_ID = as.character(df_chem$CHEM_ID)

# merge chemical annotation information with regression results
final_out %>% 
  mutate(adjP_FDR = p.adjust(P, method = "BH"), 
         adjP_Bonf = p.adjust(P, method = "bonferroni")) %>%
  left_join(df_chem, by = "CHEM_ID") -> final_out_annot

head(final_out_annot)

final_out_annot %>% 
  dplyr::select(SUPER_PATHWAY, SUB_PATHWAY, CHEMICAL_NAME, Est, SE, P, adjP_FDR, adjP_Bonf) %>% 
  filter(adjP_Bonf < 0.05) %>% 
  arrange(SUPER_PATHWAY) %>% nrow()

# save annotated results file
write.table(final_out_annot, paste0(working, "assoc_res_10k_v2/2023_glm/HELIOS_Metabolon1073_glm_logmeancIMT8124_AgeSexEthnicitySSID_2023.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

