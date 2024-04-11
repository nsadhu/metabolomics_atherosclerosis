
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

# read metabolite data
df <- fread(".../metabolite_data.csv", 
            sep = ",",  header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

# read phenotype data
pheno <- fread(".../phenotype_data.txt", 
                 sep = "\t",  header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)


# subset and factorize phenotypes 
pheno %>%
  filter(Ethnicity != "O") -> df_pheno #remove ethnicity that is not Chinese, Malay, Indian

df_pheno$Ethnicity = as.factor(df_pheno$Ethnicity)
df_pheno$Sex = as.factor(df_pheno$Sex)
df_pheno$Batch = as.factor(df_pheno$Batch)

# drop na from outcome to analyze
df_pheno %>% 
  drop_na(ln_mean_cIMT, Age, Sex, Ethnicity, Batch) -> df_pheno

# subset and sort samples in same order
df %>% 
  filter(PARENT_SAMPLE_NAME %in% df_pheno$PARENT_SAMPLE_NAME) %>% 
  arrange(PARENT_SAMPLE_NAME) -> df

table(df$PARENT_SAMPLE_NAME == df_pheno$PARENT_SAMPLE_NAME) #should be TRUE if in same order

df = df[, -1] # retain all columns except sample names

# remove metabolites with only one value i.e. column variance = 0
df %>% 
  dplyr::select_if(~var(., na.rm = TRUE) != 0) -> df

# save metabolite names
CHEM_ID <- colnames(df)
colnames(df) = NULL
df <- as.matrix(df)

# run linear regression
lm.formula <- as.formula('df_pheno$ln_mean_cIMT ~ df[, i] +
                         df_pheno$Age + df_pheno$Sex + df_pheno$Ethnicity + df_pheno$Batch')

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
df_chem <- fread(".../Chemical Annotation All.csv", 
                 header = TRUE, stringsAsFactors = FALSE, fill = TRUE, data.table = FALSE)

library(stringr)
df_chem %>%
  mutate_if(is.character, str_trim) -> df_chem

df_chem$CHEM_ID = as.character(df_chem$CHEM_ID)

# merge chemical annotation information with regression results
final_out %>% 
  mutate(adjP_FDR = p.adjust(P, method = "BH"), 
         adjP_Bonf = p.adjust(P, method = "bonferroni")) %>%
  left_join(df_chem, by = "CHEM_ID") -> final_out_annot


