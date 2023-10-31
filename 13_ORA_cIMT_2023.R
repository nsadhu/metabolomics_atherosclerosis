################################
## Over-representation analysis
###############################


# remove environment variables
rm(list = ls())

# load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)

working = "C:/Users/Nilanjana/OneDrive - Nanyang Technological University/Metabolomics/NS_Helios_Metabolon/"

# get metabolites universal set
df <- fread(paste0(working, "Batch-normImpDataAll_Jan2023_missing25_zerovar_1073metabolites_235repVis2Excl_9pc1OutliersExcl_log_std.csv"), 
            sep = ",",  header = TRUE, stringsAsFactors = FALSE, fill = TRUE, data.table = FALSE)
df[1:5, 1:5]
metabolites <- as.numeric(colnames(df[, -1]))
rm(df)

# read metabolon chemical annotation data
df_chem <- fread("C:/Users/Nilanjana/OneDrive - Nanyang Technological University/Metabolomics/HELIOS data-Metabolome/HELIOS-Metabolon_281222/Chemical Annotation All.csv", 
                 header = TRUE, stringsAsFactors = FALSE, fill = TRUE, data.table = FALSE)
head(df_chem)

library(stringr)
df_chem %>%
  mutate_if(is.character, str_trim) -> df_chem

# select metabolites set, remove unknown pathways
df_chem %>% 
  filter(CHEM_ID %in% metabolites) %>% 
  filter(!SUPER_PATHWAY %in% c("", "Partially Characterized Molecules")) %>%
  dplyr::select(CHEM_ID, CHEMICAL_NAME, SUPER_PATHWAY, SUB_PATHWAY) -> df_total

head(df_total)
table(df_total$SUPER_PATHWAY)

# read regression results
df_reg1 <- fread(paste0(working, "assoc_res_10k_v2/Nov2022_glm/HELIOS_Metabolon883_glm_logmeancIMT8124_AgeSexEthnicitySSID_Nov2022_std.txt"),
                 header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

df_reg2 <- fread(paste0(working, "assoc_res_10k_v2/Nov2022_glm/HELIOS_Metabolon883_glm_logmeancIMT8056_FullModel_Nov2022_std.txt"),
                 header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

# select significant results, remove unknown pathways
df_reg1 %>%
  filter(adjP_Bonf < 0.05) %>%
  dplyr::select(CHEM_ID, Est) %>% 
  left_join(df_reg2, by = "CHEM_ID") %>% 
  filter(P < 0.05) %>% 
  filter(sign(Est.x) == sign(Est.y)) %>% 
  filter(!SUPER_PATHWAY %in% c("Partially Characterized Molecules", "")) %>% 
  dplyr::select(CHEM_ID, CHEMICAL_NAME, SUPER_PATHWAY, SUB_PATHWAY) -> df_signif

head(df_signif)

# randomly sample signif rows without replacement from df_total for expected distribution of metabolites within pathways
df_total %>% 
  group_by(SUB_PATHWAY) %>% 
  summarise(total = n()) -> all_pathways

curr_res = NULL
set.seed(1234)

for(i in 1:10000) {
  sampled_pathway <- sample(df_total$SUB_PATHWAY, nrow(df_signif))
  
  sampled_pathway %>%
    table() %>%
    as_tibble() %>%
    dplyr::rename(SUB_PATHWAY = 1, freq = 2) %>%
    right_join(all_pathways, by = "SUB_PATHWAY") %>% 
    mutate(iter = i) %>%
    dplyr::select(iter, SUB_PATHWAY, total, freq) %>% 
    rbind(curr_res) -> curr_res
}

# get observed distribution of metabolites within pathways 
df_signif %>% 
  group_by(SUB_PATHWAY) %>% 
  summarise(hits = n()) -> signif_pathways

  # check over-representation
curr_res %>% 
  mutate(freq = replace_na(freq, 0)) %>% 
  group_by(SUB_PATHWAY) %>% 
  summarise(total = mean(total), 
            exp_values = list(freq)) %>% 
  left_join(signif_pathways, by = "SUB_PATHWAY") %>% 
  filter(is.na(hits) == FALSE & total > 5) %>% 
  ungroup() %>% 
  rowwise() %>%
  mutate(p = 1-ecdf(unlist(exp_values))(hits)) %>% 
  dplyr::select(SUB_PATHWAY, total, hits, p) -> ora_res

ora_res %>% 
  ungroup() %>% 
  mutate(fdr = p.adjust(p, method = "BH")) %>% 
  filter(p < 0.05) %>% arrange(p)

  
