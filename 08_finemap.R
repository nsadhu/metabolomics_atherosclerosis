
##############################
## Fine-mapping
#############################


# remove environment variables
rm(list = ls())

# load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(susieR)


#read metabolite summary stats
meta <- fread(".../meta_100006370.res_scale.loco.mlma", 
            header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

peak = 110371275
win = c((peak-1500000) : (peak+1500000)) #3mb window around peak snp

meta %>%
  filter(Chr == 11 & bp %in% win) -> meta_finemap

# obtain correlation matrix of snps in this region as read as R_sg10K

# finemap 
fitted_rss <- susie_rss(bhat = meta_finemap$b_finemap, shat = meta_finemap$se, R = R_sg10k, n=1867, L = 10)
summary(fitted_rss)$cs
fitted_rss$sets

i  <- fitted_rss$sets$cs[[1]]
z <- cbind(i,meta_finemap$b_finemap[i], meta_finemap$p[i], fitted_rss$pip[i])
colnames(z) <- c('position', 'beta', 'p-value', 'PIP')
z[order(z[,2], decreasing = TRUE),]

susie_plot(fitted_rss, y="PIP")
