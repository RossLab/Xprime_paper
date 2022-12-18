# R
# analysing coverage: means cov per scaffold, 2 samples
# This script is for submitting on the cluster, the output should be analysed in RStudio

install.packages("data.table",repos = "http://cran.us.r-project.org")
install.packages("dplyr",repos = "http://cran.us.r-project.org")
install.packages("zoo",repos = "http://cran.us.r-project.org")

library(data.table)
library(dplyr)
library(zoo)

# load in per_base read cov
per_base_cov <- read.table('per_base_cov_M.txt', stringsAsFactors = FALSE, header=F, fill=T)
colnames(per_base_cov) <- c('id','pos','cov')
per_base_cov$cov <- as.numeric(per_base_cov$cov)
per_base_cov$pos <- as.numeric(per_base_cov$pos)

# calc means per scaffold
ctg_means <- per_base_cov %>%
  group_by(id) %>%
  summarise(mean_cov=mean(cov))
ctg_means_final <- as.data.frame(ctg_means)

# write table
write.table(ctg_means_final, file="ctg_cov_means_XO.txt", quote=F, sep="\t", row.names=F, col.names=T)

# load in per_base read cov
per_base_cov <- read.table('per_base_cov_F.txt', stringsAsFactors = FALSE, header=F, fill=T)
colnames(per_base_cov) <- c('id','pos','cov')
per_base_cov$cov <- as.numeric(per_base_cov$cov)
per_base_cov$pos <- as.numeric(per_base_cov$pos)

# calc means per scaffold
ctg_means <- per_base_cov %>%
  group_by(id) %>%
  summarise(mean_cov=mean(cov))
ctg_means_final <- as.data.frame(ctg_means)

# write table
write.table(ctg_means_final, file="ctg_cov_means_XpX.txt", quote=F, sep="\t", row.names=F, col.names=T)
