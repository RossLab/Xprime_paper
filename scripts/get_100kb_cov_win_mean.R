# R
# analysing coverage: means across windows
# This script is for submitting on the cluster, the output should be analysed in RStudio

install.packages("data.table",repos = "http://cran.us.r-project.org")
install.packages("dplyr",repos = "http://cran.us.r-project.org")
install.packages("zoo",repos = "http://cran.us.r-project.org")

library(data.table)
library(dplyr)
library(zoo)

# load in per_base read cov (F)
per_base_cov <- read.table('per_base_cov.txt', stringsAsFactors = FALSE, header=F, fill=T)
colnames(per_base_cov) <- c('id','pos','cov')
per_base_cov$cov <- as.numeric(per_base_cov$cov)
per_base_cov$pos <- as.numeric(per_base_cov$pos)

# calc sliding window means
win_means <- per_base_cov %>%
  group_by(id) %>%
  mutate(cov.win.mean=rollapply(cov,100000,mean,by=100000, fill=NA))
win_means_no_nas <- win_means[complete.cases(win_means), ]
win_means_final <- win_means_no_nas[,c(1,2,4)]
win_means_final <- as.data.frame(win_means_final)

# write table
write.table(win_means_final, file="100kb_cov_win_means.txt", quote=F, sep="\t", row.names=F, col.names=T)
