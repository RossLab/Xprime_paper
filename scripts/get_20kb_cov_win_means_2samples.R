# R
# 6.11.21
# analysing coverage: means across windows: 20kb overlapping windows with 5kb step
# This script is for submitting on the cluster, the output should be analysed in RStudio

install.packages("data.table",repos = "http://cran.us.r-project.org")
install.packages("dplyr",repos = "http://cran.us.r-project.org")
install.packages("zoo",repos = "http://cran.us.r-project.org")

library(data.table)
library(dplyr)
library(zoo)

# load in per_base read cov (F)
per_base_cov_f <- read.table('per_base_cov_F.txt', stringsAsFactors = FALSE, header=F, fill=T)
colnames(per_base_cov_f) <- c('id','pos','cov')
per_base_cov_f$cov <- as.numeric(per_base_cov_f$cov)
per_base_cov_f$pos <- as.numeric(per_base_cov_f$pos)

# calc sliding window means
f_win_means <- per_base_cov_f %>%
  group_by(id) %>%
  mutate(cov.win.mean=rollapply(cov,20000,mean,by=20000, fill=NA))
f_win_means_no_nas <- f_win_means[complete.cases(f_win_means), ]
f_win_means_final <- f_win_means_no_nas[,c(1,2,4)]
f_win_means_final <- as.data.frame(f_win_means_final)

# write table
write.table(f_win_means_final, file="F_20-5kb_cov_means.txt", quote=F, sep="\t", row.names=F, col.names=T)

# load in per_base read cov (M)
per_base_cov_m <- read.table('per_base_cov_M.txt', stringsAsFactors = FALSE, header=F, fill=T)
colnames(per_base_cov_m) <- c('id','pos','cov')
per_base_cov_m$cov <- as.numeric(per_base_cov_m$cov)
per_base_cov_m$pos <- as.numeric(per_base_cov_m$pos)

# calc sliding window means
m_win_means <- per_base_cov_m %>%
  group_by(id) %>%
  mutate(cov.win.mean=rollapply(cov,20000,mean,by=20000, fill=NA))
m_win_means_no_nas <- m_win_means[complete.cases(m_win_means), ]
m_win_means_final <- m_win_means_no_nas[,c(1,2,4)]
m_win_means_final <- as.data.frame(m_win_means_final)
# write table
write.table(m_win_means_final, file="M_20-5kb_cov_means.txt", quote=F, sep="\t", row.names=F, col.names=T)

