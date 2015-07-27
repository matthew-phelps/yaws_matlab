## intro
rm(list = ls())
mac <- "/Users/Matthew/Google Drive/Copenhagen/DK Cholera/CPH"
pc <- 'C:/Users/wrz741/Google Drev/Yaws project/MATLAB/Uncertainty Analysis'

setwd(pc)
library (SDMTools) # used for calculating weighted standard deviations
library (MASS) # used for fiting distributions
library (vcd) # used for goodness-of-fit tests
library (ggplot2)
library (stats)
library (reshape) # for renaming variables\
library (R.matlab) # for MATLAB files
library (LaplacesDemon)

install.packages(pkgs="C:/Users/wrz741/Google Drev/R Statisical tutorials/LaplacesDemon_14.11.12.tar.gz", repos=NULL, type="source")


beta<-readMat("betaONEvect.mat")$betaONEvect
epsilon <-readMat("epsilonVect.mat")$epsilonvect


beta_500 <- as.data.frame(beta[500:nrow(beta)])
epsilon_500 <- as.data.frame(epsilon[500:nrow(epsilon)])



# 1 treatment round -------------------------------------------------------


p.interval(beta_500, HPD=F, MM=F, plot=T)
p.interval(epsilon_500, HPD=F, MM=F, plot=T)


