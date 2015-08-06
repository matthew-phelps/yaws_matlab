setwd("/Users/lauraskrip/Google Drive/Yaws project/yaws_matlab/Laura's Folder")
a <- 1
b <- 1
posterior.unif.prior <- rbeta(10000,shape1=(.855*250+a),shape2=(250-.855*250+b))
write.csv(posterior.unif.prior,'beta.csv')

#######  OLD ######
# ## Determine variance
# mean <- .855
# variance <- ((mean-.782)*sqrt(250)/1.96)^2
# 
# estBetaParams <- function(mu, var) {
#        alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
#        beta <- alpha * (1 / mu - 1)
#        return(params = list(alpha = alpha, beta = beta))
#    }
# estBetaParams(mean, variance)