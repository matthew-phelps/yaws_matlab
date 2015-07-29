## Determine variance
mean <- .855
variance <- ((mean-.782)*sqrt(250)/1.96)^2

estBetaParams <- function(mu, var) {
       alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
       beta <- alpha * (1 / mu - 1)
       return(params = list(alpha = alpha, beta = beta))
   }
estBetaParams(mean, variance)