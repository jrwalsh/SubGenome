library(mixtools)

simulate <- function(lambda=0.3, mu=c(0, 4), sd=c(1, 1), n.obs=10^5) {
  x1 <- rnorm(n.obs, mu[1], sd[1])
  x2 <- rnorm(n.obs, mu[2], sd[2])
  return(ifelse(runif(n.obs) < lambda, x1, x2))
}

x <- simulate()

model <- normalmixEM(x=x, k=2)
index.lower <- which.min(model$mu)  # Index of component with lower mean

find.cutoff <- function(proba=0.5, i=index.lower) {
  ## Cutoff such that Pr[drawn from bad component] == proba
  f <- function(x) {
    proba - (model$lambda[i]*dnorm(x, model$mu[i], model$sigma[i]) /
               (model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) + model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2])))
  }
  return(uniroot(f=f, lower=-10, upper=10)$root)  # Careful with division by zero if changing lower and upper
}

cutoffs <- c(find.cutoff(proba=0.5), find.cutoff(proba=0.75))  # Around c(1.8, 1.5)

hist(x)
abline(v=cutoffs, col=c("red", "blue"), lty=2)

#https://stats.stackexchange.com/questions/57993/how-to-explain-how-i-divided-a-bimodal-distribution-based-on-kernel-density-esti
