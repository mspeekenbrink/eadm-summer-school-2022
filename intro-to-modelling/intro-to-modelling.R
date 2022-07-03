## ----setup, include=FALSE-----------------------------------------------------
options(htmltools.dir.version = FALSE)
knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(fig.retina = 3, dev='svg', out.width = "90%", fig.align="center", message=FALSE, warning=FALSE, purl=FALSE)

## ----load-packages, include = FALSE, purl=TRUE--------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

## ---- cache=TRUE, purl = TRUE-------------------------------------------------
set.seed(20220104)
nP <- 1000; nT <- 100 # number of people and trials
eta <- .2; tau <- 1 # parameters
means <- c(-1, 1, 2, 2.5); sds <- c(3,3,1,1) # option means and SDs
simdat <- data.frame() # for storing results
for(i in 1:nP) {
  E <- matrix(0, nrow=nT+1, ncol = length(means)); D <- matrix(0, nrow=nT, ncol= length(means)); R <- rep(0.0, length=nT) # for storing results
  for(t in 1:nT) {
    # choose option
    k <- sample(1:4, size = 1, prob=exp(tau*E[t,])/sum(exp(tau*E[t,])))
    D[t,k] <- 1
    # get reward
    R[t] <- rnorm(1, means[k], sds[k])
    # update expectancies
    E[t+1,] <- E[t,]
    E[t+1,k] <- E[t,k] + eta*(R[t] - E[t,k])
  }
  E <- E[1:nT,]
  simdat <- rbind(simdat,
               data.frame(id = i, trial = 1:nT, E1 = E[,1], E2 = E[,2], E3 = E[,3], E4 = E[,4], D1 = D[,1], D2 = D[,2], D3 = D[,3], D4 = D[,4], R = R))
}

## ----likelihood-definition, purl=TRUE-----------------------------------------
# define a function to compute the likelihood
loglikelihood <- function(theta, data) {
  loglik <- 0.0
  eta <- theta[1]
  tau <- theta[2]
  for(i in unique(data$id)) {
    tdat <- subset(data, id == i)
    E <- rep(0, 4)
    for(t in 1:nrow(tdat)) {
      pD <- exp(tau*E)/sum(exp(tau*E))
      choice <- which(tdat[t,paste0("D", 1:4)] == 1)
      loglik <- loglik + log(pD[choice])
      E[choice] <- E[choice] + eta*(tdat$R[t] - E[choice])  
    }
  }
  return(loglik)
}

## ----max-likelihood-grid-search, cache=TRUE, purl=TRUE------------------------
# perform a grid search
theta <- expand.grid(eta = seq(0,1,length=20), tau = seq(0,4, length=20))
loglik <- rep(0, nrow(theta))
for(i in 1:nrow(theta)) {
  loglik[i] <- loglikelihood(as.numeric(theta[i,]), subset(simdat, id <= 20))
}

## ----Nelder-Mead, cache=TRUE, purl=TRUE---------------------------------------
# most numerical optimization routines are geared towards minimization
negloglikelihood <- function(theta, data) {
  -1*loglikelihood(theta, data)
}
optNM20 <- optim(par = c(.1, 1), fn = negloglikelihood, method = "Nelder-Mead", data=subset(simdat, id <= 20))
optNM20

## ----DEoptim, cache=TRUE, purl=TRUE-------------------------------------------
library(DEoptim)
optDEOPT20 <- DEoptim(negloglikelihood, lower = c(0,0), upper = c(1, 10), data=subset(simdat, id <= 20), control = DEoptim.control(itermax=30)) # note: you should set itermax higher!
optDEOPT20

## ----LBFGSB, cache=TRUE, purl=TRUE--------------------------------------------
optLBFGSB20 <- optim(par = c(.1, 1), fn = negloglikelihood, method = "L-BFGS-B", lower=c(0,0), upper=c(1, Inf), data=subset(simdat, id <= 20))
optLBFGSB20

## ----max-likelihood-transformed, cache=TRUE, purl=TRUE------------------------
# define a function to compute the likelihood
tnegloglikelihood <- function(theta, data) {
  eta <- 1/(1+exp(-theta[1]))
  tau <- exp(theta[2])
  negloglikelihood(c(eta, tau), data)
}
optTNM20 <- optim(par = c(log(.1/(1-.1)), log(1)), fn = tnegloglikelihood, method = "Nelder-Mead", data=subset(simdat, id <= 20))
optTNM20
c(1/(1+exp(-optTNM20$par[1])), exp(optTNM20$par[2]))

## ----perseverence-negloglikelihood, purl=TRUE---------------------------------
negloglikelihood_pers <- function(theta, data) {
  loglik <- 0.0
  eta <- theta[1]
  tau <- theta[2]
  phi <- theta[3]
  gamma <- theta[4]
  for(i in unique(data$id)) {
    tdat <- subset(data, id == i)
    E <- rep(0, 4)
    B <- rep(0, 4)
    for(t in 1:nrow(tdat)) {
      pD <- exp(tau*E + phi*B)/sum(exp(tau*E + phi*B))
      choice <- which(tdat[t,paste0("D", 1:4)] == 1)
      loglik <- loglik + log(pD[choice])
      E[choice] <- E[choice] + eta*(tdat$R[t] - E[choice])
      B <- B + gamma*(as.numeric(1:4 == choice) - B)
    }
  }
  return(-loglik)
}

## ----LBFGSB-pers, cache=TRUE, purl=TRUE---------------------------------------
optLBFGSB20_pers <- optim(par = c(.1, 1, .1, .5), fn = negloglikelihood_pers, method = "L-BFGS-B", lower=c(0,0,0,0), upper=c(1, Inf, Inf, 1), data=subset(simdat, id <= 20))
optLBFGSB20_pers

## ----fit-models-to-1-person, purl=TRUE----------------------------------------
optLBFGSB1 <- optim(par = c(.1, 1), fn = negloglikelihood, method = "L-BFGS-B", lower=c(0,0), upper=c(1, Inf), data=subset(simdat, id==10))
optLBFGSB1_pers <- optim(par = c(.1, 1, .1, .5), fn = negloglikelihood_pers, method = "L-BFGS-B", lower=c(0,0,0,0), upper=c(1, Inf, Inf, 1), data=subset(simdat, id==10))

## ----compute-information-criteria, purl=TRUE----------------------------------
AICs <- 2*optLBFGSB1$value + 2*2
AICp <- 2*optLBFGSB1_pers$value + 2*4
BICs <- 2*optLBFGSB1$value + 2*log(100)
BICp <- 2*optLBFGSB1_pers$value + 4*log(100)
matrix(c(AICs, AICp, BICs, BICp), ncol=2, dimnames=list(c("simple","perseverence"),c("AIC", "BIC")))

## ----fd-hessian, purl=TRUE----------------------------------------------------
hess <- numDeriv::hessian(tnegloglikelihood, optTNM20$par, data=subset(simdat, id <= 20))
se <- sqrt(diag(solve(hess)))
conf <- rbind(lower = optTNM20$par - qnorm(.975)*se,
                 upper = optTNM20$par + qnorm(.975)*se)
colnames(conf) <- c("eta","tau")
conf

## ---- purl=TRUE---------------------------------------------------------------
cbind(1/(1+exp(-conf[,1])), exp(conf[,2]))

