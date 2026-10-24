rm(list = ls())
library(MASS)
library(glmnet)
library(grpreg)
library(ipred)
library(randomForest)
source("utils.R")
source("gdat.R")
source("prepare.R")
source("prepare_reg_bag.R")
source("method_loc_lasso.R")
source("method_glasso.R")
source("method_hlasso.R")
source("method_prop_admm.R")
source("method_oracle.R")
source("inference.R")
source("wloss_crossfit.R")
source("simu.R")
# ---------------------------------------------------------
# Simulation settings
# ---------------------------------------------------------
K <- 4
p <- 20
n <- K*600
rho <- 0.5
treat_prob <- 0.3
balance <- TRUE

if(balance){
taskID <- paste0(
  "Balance_",
  K, "K_",
  p, "p_",
  n, "n_",
  treat_prob, "treat"
)}else{
taskID <- paste0(
  "Imbalance_",
  K, "K_",
  p, "p_",
  n, "n_",
  treat_prob, "treat"
)
}

# ---------------------------------------------------------
# True beta
# beta0 is p x K here.
# gdat() internally transposes it to K x p.
# ---------------------------------------------------------
b1 <- c(-1, 1, 1, -1,  1)
b1 <- c(b1, rep(0, p - length(b1)))

b2 <- c(-1, 1, 1, -1,  1)
b2 <- c(b2, rep(0, p - length(b2)))

b3 <- c(-1, 1, 1,  1, -1)
b3 <- c(b3, rep(0, p - length(b3)))

b4 <- c(-1, 1, 1,  1, -1)
b4 <- c(b4, rep(0, p - length(b4)))

signal <- 1

beta0 <- cbind(b1, b2, b3, b4) * signal

# ---------------------------------------------------------
# Contrast a for inference
#
# a must be K x p.
#
# Example 1:
#   psi = beta_{1,4} - beta_{3,4}
#
# True value under current beta0:
#   beta_{1,4} = 0.8
#   beta_{3,4} = 1.6
#   psi0 = -0.8
# ---------------------------------------------------------
a <- matrix(1, nrow = K, ncol = p)/sqrt(sum(beta0!=0))
 

# ---------------------------------------------------------
# Run simulation
# ---------------------------------------------------------
Nsimu <- 1000
seedstart <- 2025

simu(
  n = n,
  K = K,
  p = p,
  rho = rho,
  beta0 = beta0,
  treat_prob = treat_prob,
  balance = balance,
  nui_app = "bagging", ## code only support bagging
  seedstart = seedstart,
  Nsimu = Nsimu,
  taskID = taskID,
  a = a,
  save_raw_each = FALSE,
  raw_dir = paste0(taskID, "_RawEachSim")
)