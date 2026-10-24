rm(list = ls())
library(MASS)
library(glmnet)
library(grpreg)
library(ipred)
library(randomForest)
source("utils.R")
source("gdat.R")
source("prepare_dc.R")
source("prepare_reg_bag.R")
source("method_loc_lasso.R")
source("method_prop_admm.R")
source("method_oracle.R")
source("wloss_crossfit.R")
source("simu.R")
source("simu_nuisance_compare.R")
# ---------------------------------------------------------
# set model
# --------- varing
p <- 20
treat_prob <- 0.3
balance <- !FALSE

# --------- fixed
K <- 4
n <- K*600
rho <- 0.5

# --------- set real effects 
b1 = c(-0.6, -0.6, 0.6, 1.0,   0,   0); b1 = c(b1,rep(0,p-length(b1)))
b2 = c(-0.6, -0.6, 0.6,   0, 0.8,   0); b2 = c(b2,rep(0,p-length(b2)))
b3 = c(-0.6, -0.6, 0.6,   0, 0.8,   0); b3 = c(b3,rep(0,p-length(b3)))
b4 = c(-0.6, -0.6, 0.6,   0,   0, 1.0); b4 = c(b4,rep(0,p-length(b4)))
beta0=cbind(b1,b2,b3,b4) 
 

# --------- set taskID
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
Nsimu <- 500
seedstart <- 2025
simu_nuisance_compare(
  n = n,
  K = K,
  p = p,
  rho = rho,
  beta0 = beta0,
  treat_prob = treat_prob,
  balance = balance,
  seedstart = seedstart,
  Nsimu = Nsimu,
  taskID = taskID,
  nuisance_methods = c("known", "bagging", "rf", "glm_lasso"),
  nuisance_modes = c("pooled", "byG"),
  nfold = 10,
  use_weighted = TRUE
)