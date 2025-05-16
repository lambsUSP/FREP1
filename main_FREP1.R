library(lhs)
library(foreach)
library(doParallel)
library(plotrix)
library(readxl)
library(mcmc)
library(miscTools)
library(matrixStats)

################################################################################
# Functions ----
################################################################################

source("loading_data_experiment.R")
source("plotting_data_experiment.R")
source("tMatrix_func.R")
source("single_sim_model.R")
source("maternal_deposition.R")
source("experiment_sim_model.R")
source("obj_func.R")
source("metropolis_MCMC.R")

################################################################################
# Data ----
################################################################################

data <- loading_data_experiment()

# number of generations
numGen = nrow(data)

# genotype names
genotype = c(
  "aaBB","aaBC","aaBD","aaBE","aaCC","aaCD","aaCE","aaDD","aaDE","aaEE",
  "aABB","aABC","aABD","aABE","aACC","aACD","aACE","aADD","aADE","aAEE",
  "AABB","AABC","AABD","AABE","AACC","AACD","AACE","AADD","AADE","AAEE"
)

################################################################################
# Model simulation ----
################################################################################

# Modeling framework:
# A: vasa-Cas9
# a: no vasa-Cas9
# B: DSRed+gRNA+Q
# C: GFP+Q
# D: GFP+L
# E: NHEJ

# parameters
parameters <- data.frame(
  r.driving    = 0.8, # FREP1L -> FREP1Q
  r.NHEJ       = 0.1, # FREP1L -> NHEJ
  p.1          = 0.9
)

# model simulation for data fitting
aux.model = experiment_sim_model(
  data,
  numGen,
  genotype,
  parameters,
  N.sample = 300/2
)

model     = aux.model[[1]]
model.sto = aux.model[[2]]

################################################################################
# Visualizing Data ----
################################################################################

plotting_data_experiment(
  data,
  parameters,
  model,
  model.sto
)

obj_func(data,numGen,genotype,parameters)

################################################################################
# Model fit ----
################################################################################

# running MCMC
print("Running MCMC:")
chain = metropolis_MCMC(
  data,
  numGen,
  genotype,
  startvalue = parameters,
  iterations = 1e5,
  scale = c(0.02)
)

# acceptance ratio
burnin = 1e3
acceptance_ratio = 1-mean(duplicated(chain[-(1:burnin),]))
acceptance_ratio

# final result from MCMC
parameters_MCMC   = parameters
parameters_MCMC[] = colQuantiles(chain[burnin:nrow(chain),1:length(parameters)],probs = c(0.5))
aux.model_MCMC    = experiment_sim_model(
  data,
  numGen,
  genotype,
  parameters = parameters_MCMC,
  N.sample = 300/2
)
model_MCMC     = aux.model_MCMC[[1]]
model.sto_MCMC = aux.model_MCMC[[2]]

# Comparing likelihood values
obj_func(data,numGen,genotype,parameters)
obj_func(data,numGen,genotype,parameters_MCMC)

# Plotting final results
plotting_data_experiment(
  data       = data,
  parameters = round(parameters_MCMC,2),
  model      = model_MCMC,
  model.sto  = model.sto_MCMC
)

# Plotting chains
jpeg("FREP1_chains.jpeg", width = 5, height = 10, units = 'in', res = 100)
par(mai = c(0.5, 0.5, 0.4, 0.2), mgp = c(1.5,0.5,0))
layout(matrix(c(1:(2*length(parameters))), ncol = 2, byrow = T))
for (i in 1:length(parameters)) {
  plot(chain[,i],type = "l", ylab = names(parameters)[i], xlab = "", ylim = c(0,1.5))  
  hist(chain[,i], xlab = names(parameters)[i], main = "")
}
dev.off()

# 95% credible intervals
parameters_MCMC
round(quantile(chain[burnin:nrow(chain),1],c(0.5,0.025,0.975)),digits = 2)
round(quantile(chain[burnin:nrow(chain),2],c(0.5,0.025,0.975)),digits = 2)
round(quantile(chain[burnin:nrow(chain),3],c(0.5,0.025,0.975)),digits = 2)
