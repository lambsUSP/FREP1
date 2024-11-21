################################################################################
# Model simulation ----
################################################################################

experiment_sim_model <- function(
    data,
    numGen,
    genotype,
    parameters,
    N.sample = NULL)
{
  
  # parameters
  for (i in 1:length(parameters)) {
    assign(names(parameters)[i],as.numeric(parameters[i]))  
  }
  
  ################################################################################
  # Transition matrix ----
  ################################################################################
  
  # tMatrix for allelic elimination experiments
  tMatrix   = tMatrix_func(genotype,parameters)
  tMatrix.F = tMatrix[[1]]
  tMatrix.M = tMatrix[[2]]
  
  ################################################################################
  # model experiment simulation ----
  ################################################################################
  
  # simulations
  model = single_sim_model(
    numGen,
    genotype,
    tMatrix.F,
    tMatrix.M,
    parameters
  )
  
  # stochastic simulations
  if (!is.null(N.sample)) {
    
    N.sto     = 100
    model.sto = vector("list", length = N.sto)
    
    for (i in 1:N.sto) {
      
      # simulations
      model.sto[[i]] = single_sim_model(
        numGen,
        genotype,
        tMatrix.F,
        tMatrix.M,
        parameters,
        N.sample,
        sto=T
      )
    }    
    
    return(list(model,model.sto))
    
  }else{return(list(model))}
}
