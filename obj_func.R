###############################################################################
# objective function ---
###############################################################################

obj_func <- function(data,numGen,genotype,parameters){
  
  # model simulation
  aux.model = experiment_sim_model(
    data,
    numGen,
    genotype,
    parameters,
    N.sample = NULL
  )
  
  # model output
  model      = aux.model[[1]]
  parameters = as.numeric(parameters)
  
  # initializing ll.data
  ll.data = 0
  
  if (any(parameters>1) || any(parameters<0)) {
    
    # priors
    ll.prior = sum(
      max(-1e10, dunif(parameters[1], min = 0, max = 1, log = T)), # r.driving
      max(-1e10, dunif(parameters[2], min = 0, max = 1, log = T)), # r.NHEJ
      max(-1e10, dunif(parameters[3], min = 0, max = 1, log = T)), # r.NHEJ
      max(-1e10, dunif(parameters[1] + parameters[2], min = 0, max = 1, log = T))  # r.driving + r.NHEJ
    )
  }else{
    
    # for (i in c(1:numGen)) {
    for (i in c(2:numGen)) {
      
      ll.data = ll.data + 
        
        # cage 1
        dmultinom(x    = c(data$data_phenotype_LDonly_cage1[i],data$data_phenotype_GFPLonly_cage1[i],data$data_phenotype_LDGFPL_cage1[i]),
                  # size = data$data_screened_cage1[i], 
                  prob = c(model$LDonly[i],model$GFPLonly[i],model$LDGFPL[i]), 
                  log  = TRUE) +
        dmultinom(x    = c(data$data_allelefreq_L224_cage1[i],data$data_allelefreq_Q224_cage1[i],data$data_allelefreq_NHEJ_cage1[i]),
                  # size = data$data_genotyped_cage1[i],
                  prob = c(model$L224[i],model$Q224[i],model$NHEJ[i]),
                  log  = TRUE) +
        dbinom(   x    = data$data_phenotype_V9_cage1[i],
                  size = data$data_screened_cage1[i],
                  prob = model$V9[i],
                  log  = TRUE) +
        
        # cage 2
        dmultinom(x    = c(data$data_phenotype_LDonly_cage2[i],data$data_phenotype_GFPLonly_cage2[i],data$data_phenotype_LDGFPL_cage2[i]),
                  # size = data$data_screened_cage2[i], 
                  prob = c(model$LDonly[i],model$GFPLonly[i],model$LDGFPL[i]), 
                  log  = TRUE) +
        dmultinom(x    = c(data$data_allelefreq_L224_cage2[i],data$data_allelefreq_Q224_cage2[i],data$data_allelefreq_NHEJ_cage2[i]),
                  # size = data$data_genotyped_cage2[i],
                  prob = c(model$L224[i],model$Q224[i],model$NHEJ[i]),
                  log  = TRUE) +
        dbinom(   x    = data$data_phenotype_V9_cage2[i],
                  size = data$data_screened_cage2[i],
                  prob = model$V9[i],
                  log  = TRUE) +
        
        # cage 3
        dmultinom(x    = c(data$data_phenotype_LDonly_cage3[i],data$data_phenotype_GFPLonly_cage3[i],data$data_phenotype_LDGFPL_cage3[i]),
                  # size = data$data_screened_cage3[i], 
                  prob = c(model$LDonly[i],model$GFPLonly[i],model$LDGFPL[i]), 
                  log  = TRUE) +
        dmultinom(x    = c(data$data_allelefreq_L224_cage3[i],data$data_allelefreq_Q224_cage3[i],data$data_allelefreq_NHEJ_cage3[i]),
                  # size = data$data_genotyped_cage3[i],
                  prob = c(model$L224[i],model$Q224[i],model$NHEJ[i]),
                  log  = TRUE) +
        dbinom(   x    = data$data_phenotype_V9_cage3[i],
                  size = data$data_screened_cage3[i],
                  prob = model$V9[i],
                  log  = TRUE)
    }
    
    # priors
    ll.prior = sum(
      max(-1e10, dunif(parameters[1], min = 0, max = 1, log = T)), # r.driving
      max(-1e10, dunif(parameters[2], min = 0, max = 1, log = T)), # r.NHEJ
      max(-1e10, dunif(parameters[3], min = 0, max = 1, log = T)), # r.NHEJ
      max(-1e10, dunif(parameters[1] + parameters[2], min = 0, max = 1, log = T))  # r.driving + r.NHEJ

    )
  }
  
  # final Likelihood
  ll <- ll.data + ll.prior
  
  return(ll)
}
