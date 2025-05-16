###############################################################################
# Single model simulation ----
###############################################################################

single_sim_model <- function(
    numGen,
    genotype,
    tMatrix.F,
    tMatrix.M,
    parameters,
    N.sample=NULL,
    sto=NULL){
  
  # parameters
  for (i in 1:length(parameters)) {
    assign(names(parameters)[i],as.numeric(parameters[i]))  
  }
  
  ################################################################################
  # model variable for experiment simulation ----
  ################################################################################
  
  # variables for the model simulation
  model = data.frame(
    
    # phenotype
    LD       = c(rep(NA,numGen)),
    GFPL     = c(rep(NA,numGen)),
    V9       = c(rep(NA,numGen)),
    
    # phenotype
    LDonly   = c(rep(NA,numGen)),
    GFPLonly = c(rep(NA,numGen)),
    LDGFPL   = c(rep(NA,numGen)),
    
    
    # allele frequencies
    Q224     = c(rep(NA,numGen)),
    L224     = c(rep(NA,numGen)),
    NHEJ     = c(rep(NA,numGen))
    
  )
  
  ################################################################################
  # initial condition ----
  ################################################################################
  
  # # F0:
  # # females: vas-cas9 and GFP+L
  # # males: DSRed+gRNA+Q
  # ic.F1 = list(
  #   females = c(
  #     aaBB = 0.00,   aaBC = 0.00,   aaBD = 0.00,   aaBE = 0.00,
  #     aaCC = 0.00,   aaCD = 0.00,   aaCE = 0.00,
  #     aaDD = 0.00,   aaDE = 0.00,
  #     aaEE = 0.00,
  #     aABB = 0.00,   aABC = 0.00,   aABD = 0.00,   aABE = 0.00,
  #     aACC = 0.00,   aACD = 0.00,   aACE = 0.00,
  #     aADD = 0.00,   aADE = 0.00,
  #     aAEE = 0.00,
  #     AABB = 0.00,   AABC = 0.00,   AABD = 0.00,   AABE = 0.00,
  #     AACC = 0.00,   AACD = 0.00,   AACE = 0.00,
  #     AADD = 1.00,   AADE = 0.00,
  #     AAEE = 0.00
  #   ),
  #   males = c(
  #     aaBB = 1.00,   aaBC = 0.00,   aaBD = 0.00,   aaBE = 0.00,
  #     aaCC = 0.00,   aaCD = 0.00,   aaCE = 0.00,
  #     aaDD = 0.00,   aaDE = 0.00,
  #     aaEE = 0.00,
  #     aABB = 0.00,   aABC = 0.00,   aABD = 0.00,   aABE = 0.00,
  #     aACC = 0.00,   aACD = 0.00,   aACE = 0.00,
  #     aADD = 0.00,   aADE = 0.00,
  #     aAEE = 0.00,
  #     AABB = 0.00,   AABC = 0.00,   AABD = 0.00,   AABE = 0.00,
  #     AACC = 0.00,   AACD = 0.00,   AACE = 0.00,
  #     AADD = 0.00,   AADE = 0.00,
  #     AAEE = 0.00
  #   )
  # )
  # 
  # # matting
  # sim.F = tMatrix.F*array(outer(ic.F1$females,ic.F1$males), dim = c(length(ic.F1$females),length(ic.F1$males),length(ic.F1$females)))
  # sim.M = tMatrix.M*array(outer(ic.F1$females,ic.F1$males), dim = c(length(ic.F1$females),length(ic.F1$males),length(ic.F1$males)))
  # 
  # # maternal deposition
  # sim.F = maternal_deposition(sim.F,genotype,parameters)
  # sim.M = maternal_deposition(sim.M,genotype,parameters)
  # 
  # # deterministic simulation
  # det.F = apply(sim.F, 3, sum)
  # det.M = apply(sim.M, 3, sum)
  # 
  # # avoiding small negative numbers
  # det.F[det.F<0] = 0
  # det.M[det.M<0] = 0
  # 
  # # stochastic simulation
  # if (!is.null(sto)) {
  #   det.F = setNames(as.numeric(rmultinom(1,N.sample,prob = det.F)/N.sample), genotype)
  #   det.M = setNames(as.numeric(rmultinom(1,N.sample,prob = det.M)/N.sample), genotype)
  # }
  # 
  # # normalizing values (progeny from F0)
  # det.F = det.F/sum(det.F)
  # det.M = det.M/sum(det.M)
  # 
  # # G0: homozygous GFP+L + progeny from F0 (derived from F0 crosses)
  # det.F["aaDD"] = 1
  # det.M["aaDD"] = 1
  # 
  # # normalizing values again
  # det.F = det.F/sum(det.F)
  # det.M = det.M/sum(det.M)
  # 
  # # setting up initial condition
  # ic = list(
  #   females = det.F,
  #   males   = det.M
  # )
  
  ############################################################################
  ic = list(
    females = c(
      aaBB = 0.00,   aaBC = 0.00,   aaBD = 0.00,   aaBE = 0.00,
      aaCC = 0.00,   aaCD = 0.00,   aaCE = 0.00,
      aaDD = 0.50,   aaDE = 0.00,
      aaEE = 0.00,
      aABB = 0.00,   aABC = 0.00,   aABD = 0.50,   aABE = 0.00,
      aACC = 0.00,   aACD = 0.00,   aACE = 0.00,
      aADD = 0.00,   aADE = 0.00,
      aAEE = 0.00,
      AABB = 0.00,   AABC = 0.00,   AABD = 0.00,   AABE = 0.00,
      AACC = 0.00,   AACD = 0.00,   AACE = 0.00,
      AADD = 0.00,   AADE = 0.00,
      AAEE = 0.00
    ),
    males = c(
      aaBB = 0.00,   aaBC = 0.00,   aaBD = 0.00,   aaBE = 0.00,
      aaCC = 0.00,   aaCD = 0.00,   aaCE = 0.00,
      aaDD = 0.50,   aaDE = 0.00,
      aaEE = 0.00,
      aABB = 0.00,   aABC = 0.00,   aABD = 0.50,   aABE = 0.00,
      aACC = 0.00,   aACD = 0.00,   aACE = 0.00,
      aADD = 0.00,   aADE = 0.00,
      aAEE = 0.00,
      AABB = 0.00,   AABC = 0.00,   AABD = 0.00,   AABE = 0.00,
      AACC = 0.00,   AACD = 0.00,   AACE = 0.00,
      AADD = 0.00,   AADE = 0.00,
      AAEE = 0.00
    )
  )
  ##############################################################################
  
  # storing values
  aux.LD       = grep("B", genotype, value=TRUE)
  aux.GFP      = unique(c(grep("C", genotype, value=TRUE),grep("D", genotype, value=TRUE),grep("E", genotype, value=TRUE)))
  aux.Cas9     = grep("A", genotype, value=TRUE)
  aux.LDonly   = grep("BB", genotype, value=TRUE)
  aux.GFPLonly = c(grep("CC", genotype, value=TRUE),grep("CD", genotype, value=TRUE),grep("CE", genotype, value=TRUE),grep("DD", genotype, value=TRUE),grep("DE", genotype, value=TRUE),grep("EE", genotype, value=TRUE))
  aux.LDGFPL   = c(grep("BC", genotype, value=TRUE),grep("BD", genotype, value=TRUE),grep("BE", genotype, value=TRUE))
  
  # # only amplified the GFP allele - only receiver allele
  # aux.Q224.QQ  = c(grep("CC", genotype, value=TRUE))
  # aux.Q224.Q   = setdiff(unique(c(grep("C", genotype, value=TRUE))),aux.Q224.QQ)
  # aux.L224.LL  = grep("DD", genotype, value=TRUE)
  # aux.L224.L   = setdiff(grep("D", genotype, value=TRUE),aux.L224.LL)
  # aux.NHEJ.EE  = grep("EE", genotype, value=TRUE)
  # aux.NHEJ.E   = setdiff(grep("E", genotype, value=TRUE),aux.NHEJ.EE)
  
  # only amplified the GFP allele - entire population
  aux.Q224.QQ  = c(grep("BB", genotype, value=TRUE),grep("BC", genotype, value=TRUE),grep("CC", genotype, value=TRUE))
  aux.Q224.Q   = setdiff(unique(c(grep("B", genotype, value=TRUE),grep("C", genotype, value=TRUE))),aux.Q224.QQ)
  aux.L224.LL  = grep("DD", genotype, value=TRUE)
  aux.L224.L   = setdiff(grep("D", genotype, value=TRUE),aux.L224.LL)
  aux.NHEJ.EE  = grep("EE", genotype, value=TRUE)
  aux.NHEJ.E   = setdiff(grep("E", genotype, value=TRUE),aux.NHEJ.EE)
  
  model$LD[1]       =   sum(ic$females[aux.LD]        + ic$males[aux.LD])/2
  model$GFPL[1]     =   sum(ic$females[aux.GFP]       + ic$males[aux.GFP])/2
  model$V9[1]       =   sum(ic$females[aux.Cas9]      + ic$males[aux.Cas9])/2
  model$LDonly[1]   =   sum(ic$females[aux.LDonly]    + ic$males[aux.LDonly])/2
  model$GFPLonly[1] =   sum(ic$females[aux.GFPLonly]  + ic$males[aux.GFPLonly])/2
  model$LDGFPL[1]   =   sum(ic$females[aux.LDGFPL]    + ic$males[aux.LDGFPL])/2
  model$Q224[1]     = ( sum(ic$females[aux.Q224.QQ]   + ic$males[aux.Q224.QQ]) + sum(ic$females[aux.Q224.Q] + ic$males[aux.Q224.Q])/2 )/2
  model$L224[1]     = ( sum(ic$females[aux.L224.LL]   + ic$males[aux.L224.LL]) + sum(ic$females[aux.L224.L] + ic$males[aux.L224.L])/2 )/2
  model$NHEJ[1]     = ( sum(ic$females[aux.NHEJ.EE]   + ic$males[aux.NHEJ.EE]) + sum(ic$females[aux.NHEJ.E] + ic$males[aux.NHEJ.E])/2 )/2
  
  # norm.factor = model$Q224[1]+model$L224[1]+model$NHEJ[1]
  norm.factor = 1 # allele frequency across the entire population

  model$Q224[1]     = model$Q224[1]/norm.factor
  model$L224[1]     = model$L224[1]/norm.factor
  model$NHEJ[1]     = model$NHEJ[1]/norm.factor
  
  ################################################################################
  # model experiment simulation ----
  ################################################################################
  
  for (i in 2:numGen) {
    
    # matting
    sim.F = tMatrix.F*array(outer(ic$females,ic$males), dim = c(length(ic$females),length(ic$males),length(ic$females)))
    sim.M = tMatrix.M*array(outer(ic$females,ic$males), dim = c(length(ic$females),length(ic$males),length(ic$males)))
    
    # maternal deposition
    sim.F = maternal_deposition(sim.F,genotype,parameters)
    sim.M = maternal_deposition(sim.M,genotype,parameters)
    
    # deterministic simulation
    det.F = apply(sim.F, 3, sum)
    det.M = apply(sim.M, 3, sum)
    
    # avoiding small negative numbers
    det.F[det.F<0] = 0
    det.M[det.M<0] = 0
    
    # normalizing values
    det.F = det.F/sum(det.F)
    det.M = det.M/sum(det.M)
    
    # stochastic simulation
    if (!is.null(sto)) {
      det.F = setNames(as.numeric(rmultinom(1,N.sample,prob = det.F)/N.sample), genotype)
      det.M = setNames(as.numeric(rmultinom(1,N.sample,prob = det.M)/N.sample), genotype)
    }
    
    # storing values
    model$LD[i]       =   sum(det.F[aux.LD]        + det.M[aux.LD])/2
    model$GFPL[i]     =   sum(det.F[aux.GFP]       + det.M[aux.GFP])/2
    model$V9[i]       =   sum(det.F[aux.Cas9]      + det.M[aux.Cas9])/2
    model$LDonly[i]   =   sum(det.F[aux.LDonly]    + det.M[aux.LDonly])/2
    model$GFPLonly[i] =   sum(det.F[aux.GFPLonly]  + det.M[aux.GFPLonly])/2
    model$LDGFPL[i]   =   sum(det.F[aux.LDGFPL]    + det.M[aux.LDGFPL])/2
    model$Q224[i]     = ( sum(det.F[aux.Q224.QQ]   + det.M[aux.Q224.QQ]) + sum(det.F[aux.Q224.Q] + det.M[aux.Q224.Q])/2 )/2
    model$L224[i]     = ( sum(det.F[aux.L224.LL]   + det.M[aux.L224.LL]) + sum(det.F[aux.L224.L] + det.M[aux.L224.L])/2 )/2
    model$NHEJ[i]     = ( sum(det.F[aux.NHEJ.EE]   + det.M[aux.NHEJ.EE]) + sum(det.F[aux.NHEJ.E] + det.M[aux.NHEJ.E])/2 )/2
    
    # norm.factor = model$Q224[i]+model$L224[i]+model$NHEJ[i] # only receiver allele
    norm.factor = 1 # allele frequency across the entire population

    model$Q224[i]     = model$Q224[i]/norm.factor
    model$L224[i]     = model$L224[i]/norm.factor
    model$NHEJ[i]     = model$NHEJ[i]/norm.factor
    
    # stochastic simulation - sampling
    if (!is.null(sto)) {
      model$V9[i]       = as.numeric(rbinom(1,data$data_screened_cage1[i],prob = model$V9[i])/data$data_screened_cage1[i])
      
      aux.sampling      = as.numeric(rmultinom(1,data$data_screened_cage1[i],prob = c(model$LDonly[i],model$GFPLonly[i],model$LDGFPL[i]))/data$data_screened_cage1[i])
      model$LDonly[i]   = aux.sampling[1]
      model$GFPLonly[i] = aux.sampling[2]
      model$LDGFPL[i]   = aux.sampling[3]
      
      aux.genotyping    = as.numeric(rmultinom(1,50, prob = c(model$Q224[i],model$L224[i],model$NHEJ[i])) /50)
      model$Q224[i]     = aux.genotyping[1]
      model$L224[i]     = aux.genotyping[2]
      model$NHEJ[i]     = aux.genotyping[3]
    }
    
    # new initial condition
    ic = list(
      females = det.F,
      males   = det.M
    )
  }
  
  return(model)
}
