###############################################################################
# maternal_deposition ----
###############################################################################

maternal_deposition<- function(tMatrix,genotype,parameters){
  
  # selecting individuals
  aux.A    = grep("A",    genotype, value=TRUE)
  aux.B    = grep("B",    genotype, value=TRUE)
  aux.AB   = grep("AB",   genotype, value=TRUE)
  aux.EE   = grep("EE",   genotype, value=TRUE)
  
  aux.ABC  = grep("ABC",  genotype, value=TRUE)
  aux.ABD  = grep("ABD",  genotype, value=TRUE)
  aux.ABE  = grep("ABE",  genotype, value=TRUE)
  aux.ACC  = grep("ACC",  genotype, value=TRUE)
  aux.ACD  = grep("ACD",  genotype, value=TRUE)
  aux.ACE  = grep("ACE",  genotype, value=TRUE)
  aux.ADD  = grep("ADD",  genotype, value=TRUE)
  aux.ADE  = grep("ADE",  genotype, value=TRUE)
  
  aux.aBC  = grep("aBC",  genotype, value=TRUE)
  aux.aBD  = grep("aBD",  genotype, value=TRUE)
  aux.aBE  = grep("aBE",  genotype, value=TRUE)
  aux.aCC  = grep("aCC",  genotype, value=TRUE)
  aux.aCD  = grep("aCD",  genotype, value=TRUE)
  aux.aCE  = grep("aCE",  genotype, value=TRUE)
  aux.aDD  = grep("aDD",  genotype, value=TRUE)
  aux.aDE  = grep("aDE",  genotype, value=TRUE)
  
  # parameters
  for (i in 1:length(parameters)) {
    assign(names(parameters)[i],as.numeric(parameters[i]))  
  }
  
  r.neutral = 1 - r.driving - r.NHEJ
  
  
  ##############################################################################
  # maternal deposition ----
  ##############################################################################
  
  # saving temporally
  temp.tMatrix = tMatrix
  
  # maternal deposition conversion GFP+L (BD) -> LD (BC)
  tMatrix[,,aux.ABC]           = tMatrix[,,aux.ABC]            + p.1* r.driving   *temp.tMatrix[,,aux.ABD]
  tMatrix[,,aux.ABD]           = tMatrix[,,aux.ABD]            - p.1* r.driving   *temp.tMatrix[,,aux.ABD]
  tMatrix[aux.A,,aux.aBC]      = tMatrix[aux.A,,aux.aBC]       + p.1* r.driving   *temp.tMatrix[aux.A,,aux.aBD]
  tMatrix[aux.A,,aux.aBD]      = tMatrix[aux.A,,aux.aBD]       - p.1* r.driving   *temp.tMatrix[aux.A,,aux.aBD]
  
  # maternal deposition conversion GFP+L (BD) -> NHEJ (BE)
  tMatrix[,,aux.ABE]           = tMatrix[,,aux.ABE]            + p.1* r.NHEJ      *temp.tMatrix[,,aux.ABD]
  tMatrix[,,aux.ABD]           = tMatrix[,,aux.ABD]            - p.1* r.NHEJ      *temp.tMatrix[,,aux.ABD]
  tMatrix[aux.A,,aux.aBE]      = tMatrix[aux.A,,aux.aBE]       + p.1* r.NHEJ      *temp.tMatrix[aux.A,,aux.aBD]
  tMatrix[aux.A,,aux.aBD]      = tMatrix[aux.A,,aux.aBD]       - p.1* r.NHEJ      *temp.tMatrix[aux.A,,aux.aBD]
  
  # maternal deposition conversion GFP+L (CD) -> GFP+Q (CC)
  tMatrix[aux.B,,aux.ACC]      = tMatrix[aux.B,,aux.ACC]      + p.1* r.driving   *temp.tMatrix[aux.B,,aux.ACD]
  tMatrix[aux.B,,aux.ACD]      = tMatrix[aux.B,,aux.ACD]      - p.1* r.driving   *temp.tMatrix[aux.B,,aux.ACD]
  tMatrix[aux.AB,,aux.aCC]     = tMatrix[aux.AB,,aux.aCC]     + p.1* r.driving   *temp.tMatrix[aux.AB,,aux.aCD]
  tMatrix[aux.AB,,aux.aCD]     = tMatrix[aux.AB,,aux.aCD]     - p.1* r.driving   *temp.tMatrix[aux.AB,,aux.aCD]
  
  # maternal deposition conversion GFP+L (CD) -> NHEJ (CE)
  tMatrix[aux.B,,aux.ACE]      = tMatrix[aux.B,,aux.ACE]      + p.1* r.NHEJ      *temp.tMatrix[aux.B,,aux.ACD]
  tMatrix[aux.B,,aux.ACD]      = tMatrix[aux.B,,aux.ACD]      - p.1* r.NHEJ      *temp.tMatrix[aux.B,,aux.ACD]
  tMatrix[aux.AB,,aux.aCE]     = tMatrix[aux.AB,,aux.aCE]     + p.1* r.NHEJ      *temp.tMatrix[aux.AB,,aux.aCD]
  tMatrix[aux.AB,,aux.aCD]     = tMatrix[aux.AB,,aux.aCD]     - p.1* r.NHEJ      *temp.tMatrix[aux.AB,,aux.aCD]
  
  # maternal deposition conversion GFP+L (DE) -> NHEJ (EE)
  tMatrix[aux.B,,aux.ADE]      = tMatrix[aux.B,,aux.ADE]      - p.1 *temp.tMatrix[aux.B,,aux.ADE]
  tMatrix[aux.AB,,aux.aDE]     = tMatrix[aux.AB,,aux.aDE]     - p.1 *temp.tMatrix[aux.AB,,aux.aDE]
  
  # maternal deposition conversion GFP+L (DD) -> NHEJ (EE)
  tMatrix[aux.B,,aux.ADD]      = tMatrix[aux.B,,aux.ADD]      - p.1 *temp.tMatrix[aux.B,,aux.ADD]
  tMatrix[aux.AB,,aux.aDD]     = tMatrix[aux.AB,,aux.aDD]     - p.1 *temp.tMatrix[aux.AB,,aux.aDD]
  
  ##############################################################################
  #  Unviable individuals ----
  ##############################################################################
  
  # NHEJ/NHEJ are assumed to be unviable
  tMatrix[,,aux.EE]  = 0
  
  return(tMatrix)
}
