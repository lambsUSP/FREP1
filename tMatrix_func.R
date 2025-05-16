###############################################################################
# tMatrix ----
###############################################################################

tMatrix_func <- function(genotype,parameters){
  
  # parameters
  for (i in 1:length(parameters)) {
    assign(names(parameters)[i],as.numeric(parameters[i]))  
  }
  
  # number of offspring (equal to male and females)
  num.offs = setNames(object = numeric(length(genotype)), nm = genotype)
  
  # Allelic drive experiments
  tMatrix.F = array(data=0, dim=c(length(genotype), length(genotype), length(genotype)), dimnames=list(genotype, genotype, genotype))
  tMatrix.M = array(data=0, dim=c(length(genotype), length(genotype), length(genotype)), dimnames=list(genotype, genotype, genotype))
  
  for (i in genotype) {
    for (j in genotype) {
      
      # female gametes
      aux.F = strsplit(i, split='')[[1]]
      gam.F = as.vector(outer(aux.F[1:2],aux.F[3:4], paste0, sep=''))
      
      # male gametes
      aux.M = strsplit(j, split='')[[1]]
      gam.M = as.vector(outer(aux.M[1:2],aux.M[3:4], paste0, sep=''))
      
      # offspring (outer product)
      offspring = as.vector(outer(gam.M,gam.F, paste0, sep=''))
      
      # reorder all offspring alleles to match allele order in gtype
      offspring = vapply(strsplit(offspring, split=''),
                         function(x) {paste0(sort(x, method = 'auto'), collapse='')},
                         FUN.VALUE = character(1))
      
      # count genotypes of offspring, order according to gtype
      for (o in offspring){num.offs[o] = num.offs[o]+1}
      
      # store after normalizing
      tMatrix.F[i,j, ] = num.offs/sum(num.offs)
      tMatrix.M[i,j, ] = num.offs/sum(num.offs)
      
      # resetting the count of offspring
      num.offs[] = 0
    }
  }
  
  ################################################################################
  # Fitness cost ----
  ################################################################################
  
  # aux.BB  = grep("BB",  genotype, value=TRUE)
  # tMatrix.F[,,aux.BB] = f.BB*tMatrix.F[,,aux.BB] # LD/LD
  # tMatrix.M[,,aux.BB] = f.BB*tMatrix.M[,,aux.BB] # LD/LD
  
  ################################################################################
  # Drive action ----
  ################################################################################
  
  # Condition for drive action: Cas9 (A) + gRNA (B) + CFPL (D)
  # Drive action: GFP+L (D) into GFP+Q (C) or NHEJ (E) ----
  r.neutral = 1 - r.driving - r.NHEJ
  
  # saving temporally
  temp.tMatrix.F = tMatrix.F
  temp.tMatrix.M = tMatrix.M
  
  # females ----
  # aux.ABD[1]: aABD
  tMatrix.F[,,c("aABC")] = tMatrix.F[,,c("aABC")] + r.driving     * temp.tMatrix.F[,,c("aABD")] # GFP+L->GFP+Q    
  tMatrix.F[,,c("aABD")] = tMatrix.F[,,c("aABD")] - r.driving     * temp.tMatrix.F[,,c("aABD")]
  tMatrix.F[,,c("aABE")] = tMatrix.F[,,c("aABE")] + r.NHEJ        * temp.tMatrix.F[,,c("aABD")] # GFP+L->NHEJ    
  tMatrix.F[,,c("aABD")] = tMatrix.F[,,c("aABD")] - r.NHEJ        * temp.tMatrix.F[,,c("aABD")]
  
  # aux.ABD[2]: AABD 
  tMatrix.F[,,c("AABC")] = tMatrix.F[,,c("AABC")] + r.driving     * temp.tMatrix.F[,,c("AABD")] # GFP+L->GFP+Q    
  tMatrix.F[,,c("AABD")] = tMatrix.F[,,c("AABD")] - r.driving     * temp.tMatrix.F[,,c("AABD")]
  tMatrix.F[,,c("AABE")] = tMatrix.F[,,c("AABE")] + r.NHEJ        * temp.tMatrix.F[,,c("AABD")] # GFP+L->NHEJ    
  tMatrix.F[,,c("AABD")] = tMatrix.F[,,c("AABD")] - r.NHEJ        * temp.tMatrix.F[,,c("AABD")]
  
  
  # males ----
  # aux.ABD[1]: aABD
  tMatrix.M[,,c("aABC")] = tMatrix.M[,,c("aABC")] + r.driving     * temp.tMatrix.M[,,c("aABD")] # GFP+L->GFP+Q    
  tMatrix.M[,,c("aABD")] = tMatrix.M[,,c("aABD")] - r.driving     * temp.tMatrix.M[,,c("aABD")]
  tMatrix.M[,,c("aABE")] = tMatrix.M[,,c("aABE")] + r.NHEJ        * temp.tMatrix.M[,,c("aABD")] # GFP+L->NHEJ    
  tMatrix.M[,,c("aABD")] = tMatrix.M[,,c("aABD")] - r.NHEJ        * temp.tMatrix.M[,,c("aABD")]
  
  # aux.ABD[2]: AABD 
  tMatrix.M[,,c("AABC")] = tMatrix.M[,,c("AABC")] + r.driving     * temp.tMatrix.M[,,c("AABD")] # GFP+L->GFP+Q    
  tMatrix.M[,,c("AABD")] = tMatrix.M[,,c("AABD")] - r.driving     * temp.tMatrix.M[,,c("AABD")]
  tMatrix.M[,,c("AABE")] = tMatrix.M[,,c("AABE")] + r.NHEJ        * temp.tMatrix.M[,,c("AABD")] # GFP+L->NHEJ    
  tMatrix.M[,,c("AABD")] = tMatrix.M[,,c("AABD")] - r.NHEJ        * temp.tMatrix.M[,,c("AABD")]
  
  ################################################################################
  #  Unviable individuals ----
  ################################################################################  
  
  # NHEJ/NHEJ are assumed to be unviable
  aux.EE   = grep("EE",   genotype, value=TRUE)
  
  tMatrix.F[,,aux.EE]  = 0
  tMatrix.M[,,aux.EE]  = 0
  
  ##############################################################################
  
  return(list(tMatrix.F,tMatrix.M))
}
