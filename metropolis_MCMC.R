###############################################################################
# MCMC algorithm ---
###############################################################################

metropolis_MCMC <- function(data, numGen, genotype, startvalue, iterations, scale){
  
  # chain and its initial value (last column stores the likelihood value)
  chain      = array(dim = c(iterations+1,length(startvalue)+1))
  parameters = startvalue
  chain[1,]  = unlist(c(startvalue,ll=obj_func(data,numGen,genotype,parameters)))
  
  # progress bar
  pb = txtProgressBar(min = 0, max = iterations, style = 1, char="=")
  
  for (i in 1:iterations){
    # sample proposal from a normal distribution considering the informed scale
    proposal = pmax(0,pmin(2,rnorm(length(startvalue),
                     mean = chain[i,1:length(startvalue)], 
                     sd   = scale)))
    
    # log-likelihood of currently parameter set
    parameters[] = chain[i,1:length(startvalue)]
    ll = obj_func(data,numGen,genotype,parameters)
    
    # log-likelihood of new proposal
    parameters[] = proposal
    ll.proposal = obj_func(data,numGen,genotype,parameters)
    
    # print(round(c(i,proposal,ll.proposal,ll),2))
    
    # select or not new proposal according to the Metropolis alg rules
    probab = exp(ll.proposal - ll)
    if (runif(1) < probab){chain[i+1,] = c(proposal,ll.proposal)
    }else{chain[i+1,] = c(chain[i,])}
    
    # progress bar
    setTxtProgressBar(pb, value = i)
  }
  
  # progress bar
  close(pb)
  
  # naming the last column of loglikelihood values
  colnames(chain) = c(names(parameters),"ll")
  
  return(chain)
}

