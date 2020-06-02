controlsPcorSimulator <- function(nobs, nclusters, nnodesxcluster, pattern, 
                           low.strength, sup.strength, nhubs, degree.hubs, nOtherEdges, 
                           alpha, plus, prob, perturb.clust, mu, probSign, seed)
{
  ## prospective patterns
  ppaterns <- c("hubs", "powerLaw", "random")
  if(length(pattern)>1) 
   warning("pattern attribute length is larger than one. Only the first component will be used")
  pattern <- ppaterns[pmatch(pattern[1], ppaterns)]
  if (is.na(pattern)) stop("pattern is not well define. It must be selected from \"hubs\", \"powerLaw\" or \"random\" ")
    
  ## nobs, nclusters and nnodesxclustern  
  if(nobs < 10) stop("nobs must be larger 10")
  if(nclusters < 1) stop("nclusters is not well defined. At least one cluster is needed") 
  if(length(nnodesxcluster)!= nclusters) 
    stop("nnodesxcluster must be a vector of length defined by nclusters")
  
  ## low.strength and sup.strength
  if(low.strength < 0 | low.strength > 1) stop("low.strength must be between 0 and 1")
  if(sup.strength < 0 | sup.strength > 1) stop("sup.strength must be between 0 and 1")
  if(sup.strength < low.strength) stop("sup.strength must be larger than low.strength")
 
  ## nhubs, degree.hubs, nOtherEdges 
  if(pattern == "hubs"){
    if(length(nhubs)==1)  nhubs <- rep(nhubs,nclusters)  
    if(length(nhubs) != nclusters)  stop("nhubs must be a vector of size 1 or nclusters")
    if(any(nhubs < 0))  stop("nhubs must be larger than zero")
    if(any(nhubs - nnodesxcluster >= 0) ) stop("nhubs must contain smaller values than the cluster sizes given by attribute nnodesxcluster")
  
    if(length(degree.hubs)==1)  degree.hubs <- rep(degree.hubs,nclusters)  
    if(length(degree.hubs) != nclusters)  stop("degree.hubs must be a vector of size 1 or nclusters")
    if(any(degree.hubs <= 1))  stop("degree.hubs must be larger than one")
    if(any(degree.hubs - nnodesxcluster >= 0))  stop("degree.hubs must contain smaller values than the cluster sizes given by attribute nnodesxcluster")
  
    if(length(nOtherEdges)==1)  nOtherEdges <- rep(nOtherEdges,nclusters)  
    if(length(nOtherEdges) != nclusters)  stop("nOtherEdges must be a vector of size 1 or nclusters")
    if(any(nOtherEdges < 0))  stop("nOtherEdges must be larger than zero")
    if(any(nOtherEdges - nnodesxcluster*(nnodesxcluster-1)/2 >= 0))  stop("too many nOtherEdges are selected. Reduce the number")
  }
  
  ## alpha, plus
  if(pattern == "powerLaw"){
    if(alpha <= 0) stop("alpha must be positive")
    if(plus < 0) stop("plus must be positive")
  }
  
  ## prob
  if(pattern == "random"){
    if(length(prob) == 1)  prob <- rep(prob,nclusters)  
    if(length(prob) != nclusters)  stop("prob must be a vector of size 1 or nclusters")
    if(any(prob < 0 | prob >1))  stop("prob must be between zero and one")
  }

  ## perturb.clust, mu, prob.sign, seed
  if(length(perturb.clust) > 1)  perturb.clust <- perturb.clust[1]
  if(perturb.clust < 0 | perturb.clust >1)  stop("perturb.clust must be between zero and one")
  
  if(length(mu) == 1) mu <- rep(mu, sum(nnodesxcluster))
  if(length(mu) != sum(nnodesxcluster)) stop("length of mu must coincide with the sum of nnodesxcluster")
 
  if(length(probSign) > 1)  probSign <- probSign[1]
  if(probSign < 0 | probSign >1)  stop("probSign must be between zero and one")
  
  if(length(seed) < nclusters) stop("seed must be a vector of at least size nclusters")
  
  listAll <- list(nobs = nobs, nclusters = nclusters, nnodesxcluster = nnodesxcluster,
  			      pattern = pattern, low.strength = low.strength, sup.strength = sup.strength, 
  			      nhubs = nhubs, degree.hubs = degree.hubs, nOtherEdges = nOtherEdges, 
                  alpha = alpha, plus = plus, prob = prob, perturb.clust = perturb.clust, 
                  mu = mu, probSign = probSign, seed = seed) 
  return(listAll)
}
