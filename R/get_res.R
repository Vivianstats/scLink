

filter_genes_mean = function(x, topn = 500){
  mu = colMeans(x)
  ord = order(mu, decreasing = TRUE)[1:topn]
  genes = colnames(x)
  return(genes[ord])
}



if_dropout_scimpute = function(mat, ncores, dthre = 0.9){
  mat[mat <= log10(1.01)] = log10(1.01)
  pa = get_mix_parameters(t(mat), ncores = ncores)
  I = ncol(mat)
  J= nrow(mat)
  ### calculate cell-wise dropout rate
  droprate = sapply(1:I, function(i) {
    if(is.na(pa[i,1])) return(rep(0,J))
    wt = calculate_weight(mat[, i], pa[i, ])
    return(wt[, 1])
  })
  dropind = 1* (droprate > dthre)
  return(dropind)
}


easy.psd = function(sigma,method="perturb")
{
  if (method=="perturb")
  {
    p<-ncol(sigma)
    eig<-eigen(sigma, symmetric=T, only.values = T)
    const<-abs(min(eig$values,0))
    sigma.psd<-sigma+diag(p)*const
  }
  if (method=="npd")
  {
    eig<-eigen(sigma, symmetric=T)
    d<-pmax(eig$values,0)
    sigma.psd<-eig$vectors%*%diag(d)%*%t(eig$vectors)
  }
  return(sigma.psd)
}

est_pearson_scimpute.c = function(count, cor.p, thre_no = 20, 
                                  dthre = 0.9, ncores){
  count_dr = count
  x = log10(10^count - 1 + 1.01)
  dropind = if_dropout_scimpute(x, ncores = ncores, dthre)
  count_dr[dropind == 1] = NA
  cor_complete = cor(count_dr, use = "pairwise.complete.obs")
  cor_complete[is.na(cor_complete)] = cor.p[is.na(cor_complete)]
  
  count_dr = 1 * (!is.na(count_dr))
  nobs = t(count_dr) %*% count_dr
  cor_complete[nobs < thre_no] = cor.p[nobs < thre_no]
  return(cor_complete)
}

get_res_wGL_quic = function(covobs, weights, genes, nobs, l1, ncores){
  theta_list = mclapply(l1, function(lda){
    message(paste0("lambda = ", lda))
    res = QUIC(covobs, rho = lda * weights, msg = 0)
    X = res$X
    rownames(X) = colnames(X) = genes
    return(X)
  }, mc.cores = ncores)

  res = lapply(1:length(l1), function(i){
    Sigma = theta_list[[i]]
    adj = 1 * (Sigma != 0)
    
    p = nrow(Sigma)
    loglik = -p * log(2*pi) + log(det(Sigma)) - sum(diag(covobs %*% Sigma))
    loglik = loglik * nobs 
    nedge = (sum(adj)-p)/2
    bic = log(nobs) * nedge - loglik
    return(list(adj = adj, Sigma = Sigma, nedge = nedge, 
                bic = bic, lambda = l1[i]))
  })
  return(res)
}


est_graph = function(rcnt, thre_no = 8, l1, ncores = 1){
  cor_pearson = cor(rcnt)
  x.q=apply(rcnt,2,sd)
  ### calculate pairwise complete cor --------------------------------------------------
  cor_pearson.c = est_pearson_scimpute.c(rcnt, ncores = ncores, cor.p = cor_pearson, thre_no = thre_no)
  weights_p = 1-abs(cor_pearson.c)
  ### calculate pairwise complete cov --------------------------------------------------
  cov_pearson.c=diag(x.q)%*%cor_pearson.c%*%diag(x.q)
  cov_pearson.c = easy.psd(cov_pearson.c,method="perturb")
  ### wGL4 -----------------------------------------------------------------------------
  nobs = nrow(rcnt)
  genes = colnames(rcnt)
  res_seq = get_res_wGL_quic(covobs = cov_pearson.c, weights = weights_p, 
                             genes, nobs, l1 = l1, ncores = 1)
  return(list(res_seq=res_seq, cor_pearson.c = cor_pearson.c, cor_pearson = cor_pearson))
}




