if_zero_inflated = function(mat){
  zi = sapply(1:ncol(mat), function(j){
    #print(j)
    x = mat[,j]
    if(min(x) > 0) return(FALSE)
    fitobs = fitdist(x, "nbinom", method = "mme")
    bic0 = fitobs$aic
    
    a = fitobs$estimate["size"]
    b = fitobs$estimate["mu"]
    names(a) = names(b) = NULL
    start = list(k=a, lambda=b, omega=0.1)
    
    fitobs_zero = try(fitdist(x, "zinb", start = start, discrete = TRUE), silent = TRUE)
    ### !!!!! need further check
    if(class(fitobs_zero) == "try-error"){
      if(mean(x <= 1) > 0.5) return(TRUE)
      return(FALSE)
    }
    bic1 = fitobs_zero$aic
    if(fitobs_zero$estimate["omega"] > 0 & bic1 < bic0) return(TRUE)
    return(FALSE)
  })
  return(zi)
}

# mat: row-cells, column-genes
if_dropout_scimpute = function(mat, ncores, thre = 0.9){
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
  dropind = 1* (droprate > 0.9)
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


### covariance matrix - Quadrant
quadrant.transformed = function(x)
{
  x.m=apply(x,2,median)
  x=sweep(x,2,x.m)
  x.s=sign(x)
  # x.q=apply(x,2,Qn)
  x.q=apply(x,2,sd)
  cor.quadrant=cor(x.s)
  sigma.quadrant=diag(x.q)%*%cor.quadrant%*%diag(x.q)
  return(sigma.quadrant)
}

### covariance matrix - Spearman
spearman.transformed<-function(x)
{
  # x.q=apply(x,2,Qn)
  x.q=apply(x,2,sd)
  cor.sp=cor(x, method="spearman")
  sigma.sp=diag(x.q)%*%cor.sp%*%diag(x.q)
  return(sigma.sp)
}

### covariance matrix - Gaussian rank
Grank<-function(x)
{
  n=nrow(x)
  #x.q=apply(x,2,Qn)
  x.q=apply(x,2,sd)
  x.r=apply(x,2,rank)
  cor.Grank=cor(qnorm(x.r/(n+1)))
  sigma.gr=diag(x.q)%*%cor.Grank%*%diag(x.q)
  return(sigma.gr)
}


est_pearson.c = function(count, cor.p, thre_no = 10){
  count_dr = count
  count_dr[count_dr <= log10(2)] = NA
  cor_complete = cor(count_dr, use = "pairwise.complete.obs")
  cor_complete[is.na(cor_complete)] = cor.p[is.na(cor_complete)]
  
  count_dr = 1 * (!is.na(count_dr))
  nobs = t(count_dr) %*% count_dr
  cor_complete[nobs < thre_no] = cor.p[nobs < thre_no]
  return(cor_complete)
}

est_pearson_zinb.c = function(count, cor.p, zi, thre_no = 10){
  count_dr = count
  # count_dr[count_dr <= log10(2)] = NA
  for(j in which(zi)){
    count_dr[count_dr[,j]==0,j] = NA
  }
  cor_complete = cor(count_dr, use = "pairwise.complete.obs")
  cor_complete[is.na(cor_complete)] = cor.p[is.na(cor_complete)]
  
  count_dr = 1 * (!is.na(count_dr))
  nobs = t(count_dr) %*% count_dr
  cor_complete[nobs < thre_no] = cor.p[nobs < thre_no]
  return(cor_complete)
}


est_pearson_scimpute.c = function(count, cor.p, thre_no = 10, thre = 0.9){
  count_dr = count
  dropind = if_dropout_scimpute(count, ncores = 1, thre)
  count_dr[dropind == 1] = NA
  cor_complete = cor(count_dr, use = "pairwise.complete.obs")
  cor_complete[is.na(cor_complete)] = cor.p[is.na(cor_complete)]
  
  count_dr = 1 * (!is.na(count_dr))
  nobs = t(count_dr) %*% count_dr
  cor_complete[nobs < thre_no] = cor.p[nobs < thre_no]
  return(cor_complete)
}


est_weight_pearson_zinb = function(count, thre_no = 10){
  ### assuming log-sclae count
  count_dr = count
  zi = if_zero_inflated(10^count_dr-1)
  for(j in 1:ncol(count_dr)){
    if(zi[j]) count_dr[ count_dr[j]==0 , j] = NA
  }
  # count_dr[count_dr <= log10(2)] = NA
  cor_all = cor(count)
  cor_complete = cor(count_dr, use = "pairwise.complete.obs")
  cor_complete[is.na(cor_complete)] = cor_all[is.na(cor_complete)]
  
  count_dr = 1 * (!is.na(count_dr))
  nobs = t(count_dr) %*% count_dr
  cor_complete[nobs < thre_no] = cor_all[nobs < thre_no]
  weights = 1-abs(cor_complete)
  return(weights)
}

est_weight_pearson = function(count, thre_no = 10){
  ### assuming log-scale count
  count_dr = count
  count_dr[count_dr <= log10(2)] = NA
  cor_all = cor(count)
  cor_complete = cor(count_dr, use = "pairwise.complete.obs")
  cor_complete[is.na(cor_complete)] = cor_all[is.na(cor_complete)]
  
  count_dr = 1 * (!is.na(count_dr))
  nobs = t(count_dr) %*% count_dr
  cor_complete[nobs < thre_no] = cor_all[nobs < thre_no]
  weights = 1-abs(cor_complete)
  return(weights)
}

### use scImpute
est_weight_spearman = function(count, thre_no = 6){
  ### assuming log-sclae count
  count_dr = count
  #count_dr[count_dr == 0] = NA
  count_dr[count_dr <= log10(2)] = NA
  cor_all = cor(count)
  cor_complete = cor(count_dr, use = "pairwise.complete.obs", method = "spearman")
  cor_complete[is.na(cor_complete)] = cor_all[is.na(cor_complete)]
  
  count_dr = 1 * (!is.na(count_dr))
  nobs = t(count_dr) %*% count_dr
  cor_complete[nobs < thre_no] = cor_all[nobs < thre_no]
  weights = 1-abs(cor_complete)
  return(weights)
}


### use scImpute
est_weight_gaussian = function(count, thre_no = 6){
  ### assuming log-sclae count
  count_dr = count
  #count_dr[count_dr == 0] = NA
  count_dr[count_dr <= log10(2)] = NA
  
  n=nrow(count_dr)
  count_dr=apply(count_dr,2,rank, na.last = "keep")
  count_dr=qnorm(count_dr/(n+1))
  
  cor_all = cor(count)
  cor_complete = cor(count_dr, use = "pairwise.complete.obs", method = "spearman")
  cor_complete[is.na(cor_complete)] = cor_all[is.na(cor_complete)]
  
  count_dr = 1 * (!is.na(count_dr))
  nobs = t(count_dr) %*% count_dr
  cor_complete[nobs < thre_no] = cor_all[nobs < thre_no]
  weights = 1-abs(cor_complete)
  return(weights)
}


### use scImpute
est_weight_quadrant = function(count, thre_no = 6){
  ### assuming log-sclae count
  count_dr = count
  #count_dr[count_dr == 0] = NA
  count_dr[count_dr <= log10(2)] = NA
  
  x.m=apply(count_dr,2,median,na.rm = TRUE)
  count_dr=sweep(count_dr,2,x.m)
  count_dr=sign(count_dr)
  
  cor_all = cor(count)
  cor_complete = cor(count_dr, use = "pairwise.complete.obs", method = "spearman")
  cor_complete[is.na(cor_complete)] = cor_all[is.na(cor_complete)]
  
  count_dr = 1 * (!is.na(count_dr))
  nobs = t(count_dr) %*% count_dr
  cor_complete[nobs < thre_no] = cor_all[nobs < thre_no]
  weights = 1-abs(cor_complete)
  return(weights)
}


### glasso
get_res_GL = function(covobs, nobs, l1, adj_true){
  theta_list = QUIC(covobs, rho = 1, path = l1, msg = 0)
  
  accuracy = t(sapply(1:length(l1), function(i){
    theta = theta_list$X[,,i]
    theta[abs(theta) < 1e-5] = 0
    p = nrow(theta)
    
    loglik = -p * log(2*pi) + log(det(theta)) - sum(diag(covobs %*% theta))
    loglik = loglik * nobs * 0.5
    pn = sum(theta!=0 & (row(theta) > col(theta)))
    bic = log(nobs) * pn - 2*loglik
    
    ### precision and recall
    adj_est = 1 * (theta != 0)
    truepos = sum(adj_true == 1 & adj_est == 1) - p
    if(sum(adj_est == 1) == p){ precision = 1
    }else{ precision = truepos/(sum(adj_est == 1)-p) }
    recall = truepos/(sum(adj_true == 1)-p)
    ### FP rates
    falsepos = sum(adj_est == 1 & adj_true == 0)
    FP = falsepos/sum(adj_true == 0)
    
    metric = c(pn, bic, precision, recall, FP, l1[i])
    return(metric)
  }))
  return(list(theta_list = theta_list, accuracy = accuracy))
}


### glasso (QUIC) + scDesign weights
get_res_wGL_quic = function(covobs, weights, nobs, l1, adj_true){
  theta_list = lapply(l1, function(lda){
    res = QUIC(covobs, rho = lda * weights, msg = 0)
    return(res$X)
  })
  
  accuracy = t(sapply(1:length(theta_list), function(i){
    theta = theta_list[[i]]
    diag(theta) = 1
    theta[abs(theta) < 1e-5] = 0
    
    if(sum(is.na(theta)) > 0) return(rep(NA, 6))
    
    p = nrow(theta)
    loglik = -p * log(2*pi) + log(det(theta)) - sum(diag(covobs %*% theta))
    loglik = loglik * nobs * 0.5
    pn = sum(theta!=0 & (row(theta) > col(theta)))
    bic = log(nobs) * pn - 2*loglik
    
    ### precision and recall
    adj_est = 1 * (theta != 0)
    truepos = sum(adj_true == 1 & adj_est == 1) - p
    if(sum(adj_est == 1) == p){ precision = 1
    }else{ precision = truepos/(sum(adj_est == 1)-p) }
    recall = truepos/(sum(adj_true == 1)-p)
    ### FP rates
    falsepos = sum(adj_est == 1 & adj_true == 0)
    FP = falsepos/sum(adj_true == 0)
    
    metric = c(pn, bic, precision, recall, FP, l1[i])
    return(metric)
  }))
  return(list(theta_list = theta_list, accuracy = accuracy))
}




get_res_cors = function(count, adj_true, method = "pearson", cmat = NULL, thre_no = 6){
  if(is.null(cmat)){
    if(method == "pearson"){
      cmat = cor(count, method = "pearson")
    }
    if(method == "spearman"){
      cmat = cor(count, method = "spearman")
    }
    if(method == "gaussian"){
      n=nrow(count)
      x.r=apply(count,2,rank)
      cmat=cor(qnorm(x.r/(n+1)))
    }
    if(method == "quadrant"){
      x.m=apply(count,2,median)
      x=sweep(count,2,x.m)
      x.s=sign(x)
      cmat=cor(x.s)
    }
    
    if(method == "pearson.c"){
      count_dr = count
      #count_dr[count_dr == 0] = NA
      count_dr[count_dr <= log10(2)] = NA
      cor_all = cor(count)
      cor_complete = cor(count_dr, use = "pairwise.complete.obs")
      cor_complete[is.na(cor_complete)] = cor_all[is.na(cor_complete)]
      count_dr = 1 * (!is.na(count_dr))
      nobs = t(count_dr) %*% count_dr
      cor_complete[nobs < thre_no] = cor_all[nobs < thre_no]
      cmat = cor_complete
    }
    if(method == "spearman.c"){
      count_dr = count
      #count_dr[count_dr == 0] = NA
      count_dr[count_dr <= log10(2)] = NA
      cor_all = cor(count)
      cor_complete = cor(count_dr, use = "pairwise.complete.obs", method = "spearman")
      cor_complete[is.na(cor_complete)] = cor_all[is.na(cor_complete)]
      count_dr = 1 * (!is.na(count_dr))
      nobs = t(count_dr) %*% count_dr
      cor_complete[nobs < thre_no] = cor_all[nobs < thre_no]
      cmat = cor_complete
    }
    if(method == "gaussian.c"){
      count_dr = count
      count_dr[count_dr <= log10(2)] = NA
      
      n=nrow(count_dr)
      count_dr=apply(count_dr,2,rank, na.last = "keep")
      count_dr=qnorm(count_dr/(n+1))
      
      cor_all = cor(count)
      cor_complete = cor(count_dr, use = "pairwise.complete.obs", method = "spearman")
      cor_complete[is.na(cor_complete)] = cor_all[is.na(cor_complete)]
      
      count_dr = 1 * (!is.na(count_dr))
      nobs = t(count_dr) %*% count_dr
      cor_complete[nobs < thre_no] = cor_all[nobs < thre_no]
      cmat = cor_complete
    }
    if(method == "quadrant.c"){
      count_dr = count
      #count_dr[count_dr == 0] = NA
      count_dr[count_dr <= log10(2)] = NA
      
      x.m=apply(count_dr,2,median,na.rm = TRUE)
      count_dr=sweep(count_dr,2,x.m)
      count_dr=sign(count_dr)
      
      cor_all = cor(count)
      cor_complete = cor(count_dr, use = "pairwise.complete.obs", method = "spearman")
      cor_complete[is.na(cor_complete)] = cor_all[is.na(cor_complete)]
      
      count_dr = 1 * (!is.na(count_dr))
      nobs = t(count_dr) %*% count_dr
      cor_complete[nobs < thre_no] = cor_all[nobs < thre_no]
      cmat = cor_complete
    }
  }

  seqs = seq(1e-3,1-1e-3,length.out = 100)
  accuracy = t(sapply(seqs, function(thre){
    theta = cmat
    theta = (abs(theta) >= thre) * 1
    # theta = pvalmat
    # diag(theta) = 0
    # theta = (theta < alpha) * 1
    
    p = nrow(theta); nobs = nrow(count)
    pn = sum(theta!=0 & (row(theta) > col(theta)))
    ### precision and recall
    adj_est = theta
    truepos = sum(adj_true == 1 & adj_est == 1) - p
    if(sum(adj_est == 1) == p){ precision = 1
    }else{ precision = truepos/(sum(adj_est == 1)-p) }
    recall = truepos/(sum(adj_true == 1)-p)
    ### FP rates
    falsepos = sum(adj_est == 1 & adj_true == 0)
    FP = falsepos/sum(adj_true == 0)
    
    metric = c(pn, 0, precision, recall, FP, 0)
    return(metric)
  }))
  return(list(theta_list = cmat, accuracy = accuracy))
}

#path = paste0(rdir, "PIDC/weight-p20rho-0.2b-1.txt")
get_res_pidc = function(path, adj_true){
  cmat = read.table(path, header = FALSE)
  cmat = as.matrix(cmat)
  diag(cmat) = max(cmat)
  
  seqs = seq(1e-3,max(cmat)-1e-3,length.out = 100)
  accuracy = t(sapply(seqs, function(thre){
    theta = cmat
    theta = (abs(theta) >= thre) * 1
    
    p = nrow(theta); nobs = nrow(count)
    pn = sum(theta!=0 & (row(theta) > col(theta)))
    ### precision and recall
    adj_est = theta
    truepos = sum(adj_true == 1 & adj_est == 1) - p
    if(sum(adj_est == 1) == p){ precision = 1
    }else{ precision = truepos/(sum(adj_est == 1)-p) }
    recall = truepos/(sum(adj_true == 1)-p)
    ### FP rates
    falsepos = sum(adj_est == 1 & adj_true == 0)
    FP = falsepos/sum(adj_true == 0)
    
    metric = c(pn, 0, precision, recall, FP, 0)
    return(metric)
  }))
  return(list(theta_list = cmat, accuracy = accuracy))
}

get_res_genie3 = function(count, adj_true){
  cmat = GENIE3(t(count))
  diag(cmat) = max(cmat)
  seqs = seq(1e-3,max(cmat)-1e-3,length.out = 100)
  accuracy = t(sapply(seqs, function(thre){
    theta = cmat
    theta = (abs(theta) >= thre) * 1
    
    p = nrow(theta); nobs = nrow(count)
    pn = sum(theta!=0 & (row(theta) > col(theta)))
    ### precision and recall
    adj_est = theta
    truepos = sum(adj_true == 1 & adj_est == 1) - p
    if(sum(adj_est == 1) == p){ precision = 1
    }else{ precision = truepos/(sum(adj_est == 1)-p) }
    recall = truepos/(sum(adj_true == 1)-p)
    ### FP rates
    falsepos = sum(adj_est == 1 & adj_true == 0)
    FP = falsepos/sum(adj_true == 0)
    
    metric = c(pn, 0, precision, recall, FP, 0)
    return(metric)
  }))
  return(list(theta_list = cmat, accuracy = accuracy))
}











### Pearson's correlation
get_res_pearson = function(count, adj_true){
  #obs = rcorr(count, type = "pearson")
  #pvalmat = obs$P
  cmat = cor(count, method = "pearson")
  seqs = c(seq(1e-5,0.1,length.out = 8), seq(0.15,1,length.out = 8))
  
  accuracy = t(sapply(seqs, function(thre){
    theta = cmat
    theta = (abs(theta) >= thre) * 1
    # theta = pvalmat
    # diag(theta) = 0
    # theta = (theta < alpha) * 1
    
    p = nrow(theta); nobs = nrow(count)
    pn = sum(theta!=0 & (row(theta) > col(theta)))
    ### precision and recall
    adj_est = theta
    truepos = sum(adj_true == 1 & adj_est == 1) - p
    if(sum(adj_est == 1) == p){ precision = 1
    }else{ precision = truepos/(sum(adj_est == 1)-p) }
    recall = truepos/(sum(adj_true == 1)-p)
    ### FP rates
    falsepos = sum(adj_est == 1 & adj_true == 0)
    FP = falsepos/sum(adj_true == 0)
    
    metric = c(pn, 0, precision, recall, FP, 0)
    return(metric)
  }))
  return(list(theta_list = cmat, accuracy = accuracy))
}

### Spearman's correlation
get_res_spearman = function(count, adj_true){
  #obs = rcorr(count, type = "spearman")
  #pvalmat = obs$P
  cmat = cor(count, method = "spearman")
  seqs = c(seq(1e-5,0.12,length.out = 8), seq(0.15,1,length.out = 8))
  accuracy = t(sapply(seqs, function(thre){
    theta = cmat
    # diag(theta) = 0
    theta = (abs(theta) >= thre) * 1
    
    p = nrow(theta); nobs = nrow(count)
    pn = sum(theta!=0 & (row(theta) > col(theta)))
    ### precision and recall
    adj_est = theta
    truepos = sum(adj_true == 1 & adj_est == 1) - p
    if(sum(adj_est == 1) == p){ precision = 1
    }else{ precision = truepos/(sum(adj_est == 1)-p) }
    recall = truepos/(sum(adj_true == 1)-p)
    ### FP rates
    falsepos = sum(adj_est == 1 & adj_true == 0)
    FP = falsepos/sum(adj_true == 0)
    
    metric = c(pn, 0, precision, recall, FP, 0)
    return(metric)
  }))
  return(list(theta_list = cmat, accuracy = accuracy))
}

### Gaussian rank correlation
get_res_gaussian = function(count, adj_true){
  n=nrow(count)
  x.r=apply(count,2,rank)
  cmat=cor(qnorm(x.r/(n+1)))
  seqs = c(seq(1e-5,0.12,length.out = 8), seq(0.15,1,length.out = 8))
  accuracy = t(sapply(seqs, function(thre){
    theta = cmat
    # diag(theta) = 0
    theta = (abs(theta) >= thre) * 1
    
    p = nrow(theta); nobs = nrow(count)
    pn = sum(theta!=0 & (row(theta) > col(theta)))
    ### precision and recall
    adj_est = theta
    truepos = sum(adj_true == 1 & adj_est == 1) - p
    if(sum(adj_est == 1) == p){ precision = 1
    }else{ precision = truepos/(sum(adj_est == 1)-p) }
    recall = truepos/(sum(adj_true == 1)-p)
    ### FP rates
    falsepos = sum(adj_est == 1 & adj_true == 0)
    FP = falsepos/sum(adj_true == 0)
    
    metric = c(pn, 0, precision, recall, FP, 0)
    return(metric)
  }))
  return(list(theta_list = cmat, accuracy = accuracy))
}

### Quadrant correlation
get_res_quadrant = function(count, adj_true){
  x.m=apply(count,2,median)
  x=sweep(count,2,x.m)
  x.s=sign(x)
  cmat=cor(x.s)
  seqs = c(seq(1e-5,0.12,length.out = 8), seq(0.15,1,length.out = 8))
  accuracy = t(sapply(seqs, function(thre){
    theta = cmat
    # diag(theta) = 0
    theta = (abs(theta) >= thre) * 1
    
    p = nrow(theta); nobs = nrow(count)
    pn = sum(theta!=0 & (row(theta) > col(theta)))
    ### precision and recall
    adj_est = theta
    truepos = sum(adj_true == 1 & adj_est == 1) - p
    if(sum(adj_est == 1) == p){ precision = 1
    }else{ precision = truepos/(sum(adj_est == 1)-p) }
    recall = truepos/(sum(adj_true == 1)-p)
    ### FP rates
    falsepos = sum(adj_est == 1 & adj_true == 0)
    FP = falsepos/sum(adj_true == 0)
    
    metric = c(pn, 0, precision, recall, FP, 0)
    return(metric)
  }))
  return(list(theta_list = cmat, accuracy = accuracy))
}