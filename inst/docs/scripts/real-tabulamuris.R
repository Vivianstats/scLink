options(stringsAsFactors = FALSE)

ncores = 3

### setup ------------------------------------------------------------
library(parallel)
library(igraph)
library(matrixcalc)
library(MASS)
library(cglasso)
library(QUIC)


source("./get_mix_parameters.R")
source("./funs-get-res-realdata-wGL4.R")


rdir = "./"
dir.create(rdir, recursive = TRUE)

theta2pcor = function(theta){
  denom = sqrt(diag(theta))
  pcor = sweep(theta, MARGIN = 1, denom, FUN = "/")
  pcor = sweep(pcor, MARGIN = 2, denom, FUN = "/")
  return(-pcor)
}

### setup pars ------------------------------------------------------------
thre_edge_per = 0.05

### read and process data --------------------------------------------

data_by_ont_pro = readRDS(paste0(rdir, "data_by_ont_pro.rds"))
all_types = names(data_by_ont_pro)



### scLink computation
l1 = seq(1, 0.05, -0.05)
thre_no = 10
for(type in all_types){
  print(type)
  cnt = data_by_ont_pro[[type]]
  nc = ncol(cnt)
  res = est_graph(cnt, thre_no = thre_no, l1 = l1, ncores = ncores)
  saveRDS(res, file = paste0(rdir, "res_", type, ".rds"))
}

### extract results
p = 500
res_all = lapply(all_types, function(type){
  print(type)
  res = readRDS(paste0(rdir, "res_", type, ".rds"))
  is = which(sapply(1:length(l1), function(i){
    res$res_seq$accuracy[i,1] < thre_edge_per * p * (p-1) * 0.5
  }))
  if(length(is) == 0)return(NULL)
  istar = max(is)
  nedge = res$res_seq$accuracy[istar,1]
  theta = res$res_seq$theta_list[[istar]]
  pcor = theta2pcor(theta)
  pcor[abs(pcor) < 1e-5] = 0
  adjmat = 1 * (abs(theta) >= 1e-5)
  return(list(nedge = nedge, theta = theta, pcor = pcor, adjmat = adjmat))
})
names(res_all) = all_types
saveRDS(res_all, file = paste0(rdir, "res_all.rds"))
