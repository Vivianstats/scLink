
ncores = 3

rdir = "./"
dir.create(rdir, recursive = TRUE)


library(igraph)

library(pheatmap)
library(gridExtra)
library(grid)

library(RColorBrewer)
library(colorspace)
library(dplyr)

library(clusterProfiler)
library(AnnotationHub)
library(quanteda)
library(msigdbr)

library(QUIC)

source("./get_mix_parameters.R")
source("./funs-get-res-realdata-wGL4.R")

m2v = function(mat){mat[upper.tri(mat)]}
theta2pcor = function(theta){
  denom = sqrt(diag(theta))
  pcor = sweep(theta, MARGIN = 1, denom, FUN = "/")
  pcor = sweep(pcor, MARGIN = 2, denom, FUN = "/")
  return(-pcor)
}


### ------------------------------------------------------------------------------
### read data --------------------------------------------------------------------
### ------------------------------------------------------------------------------

data  = readRDS('./Azizi2018/countByPatient.rds')
data2 = readRDS('./Azizi2018/countByPatient-NORMAL.rds')
count.patient = data$BC3$count
count.normal = data2$BC3$count
dim(count.patient)
dim(count.normal)
universe.genes = rownames(count.patient)


rcnt = count.patient
thre_no = 10
thre_edge_per = 0.05
scale.factor = 1e5


### ------------------------------------------------------------------------------
### normal cells --------------------------------------------------------------------
### ------------------------------------------------------------------------------

### get normalized count matrix
rcnt = t(count.normal)
cs = rowSums(rcnt)
rcnt = sweep(rcnt, 1, scale.factor/cs, FUN = "*")
### filter by genes
keep.genes = filter_genes_mean(rcnt, topn = 500)
# genes2 = filter_genes_mean(rcnt, topn = 300)
cnt.normal = rcnt[, keep.genes]
cnt.normal = log10(cnt.normal + 1)

l1 = seq(1.2, 0.5, -0.1)
res_seq_normal = est_graph(cnt.normal, thre_no = 10, l1 = l1, ncores = ncores)
p = 500
istar = max(which(sapply(1:length(l1), function(i){
  res_seq_normal$res_seq$accuracy[i,1] < thre_edge_per * p * (p-1) * 0.5
})))
theta.normal = res_seq_normal$res_seq$theta_list[[istar]]
pcor.normal = theta2pcor(theta.normal)
adjmat.normal = 1 * (abs(theta.normal) >= 1e-5)


### ------------------------------------------------------------------------------
### tumor cells --------------------------------------------------------------------
### ------------------------------------------------------------------------------

rcnt = t(count.patient)
cs = rowSums(rcnt)
rcnt = sweep(rcnt, 1, scale.factor/cs, FUN = "*")
cnt.pt = rcnt[, keep.genes]
cnt.pt = log10(cnt.pt + 1)

l1 = seq(1.2, 0.5, -0.1)
res_seq_pt = est_graph(cnt.pt, thre_no = 10, l1 = l1, ncores = ncores)
p = 500
istar = max(which(sapply(1:length(l1), function(i){
  res_seq_pt$res_seq$accuracy[i,1] < thre_edge_per * p * (p-1) * 0.5
})))

theta.pt = res_seq_pt$res_seq$theta_list[[istar]]
pcor.pt = theta2pcor(theta.pt)
adjmat.pt = 1 * (abs(theta.pt) >= 1e-5)

