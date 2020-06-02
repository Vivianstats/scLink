library(ggplot2); theme_set(theme_classic())
library(ggrepel)

library(QUIC)
library(pheatmap)


library(clusterProfiler)
library(AnnotationHub)
library(quanteda)
library(msigdbr)

source("./get_mix_parameters.R")
source("./funs-get-res-realdata-wGL4.R")
rdir = "./"
thre_no = 10
ncores = 2

count = readRDS("./ESC/chu.rds")
labels = readRDS("./ESC/chu-label.rds")

count = t(count)
count = count[, colSums(count) > 0]
genes = colnames(count)
genes = strsplit(genes, split = "-")
genes = sapply(genes, function(x) x[1])
count = count[, genes != "ERCC"]
cs = rowSums(count)
scale.factor = 1e6
rcnt = sweep(count, 1, scale.factor/cs, FUN = "*")

### filter by genes -------------------------------------------
filter_genes_median = function(x, topn = 1500, frac_thre = 0.1){
  frac = colMeans(x > 0)
  x[x == 0] = NA
  mu = apply(x, 2, median, na.rm = TRUE)
  genes = colnames(x)
  mu = mu[frac > frac_thre]
  genes = genes[frac > frac_thre]
  
  ord = order(mu, decreasing = TRUE)[1:topn]
  return(genes[ord])
}
keep.genes = filter_genes_median(rcnt, topn = 1000, frac_thre = 0.1)
rcnt = rcnt[, keep.genes]
rcnt = log10(rcnt + 1)

### Estimate scLink's correlation -----------------------------

# for(t in times[-1]){
#   print(t)
#   l1 = seq(0.2, 0.01, -0.01)
#   res_seq = est_graph(rcnt[labels == t,], thre_no = 10, l1 = l1, ncores = ncores)
#   saveRDS(res_seq, file = paste0(rdir, "res_seq_", t, ".rds"))
# }
for(t in times){
  print(t)
  l1 = seq(0.2, 0.01, -0.01)
  res_seq = readRDS(paste0(rdir, "res_seq_", t, ".rds"))
  print(res_seq$res_seq$accuracy)
}

# 0.01 * 1000 * 999/2
res_seq = readRDS(paste0(rdir, "res_seq_", "00h", ".rds"))
theta = res_seq$res_seq$theta_list[[13]]
adj_start = 1 * (theta !=0)

res_seq = readRDS(paste0(rdir, "res_seq_", "96h", ".rds"))
theta = res_seq$res_seq$theta_list[[4]]
adj_end = 1 * (theta !=0)
