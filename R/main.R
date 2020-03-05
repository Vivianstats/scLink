#' Pre-process data for scLink
#' @param count A full gene count matrix with rows representing cells and columns representing genes.
#' Gene names are given as column names.
#' @param scale.factor A number specifying the sclae factor used for library size normalization.
#' Defaults to 1e6.
#' @param filter.genes A Boolean specifying whether scLink should select genes based on mean expression.
#' When set to \code{FALSE}, users need to speicfy a set of genes to be used for network construction with
#' \code{gene.names}. When set to \code{TRUE}, scLink will select genes based on their mean expression,
#' and users need to specify the number of genes to be selected with \code{n}.
#' @param gene.names A character vector specifying the genes used for network construction.
#' Only needed when \code{filter.genes = FALSE}.
#' @param n An integer specifying the number of genes to be selected by scLink (defaults to 500).
#' Only needed when \code{filter.genes = TRUE}.
#' @return A transformed and normalized gene expression matrix that can be used for correlation calculation
#' and network construction.
#' @export
#' @author Wei Vivian Li, \email{vivian.li@rutgers.edu}
sclink_norm = function(count, scale.factor = 1e6, filter.genes = FALSE, gene.names = NULL, n = 500){
  if(filter.genes == FALSE){
    if(is.null(gene.names)) {stop("Need to supply gene.names when setting filter.genes to FALSE!")}
  }
  if(filter.genes == TRUE){
    message("Setting filter.genes = TRUE! \n
             scLink will select the top n expressed genes!")
  }
  count = as.matrix(count)
  count = sweep( count, MARGIN = 1, scale.factor/rowSums( count), FUN = "*")
  count = log10(count + 1)
  if(filter.genes == FALSE){
    m = sum(gene.names %in% colnames(count))
    if(m < 3){warning("Less than 3 genes in gene.names match the column names!")}
    count = count[, gene.names, drop = FALSE]
  }else{
    keep.genes = filter_genes_mean(count, topn = n)
    count = count[, keep.genes]
  }
  return(count)
}


#' Calculate scLink's correlation matrix
#' @param expr A gene expression matrix with rows representing cells and columns representing genes.
#' Gene names are given as column names. Can be the output of \code{sclink_norm} or user constructed gene
#' expression matrices.
#' @param ncores Number of cores if using parallel computation.
#' @param nthre An integer specifying a threshold on the number of complete observations.
#' Defaults to 20.
#' @param dthre A number specifying the threshold on dropout probabilities. Defaults to 0.9.
#' @return A correlation matrix for gene co-expression relationships.
#' @export
#' @importFrom parallel mclapply
#' @importFrom stats cor dgamma dnorm pchisq sd uniroot
#' @author Wei Vivian Li, \email{vivian.li@rutgers.edu}
sclink_cor = function(expr, ncores, nthre = 20, dthre = 0.9){
  cor_pearson = cor(expr)
  x.q = apply(expr, 2, sd)
  cor_pearson.c = est_pearson_scimpute.c(expr, ncores = ncores, cor.p = cor_pearson,
                                         thre_no = nthre, dthre = dthre)
  return(cor_pearson.c)
}


#' Infer gene co-expression networks
#' @param expr A gene expression matrix with rows representing cells and columns representing genes.
#' Gene names are given as column names. Can be the output of \code{sclink_norm} or user constructed gene
#' expression matrices.
#' @param ncores Number of cores if using parallel computation.
#' @param lda A vector specifying a sequence of lambda values to be used in the penalized likelihood.
#' @param nthre An integer specifying a threshold on the number of complete observations.
#' Defaults to 20.
#' @param dthre A number specifying the threshold on dropout probabilities. Defaults to 0.9.
#' @return A list for gene co-expression relationships. The list contains a \code{cor} element for
#' scLink's correlation matrix and a \code{summary} element for the gene networks. \code{summary} is a list
#' with each element corresponding to the result of one lambda value. Each element of \code{summary}
#' contains the following information:
#' \describe{
#'   \item{adj:}{the adjacency matrix specifying the gene-gene edges;}
#'   \item{Sigma:}{the estimated concentration matrix;}
#'   \item{nedge:}{number of edges in the gene network;}
#'   \item{bic:}{BIC score;}
#'   \item{lambda:}{value of lambda in the penalty.}
#' }
#' @export
#' @importFrom QUIC QUIC
#' @importFrom parallel mclapply
#' @importFrom stats cor dgamma dnorm pchisq sd uniroot
#' @author Wei Vivian Li, \email{vivian.li@rutgers.edu}
sclink_net = function(expr, ncores, lda = seq(1, 0.1, -0.05),
                      nthre = 20, dthre = 0.9){
  if(ncol(expr) > 1000){warning("For large expression matrix, \n
                                 computation may take a long time. \n
                                 Consider first using the scLink_cor function to check correlation structures.")}
  cor_pearson = cor(expr)
  x.q = apply(expr, 2, sd)
  cor_pearson.c = est_pearson_scimpute.c(expr, ncores = ncores, cor.p = cor_pearson,
                                         thre_no = nthre, dthre = dthre)

  weights_p = 1 - abs(cor_pearson.c)
  cov_pearson.c= diag(x.q) %*% cor_pearson.c %*% diag(x.q)
  cov_pearson.c = easy.psd(cov_pearson.c,method="perturb")

  nobs = nrow(expr)
  genes = colnames(expr)
  message("Constructing gene networks ...")
  res_seq = get_res_wGL_quic(covobs = cov_pearson.c, weights = weights_p,
                             genes, nobs, l1 = lda, ncores = ncores)

  return(list(summary = res_seq, cor = cor_pearson.c))
}




