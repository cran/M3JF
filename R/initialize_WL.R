#' Initialize the shared sub-matrix E and modality specific sub-matrices list Hi
#' @param WL a list of multiple modality data matrices
#' @param k the cluster number
#'
#' @return res, a list of length N+3, where N is the number of data modality. the first N elements are the modality specific sub-matrices list Hi, the (N+1) element is the shared sub-matrix E, the last two elements are the loss defined on the shared sub-matrix E and modality specific sub-matrices list Hi.
#' @export
#'
#' @examples
#' library(InterSIM)
#' sim.data <- InterSIM(n.sample=500, cluster.sample.prop = c(0.20,0.30,0.27,0.23),
#' delta.methyl=5, delta.expr=5, delta.protein=5,p.DMP=0.2, p.DEG=NULL,
#' p.DEP=NULL,sigma.methyl=NULL, sigma.expr=NULL, sigma.protein=NULL,cor.methyl.expr=NULL,
#' cor.expr.protein=NULL,do.plot=FALSE, sample.cluster=TRUE, feature.cluster=TRUE)
#' sim.methyl <- sim.data$dat.methyl
#' sim.expr <- sim.data$dat.expr
#' sim.protein <- sim.data$dat.protein
#' temp_data <- list(sim.methyl, sim.expr, sim.protein)
#' init_list <- initialize_WL(temp_data,k=4)
initialize_WL <- function(WL,k)
{
  NN <- length(WL)
  res <- vector("list",(NN+2))
  E0 <- matrix(0,nrow(WL[[1]]),k)
  for (i in 1:NN)
  {
    SVD_res <- svd(WL[[i]])
    dd <- SVD_res$d
    res[[i]] <- (diag(dd)%*%t(SVD_res$v))[2:(k+1),]
    E0 <- SVD_res$u[,2:(k+1)]
  }
  res[[(NN+1)]] <- E0
  return(res)
}
