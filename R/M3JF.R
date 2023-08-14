#' Multi-Modal Matrix Joint Factorization
#' @param WL a list of multiple modality data matrices
#' @param lambda the parameter to set the relative weight of the group sparse constraint
#' @param theta threshold for the stopping criteria
#' @param k cluster number
#'
#' @return result, a list of 3 elements, the first element is a list comprising the
#' shared sub-matrix and the modality specific sub-matrices. The second element is
#' a vector of the clustering result. The third element is a vector of the cost in
#' each step during optimization.
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
#' M3JF_res <- M3JF(temp_data,k=4)
M3JF <- function(WL,lambda=0.01,theta=10^-6,k){
  if (!is.list(WL))
  {
    stop('Error:please provide a data list by WL')
  }
  if (is.na(k))
  {
    stop('Error:please provide a cluster number by k')
  }
  N <- length(WL)
  init_list <- initialize_WL(WL,k)
  divergence <- 1
  iter_num <- 1
  pre_cost <- 0
  cost_list <- c()
  while (divergence>theta)
  {
    update_H_list <- update_H(WL,init_list)
    update_E_list <- update_E(WL,update_H_list,lambda)
    init_list <- update_E_list
    new_cost <- cost(WL,init_list,lambda)
    divergence <- abs((new_cost-pre_cost)/new_cost)
    pre_cost <- new_cost
    cost_list <- c(cost_list,new_cost)
    message(iter_num)
    message(divergence)
    iter_num <- iter_num+1
    if (iter_num>500){
      break
    }
  }
  E <- init_list[[(N+1)]]
  clu_res <- kmeanspp(E,k)
  cluster_res <- clu_res$cluster
  result <- vector("list",3)
  names(result) <- c("sub_matrices","cluster_res",'cost')
  sub_result <- init_list
  names(sub_result) <- c(paste("sub_matrix_H",c(1:N),sep='_'),'E','E_loss')
  result[[1]] <- sub_result
  result[[2]] <- cluster_res
  result[[3]] <- cost_list
  return(result)
}
