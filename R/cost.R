#' Calculate the cost defined by the objective function
#' @param WL a list of multiple modality data matrices
#' @param init_list a list of the initialized modality specific sub-matrices list Hi and shared sub-matrix E
#' @param lambda a parameter to set the relative weight of the L1,infinity norm defined on sub-matrices list E
#'
#' @return res, a real data value of the cost
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
#' update_H_list <- update_H(temp_data,init_list)
#' lambda <- 0.01
#' update_E_list <- update_E(temp_data,update_H_list,lambda)
#' new_cost <- cost(temp_data,update_E_list,lambda)
cost <- function(WL,init_list,lambda)
{
  res <- 0
  NN <- length(WL)
  E1 <- init_list[[(NN+1)]]
  for (i in 1:NN)
  {
    matr_loss <- WL[[i]]-E1%*%init_list[[i]]
    res <- res+(norm(matr_loss,"F")^2)/2
  }
  res <- res+lambda*init_list[[(NN+2)]]
  return(res)
}
