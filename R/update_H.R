#' Update sub-matrices list Hi
#' @param WL a list of multiple modality data matrices
#' @param init_list a list of the initialized modality specific sub-matrices list Hi and shared sub-matrix E
#' @param lambda a parameter to set the relative weight of the L1,infinity norm defined on sub-matrices list Hi
#'
#' @return update_H_list, the data list init_list with the modality specific sub-matrices list Hi updated.
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
update_H <- function(WL,init_list)
{
  N1 <- length(init_list)-2
  E0 <- init_list[[(N1+1)]]
  dominators <- 0
  nominators <- 0
  for (b in 1:N1)
  {
    nominators <- t(E0)%*%WL[[b]]
    dominators <- t(E0)%*%E0
    dominators <- ginv(dominators)
    init_list[[b]] <- dominators%*%nominators
  }
  return(init_list)
}
