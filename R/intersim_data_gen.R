#' Generate the simulated dataset with three modalities with the package InterSIM
#' @param prop proportion of samples for each cluster
#' @param n_sample the number of samples
#'
#' @return res, a list of length 2, where the first element is a list of simulated data, while the second element is a vector indicating the true label of each sample.
#' @importFrom InterSIM InterSIM
#' @export
#'
#' @examples
#' library(dplyr)
#' library(InterSIM)
#' intersim_data <- intersim_data_gen(prop=c(0.20,0.30,0.27,0.23), n_sample=500)


intersim_data_gen <- function(prop=c(0.20,0.30,0.27,0.23), n_sample=500){
  sim.data <- InterSIM(n.sample=n_sample, cluster.sample.prop = prop,
                       delta.methyl=5, delta.expr=5, delta.protein=5,
                       p.DMP=0.2, p.DEG=NULL, p.DEP=NULL,
                       sigma.methyl=NULL, sigma.expr=NULL, sigma.protein=NULL,
                       cor.methyl.expr=NULL, cor.expr.protein=NULL,
                       do.plot=FALSE, sample.cluster=TRUE, feature.cluster=TRUE)
  sim.methyl <- sim.data$dat.methyl
  sim.expr <- sim.data$dat.expr
  sim.protein <- sim.data$dat.protein
  temp_data <- list(sim.methyl, sim.expr, sim.protein)
  temp_id <- sim.data$clustering.assignment$cluster.id
  res <- list(temp_data,temp_id)
  return(res)
}
