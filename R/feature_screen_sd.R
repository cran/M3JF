#' Screen the cluster related features via hypergeometric test p value and distribution standard derivation
#' @param feature_list a data list, which is the output of feature_selection function
#' @param sig_num the number of significant features for each cluster
#'
#' @return selected_features, a list the same long as the cluster number, each element is a sub-list with two vectors, one for the over-expressed features, one for the under-expressed features for the current cluster
#' @export
#' @importFrom stats sd
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
#' feature_list <- feature_selection(temp_data[[1]],M3JF_res$cluster_res,z_score=TRUE,
#' upper_bound=1, lower_bound=-1)
#' selected_features <- feature_screen_sd(feature_list,sig_num=20)
feature_screen_sd <- function(feature_list,sig_num=20){
  cluster_num <- length(feature_list)-2
  feature_num <- length(feature_list[[1]][[1]])
  selected_features <- vector('list',cluster_num)
  over_matr <- feature_list[[(cluster_num+1)]]
  under_matr <- feature_list[[(cluster_num+2)]]
  for (i in 1:cluster_num){
    this_feature_list <- vector('list',2)
    over_list <- c()
    under_list <- c()
    over_diff <- vector('numeric',feature_num)
    under_diff <- vector('numeric',feature_num)
    for (j in 1:feature_num){
      over_diff[j] <- sd(over_matr[j,])
      under_diff[j] <- sd(under_matr[j,])
    }
    this_feature_list[[1]] <- names(which(feature_list[[i]][[1]][sort(over_diff,decreasing = TRUE,index.return=TRUE)$ix]<0.05)[1:sig_num])
    this_feature_list[[2]] <- names(which(feature_list[[i]][[2]][sort(under_diff,decreasing = TRUE,index.return=TRUE)$ix]<0.05)[1:sig_num])
    selected_features[[i]] <- this_feature_list
  }
  return(selected_features)
}
