#' Select the cluster related features via hypergeometric test
#' @param X the feature matrix to be analyzed, with rows as samples and columns as features
#' @param clusters the numeric cluster results with number specifying the cluster
#' @param z_score a binary value to specify whether to calculate z-score for X first
#' @param upper_bound values larger than this value should be treated as over-expressed
#' @param lower_bound values smaller than this value should be treated as under-expressed
#' @param p.adjust.method the p value adjust method, defalut as 'BH'
#'
#' @return results, a list, which is as long as (cluster number+2), with the first (cluster number) element
#' as two sub-list, each composing a feature vector and a FDR vector. The last two elements are two matrices,
#' one is the matrix representing the fraction of over-express samples in each cluster for each features , and
#' the other represents that of under-express.
#' @export
#' @importFrom stats phyper p.adjust
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
feature_selection <- function(X, clusters, z_score=FALSE, upper_bound, lower_bound, p.adjust.method='BH'){
  sample_num <- nrow(X)
  feature_num <- ncol(X)
  uni_cluster <- unique(clusters)
  cluster_num <- length(uni_cluster)
  if (sample_num!=length(clusters)){
    stop('Error:rows of X should be the same as the length of clusters')
  }
  if (z_score){
    X <- apply(X,2,scale)
  }
  results <- vector('list',(cluster_num+2))
  for (i in 1:cluster_num){
    message(i)
    this_result <- vector('list',2)
    # over_result <- vector('list',2)
    over_matrix <- matrix(0,feature_num,6)
    over_matrix[,1] <- colnames(X)
    over_matrix[,2] <- sample_num
    this_index <- which(clusters==i,TRUE)
    this_cluster_num <- length(this_index)
    over_matrix[,3] <- this_cluster_num
    # binarize for over-expressed test
    new_X <- X
    new_X[which(new_X<upper_bound,TRUE)] <- 0
    new_X[which(new_X>=upper_bound,TRUE)] <- 1
    # hypergenomic test
    for (j in 1:feature_num){
      total_over_num <- length(which(new_X[,j]==1,TRUE))
      over_matrix[j,4] <- total_over_num
      this_over_num <- length(which(new_X[this_index,j]==1,TRUE))
      over_matrix[j,5] <- this_over_num
      over_matrix[j,6] <- phyper((as.numeric(over_matrix[j,5])-1),as.numeric(over_matrix[j,4]),(as.numeric(over_matrix[j,2])-as.numeric(over_matrix[j,4])),as.numeric(over_matrix[j,3]),lower.tail = FALSE)
    }
    over_p_values <- as.numeric(over_matrix[,6])
    names(over_p_values) <- over_matrix[,1]
    over_p_adjust <- p.adjust(over_p_values,method = p.adjust.method)
    # over_result[[1]] <- over_p_adjust


    # under_result <- vector('list',2)
    under_matrix <- matrix(0,feature_num,6)
    under_matrix[,1] <- colnames(X)
    under_matrix[,2] <- sample_num
    this_index <- which(clusters==i,TRUE)
    this_cluster_num <- length(this_index)
    under_matrix[,3] <- this_cluster_num
    # binarize for over-expressed test
    new_X <- X
    new_X[which(new_X>lower_bound,TRUE)] <- 0
    new_X[which(new_X<=lower_bound,TRUE)] <- 1
    # hypergenomic test
    for (j in 1:feature_num){
      total_under_num <- length(which(new_X[,j]==1,TRUE))
      under_matrix[j,4] <- total_under_num
      this_under_num <- length(which(new_X[this_index,j]==1,TRUE))
      under_matrix[j,5] <- this_under_num
      under_matrix[j,6] <- phyper((as.numeric(under_matrix[j,5])-1),as.numeric(under_matrix[j,4]),(as.numeric(under_matrix[j,2])-as.numeric(under_matrix[j,4])),as.numeric(under_matrix[j,3]),lower.tail = FALSE)
    }
    under_p_values <- as.numeric(under_matrix[,6])
    names(under_p_values) <- under_matrix[,1]
    under_p_adjust <- p.adjust(under_p_values,method = p.adjust.method)
    # under_result[[1]] <- under_p_adjust

    this_result[[1]] <- over_p_adjust
    this_result[[2]] <- under_p_adjust
    results[[i]] <- this_result
  }
  over_fraction <- matrix(0,feature_num,cluster_num)
  new_X <- X
  new_X[which(new_X<upper_bound,TRUE)] <- 0
  new_X[which(new_X>=upper_bound,TRUE)] <- 1
  for (j in 1:feature_num){
    for (i in 1:cluster_num){
      this_index <- which(clusters==i,TRUE)
      this_total <- length(this_index)
      this_over <- length(which(new_X[this_index,j]==1,TRUE))
      over_fraction[j,i] <- this_over/this_total
    }
  }
  results[[(cluster_num+1)]] <- over_fraction

  under_fraction <- matrix(0,feature_num,cluster_num)
  new_X <- X
  new_X[which(new_X>lower_bound,TRUE)] <- 0
  new_X[which(new_X<=lower_bound,TRUE)] <- 1
  for (j in 1:feature_num){
    for (i in 1:cluster_num){
      this_index <- which(clusters==i,TRUE)
      this_total <- length(this_index)
      this_under <- length(which(new_X[this_index,j]==1,TRUE))
      under_fraction[j,i] <- this_under/this_total
    }
  }
  results[[(cluster_num+2)]] <- under_fraction
  return(results)
}
