#' Update sub-matrix E
#' @param WL a list of multiple modality data matrices
#' @param init_list a list of the initialized modality specific sub-matrices list Hi and shared sub-matrix E
#'
#' @return update_E_list, the data list init_list with the shared sub-matrix E updated.
#' @export
#' @importFrom MASS ginv
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
update_E <- function(WL,init_list,lambda)
{
  N2 <- length(init_list)-2
  E1 <- init_list[[(N2+1)]]
  new_E <- E1
  k <- ncol(E1)
  loss <- 0
  R_ij_list <- vector("list",N2)
  H_ij_list <- vector("list",N2)
  for (j in 1:k){
    for (i in 1:N2){
      R_ij <- WL[[i]]
      H_i <- init_list[[i]]
      for (l in 1:k){
        if (l!=j){
          R_ij <- R_ij-as.matrix(E1[,l],nrow(E1),1)%*%t(as.matrix(H_i[l,],1,ncol(H_i)))
        }
        else{
          H_ij <- t(as.matrix(H_i[j,],1,ncol(H_i)))
        }
      }
      R_ij_list[[i]] <- R_ij
      H_ij_list[[i]] <- H_ij
    }
    new_R_ij <- Reduce(cbind,R_ij_list)
    new_H_ij <- Reduce(cbind,H_ij_list)
    shrink <- (norm(new_H_ij,'2')^2)
    if (shrink<=.Machine$double.eps){
      loss <- loss+max(abs(new_E[,j]))
    }
    else{
      yy <- new_R_ij%*%t(new_H_ij)
      yy1 <- yy/(shrink)
      lambda1 <- lambda/shrink
      yy2 <- yy/lambda
      abs_y1 <- sort(abs(yy2),decreasing = T)
      count00 <- 0
      for (k1 in 1:length(abs_y1))
      {
        aa1<- (sum(abs_y1[1:k1])-1)/k1
        if (aa1<abs_y1[k1])
        {
          count00 <- k1
        }
      }
      a <- 1
      if (sum(abs_y1)<1)
      {
        a <- a+1
      }
      else
      {
        tao <- (sum(abs_y1[1:count00])-1)/count00
        for (k3 in 1:length(abs_y1))
        {
          if (yy2[k3]>=tao)
          {
            yy1[k3] <- tao*lambda1
          }
          else if (yy2[k3]<=(-tao))
          {
            yy1[k3] <- (-tao)*lambda1
          }
          else
          {
            yy1[k3] <- yy1[k3]
          }
        }

      }
      new_E[,j] <- yy1
      loss <- loss+max(abs(yy1))
    }
  }

  init_list[[(N2+1)]] <- new_E
  init_list[[(N2+2)]] <- loss
  return(init_list)
}
