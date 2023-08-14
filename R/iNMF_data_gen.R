#' Generate the simulated dataset with three modalities as the work iNMF
#' @param Xs_dim_list a list of data matrix dimensions for multiple modality data
#' @param mod_dim_list a list of the dimensions of each cluster and their features
#' @param e_u the level of uniform noise
#' @param e_s signal to noise ratio
#' @param e_h block adding probability
#'
#' @return res, a list of length 2, where the first element is a list of simulated data, while the second element is a vector indicating the true label of each sample.
#' @importFrom stats rbinom rbeta runif
#' @export
#'
#' @examples
#' library(dplyr)
#' iNMF_data <- iNMF_data_gen(Xs_dim_list=list(c(100,100),c(100,100),c(100,100)),mod_dim_list=list(matrix(c(20,30,20,30,20,30,20,30),4,2),matrix(c(20,20,30,30,20,30,20,30),4,2),matrix(c(26,24,26,24,20,30,20,30),4,2)),e_u=0.15, e_s=0.9, e_h=0)

iNMF_data_gen <- function(Xs_dim_list=list(c(100,100),c(100,100),c(100,100)),mod_dim_list=list(matrix(c(20,30,20,30,20,30,20,30),4,2),matrix(c(20,20,30,30,20,30,20,30),4,2),matrix(c(26,24,26,24,20,30,20,30),4,2)),e_u=0.15, e_s=0.9, e_h=0){
  ### Xs_dim_list specifies the matrix dimensions of multiple modality matrix Xs, with each element as
  ### Xs_dim_list[[i]]=[nrow(Xs[[i]]),ncol(Xs[[i]])]=[N,Ms[i]], i=1,...,K

  # the perturbation generation function used in iNMF data generation
  ptb_gen <- function(dims, mod_dims, rate, sign=1){
    n_mod <- nrow(mod_dims)
    matr <- matrix(0,dims[1],dims[2])
    r_idx <- 0
    for (i in 1:n_mod){
      c_idx <- 0
      for (j in 1:n_mod){
        rows <- mod_dims[i,1]
        cols <- mod_dims[j,2]
        if (rbinom(1,1,prob = rate)==1){
          if (rbinom(1,1,prob = 0.5)==1){
            matr[(r_idx+1):(r_idx+rows/2),(c_idx+1):(c_idx+cols/2)] <- matrix(rbinom(1,size=rows*cols/2,prob =sign ),rows/2,cols)*2-1
          }
          else{
            matr[(r_idx+rows/2+1):(r_idx+rows),(c_idx+1):(c_idx+cols)] <- matrix(rbinom(1,size=(rows-rows/2)*cols,prob =sign ),(rows-rows/2),cols)*2-1
          }
        }
        c_idx <- c_idx + cols
      }
      r_idx <- r_idx + rows
    }
    return(matr)
  }

  # mod_dim_list specifies the module dimensions
  # 1(a)
  K <- length(Xs_dim_list)
  N <- Xs_dim_list[[1]][1]
  Ms <- c()
  for (k in 1:K){
    Ms <- c(Ms,Xs_dim_list[[k]][2])
  }
  D <- nrow(mod_dim_list[[1]])
  e_u <- max(e_u, 0.01)
  W <- matrix(0,N,D)
  Hs_list <- vector('list',K)
  for (k in 1:K){
    Hs_list[[k]] <- matrix(0,D,Ms[k])
  }
  for (k in 1:K){
    r_idx <- 0
    c_idx <- 0
    for (d in 1:D){
      rows <- mod_dim_list[[k]][d,1]
      cols <- mod_dim_list[[k]][d,2]
      W[(r_idx+1):(r_idx+rows),d] <- 1
      Hs_list[[k]][d,(c_idx+1):(c_idx+cols)] <- 1
      r_idx <- r_idx+rows
      c_idx <- c_idx+cols
    }
  }

  new_mod_dim <- vector("list",K)
  for (k in 1:K){
    new_matr <- matrix(1,nrow(mod_dim_list[[k]]),2)
    for (j in 1:nrow(mod_dim_list[[k]])){
      new_matr[j,1] <- mod_dim_list[[k]][j,1]
    }
    new_mod_dim[[k]] <- new_matr
  }
  # 1(b)
  Vs <- vector("list",K)
  Xs <- vector("list",K)
  # 2(b)
  W <- W*matrix(rbeta(N*D,2,2),N,D)*2
  for (k in 1:K){
    # 2(a) 2(b)
    Vs[[k]] <- ptb_gen(c(N,D),new_mod_dim[[k]],e_h)
    # 2(c)
    Vs[[k]] <- Vs[[k]]*matrix(rbeta(N*D,2,2),N,D)*2
    # 1(b)
    Hs_list[[k]] <- Hs_list[[k]]*matrix(rbeta(D*Ms[k],2,2),D,Ms[k])*2
    # 1(c) 2(d)
    Xs[[k]] <- (W+Vs[[k]])%*%Hs_list[[k]]
  }
  # 3(a)
  for (k in 1:K){
    for (n in 1:N){
      for (m in 1:Ms[k]){
        if (rbinom(1,1,prob = e_s)==1){
          if (Xs[[k]][n,m]>0){
            Xs[[k]][n,m] <- 0
          }
          else{
            Xs[[k]][n,m] <- (rbeta(1,2,2)*2)^2
          }
        }
      }
    }
  }
  # 3(b)
  for (k in 1:K){
    Xs[[k]] <- abs(Xs[[k]])+matrix(runif(N*Ms[k],-e_u,e_u),N,Ms[k])
    Xs[[k]] <- t(Xs[[k]])
  }

  # return(Xs)
  temp_data <- Xs
  temp_id <- rep(c(1:nrow(mod_dim_list[[1]])),mod_dim_list[[1]][,2])
  res <- list(temp_data,temp_id)
  return(res)
}
