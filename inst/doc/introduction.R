## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  #Install
#  install.packages('M3JF')
#  #Load
#  library(M3JF)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages('/path/to/file/M3JF.tar.gz',repos=NULL,type="source")

## ----eval=FALSE---------------------------------------------------------------
#  library(InterSIM)
#  
#  sim.data <- InterSIM(n.sample=500, cluster.sample.prop = c(0.20,0.30,0.27,0.23),
#                       delta.methyl=5, delta.expr=5, delta.protein=5,p.DMP=0.2, p.DEG=NULL,
#                       p.DEP=NULL,sigma.methyl=NULL, sigma.expr=NULL,
#                       sigma.protein=NULL,cor.methyl.expr=NULL,
#                       cor.expr.protein=NULL,do.plot=FALSE, sample.cluster=TRUE,
#                       feature.cluster=TRUE)
#  sim.methyl <- sim.data$dat.methyl
#  sim.expr <- sim.data$dat.expr
#  sim.protein <- sim.data$dat.protein
#  data_list <- list(sim.methyl, sim.expr, sim.protein)

## ----eval=FALSE---------------------------------------------------------------
#  truelabel = sim.data$clustering.assignment$cluster.id

## ----eval=FALSE---------------------------------------------------------------
#  #Build similarity matrices for your data with SNFtool
#  library(SNFtool)
#  library(dplyr)
#  WL_dist1 <- lapply(data_list,function(x){
#    dd <- x%>%as.matrix
#    w <- dd %>% dist2(dd) %>% affinityMatrix(K = 10, sigma = 0.5)
#  })
#  #Assign the interval of k according to your data
#  k_list = 2:10
#  #Initialize the varible
#  clu_eval <- RotationCostBestGivenGraph(W,k_list)
#  #The most proper is the one with minimal rotation cost
#  best_k = k_list[which.min(clu_eval)]

## ----eval=FALSE---------------------------------------------------------------
#  #Assign the parameters
#  lambda = 0.01
#  theta = 10^-6
#  k = best_k
#  res = M3JF(data_list,lambda,theta,k)

## ----eval=FALSE---------------------------------------------------------------
#  library(SNFtool)
#  #Calculate the NMI of *M3JF*
#  M3JF_res = M3JF(data_list,lambda,theta,k)
#  M3JF_cluster = M3JF_res$clusters
#  M3JF_NMI = cal_NMI(true_label,M3JF_cluster)
#  #Calculate the ARI of *M3JF*
#  library(mclust)
#  M3JF_ARI = adjustedRandIndex(true_label,M3JF_cluster)

