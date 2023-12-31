% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iNMF_data_gen.R
\name{iNMF_data_gen}
\alias{iNMF_data_gen}
\title{Generate the simulated dataset with three modalities as the work iNMF}
\usage{
iNMF_data_gen(Xs_dim_list=list(c(100,100),c(100,100),c(100,100)),
mod_dim_list=list(matrix(c(20,30,20,30,20,30,20,30),4,2),
matrix(c(20,20,30,30,20,30,20,30),4,2),
matrix(c(26,24,26,24,20,30,20,30),4,2)),e_u=0.15, e_s=0.9, e_h=0)
}
\arguments{
\item{Xs_dim_list}{a list of data matrix dimensions for multiple modality data}

\item{mod_dim_list}{a list of the dimensions of each cluster and their features}

\item{e_u}{the level of uniform noise}

\item{e_s}{signal to noise ratio}

\item{e_h}{block adding probability}
}
\value{
res, a list of length 2, where the first element is a list of simulated data, while the second element is a vector indicating the true label of each sample.
}
\description{
Generate the simulated dataset with three modalities as the work iNMF
}
\examples{
iNMF_data <- iNMF_data_gen(Xs_dim_list=list(c(100,100),c(100,100),c(100,100)),
mod_dim_list=list(matrix(c(20,30,20,30,20,30,20,30),4,2),
matrix(c(20,20,30,30,20,30,20,30),4,2),
matrix(c(26,24,26,24,20,30,20,30),4,2)),e_u=0.15, e_s=0.9, e_h=0)
}
