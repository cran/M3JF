#' Generate the simulated dataset with three modalities with the package crimmix
#' @param nclust number of clusters
#' @param n_byClust number of samples per cluster
#' @param feature_nums number of features in each modality
#' @param noises percentage of noise adding to each modality
#' @param props proportion of cluster related features in each modality
#'
#' @return res, a list of length 2, where the first element is a list of simulated data, while the second element is a vector indicating the true label of each sample.
#' @export
#'
#' @examples
#' crimmix_data <- crimmix_data_gen(nclust=4, n_byClust=c(10,20,5,25), feature_nums=c(1000,500,5000), noises=c(0.5,0.01,0.3),props=c(0.005,0.01,0.02))


crimmix_data_gen <- function(nclust=4, n_byClust=c(10,20,5,25), feature_nums=c(1000,500,5000), noises=c(0.5,0.01,0.3),props=c(0.005,0.01,0.02)){
  # for the gaussion distribution
  means <- rep(2, nclust)
  sds <- rep(1, nclust)
  params_norms <- mapply(function(m, sd) {
    c(mean=m, sd=sd)}
    ,means,sds,SIMPLIFY=FALSE)
  sigma <- noises[1]
  dat1 <- simulateY(nclust = nclust,
                    n_byClust = n_byClust,J=feature_nums[1],
                    flavor="normal",
                    prop = props[1],
                    params = params_norms,
                    noise = sigma)

  # for the binary distribution
  params_bin <- list(c(p=0.6))
  p_noise <- noises[2]
  dat2 <- simulateY(nclust = nclust,
                    n_byClust = n_byClust,J=feature_nums[2],
                    flavor = "binary",
                    params = params_bin,
                    prop = props[2], noise = p_noise)

  #for the beta-like distribution
  params_beta <- list(c(mean1=-2,mean2=2,sd1=0.5,sd2=0.5))
  sd <- noises[3]
  dat3 <- simulateY(nclust = nclust,
                    n_byClust = n_byClust,J=feature_nums[3],
                    flavor = "beta",params = params_beta,
                    prop = props[3],noise = sd)

  data1 <- dat1$data
  data2 <- dat2$data
  data3 <- dat3$data
  temp_data <- list(data1, data2, data3)
  temp_id <- rep(c(1:nclust),n_byClust)
  res <- list(temp_data,temp_id)
  return(res)
}
