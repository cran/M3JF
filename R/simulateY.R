#' The simulation data generation function from R package CrIMMix, importing the package is so annoying, we only include the simulateY function here.
#' @param nclust number of clusters
#' @param n_byClust number of samples per cluster
#' @param feature_nums number of features in each modality
#' @param prop proportion of cluster related features
#' @param noise percentage of noise adding to each modality
#' @param flavor a vector indicating the data type
#' @param params a list indicating the mean and standard derivation of the simulated data
#'
#' @return data_res, the simulation data list.
#' @export
#' @importFrom stats rbinom rnorm
#' @examples
#' nclust=4
#' n_byClust=c(10,20,5,25)
#' feature_nums=c(1000,500,5000)
#' noises=c(0.5,0.01,0.3)
#' props=c(0.005,0.01,0.02)
#' means <- rep(2, nclust)
#' sds <- rep(1, nclust)
#' params_norms <- mapply(function(m, sd) {c(mean=m, sd=sd)},means,sds,SIMPLIFY=FALSE)
#' sigma <- noises[1]
#' dat1 <- simulateY(nclust = nclust,n_byClust = n_byClust,J=feature_nums[1],flavor="normal",prop = props[1],params = params_norms[[1]],noise = sigma)
simulateY <- function(nclust = 4, n_byClust = c(10,20,5,25), J=1000, prop = 0.01, noise = 0.1,
                       flavor = c("normal", "beta", "binary"), params = list(c(mean = 1,
                                                                               sd = 1)))
{
  if (!is.numeric(nclust)) {
    stop("nclust must be an integer")
  }
  if (length(n_byClust) == 1) {
    n_byClust <- rep(n_byClust, nclust)
  }
  if (length(n_byClust) != nclust) {
    stop(sprintf("length of n_byClust must be equal to nclust=%s",nclust))
  }
  n = sum(n_byClust)
  if (length(prop) == 1) {
    prop <- rep(prop, nclust)
  }
  if (length(prop) != nclust) {
    stop(sprintf("length of prop must be equal to nclust=%s",nclust))
  }
  if (length(params) == 1) {
    params <- replicate(nclust, params, simplify = TRUE)
  }
  if (length(params) != nclust) {
    stop(sprintf("number of row of params must be equal to nclust=%s",nclust))
  }
  flavor <- match.arg(flavor)
  true.clusters = mapply(function(clust, n) {
    c(rep(sprintf("C%s", clust), n))
  }, 1:nclust, n_byClust, SIMPLIFY = FALSE) %>% unlist
  flavor <- match.arg(flavor)
  if (flavor == "beta") {
    ep <- c("mean1", "mean2", "sd1", "sd2")
    p_names <- params[[1]] %>% names()
    mm <- match(ep, p_names)
    if (any(is.na(mm))) {
      str <- sprintf("('%s')", paste(ep, collapse = "','"))
      stop("Argument 'params' should for normal contain columns named ",str)
    }
  }
  if (flavor == "normal") {
    ep <- c("mean", "sd")
    p_names <- params[[1]] %>% names()
    mm <- match(ep, p_names)
    if (any(is.na(mm))) {
      str <- sprintf("('%s')", paste(ep, collapse = "','"))
      stop("Argument 'params' for beta should contain columns named ",str)
    }
  }
  if (flavor == "binary") {
    ep <- c("p")
    p_names <- params[[1]] %>% names()
    mm <- match(ep, p_names)
    if (any(is.na(mm))) {
      str <- sprintf("('%s')", paste(ep, collapse = "','"))
      stop("Argument 'params' for binary should contain columns named ",str)
    }
  }
  data_sim = matrix(rnorm(n * J, mean = 0, sd = noise),nrow = n,ncol = J)
  simulate_data <- switch(flavor, normal = simulate_norm, beta = simulate_beta,binary = simulate_binary)
  data_sim <- mapply(simulate_data, n_byClust, round(prop * J), params, J, SIMPLIFY = FALSE)
  data <- do.call(rbind, lapply(data_sim, function(dd) dd$C))
  if (flavor == "binary") {
    data <- data + matrix(rbinom(n * J, 1, noise), nrow = n,ncol = J)
    data <- pmin(data, 1)
  }
  else if (flavor == "beta") {
    rev.logit <- function(x) 1/(1 + exp(-x + matrix(rnorm(n*J, mean = 0, sd = noise), nrow = n, ncol = J)))
    data <- rev.logit(data)
  }
  else {
    data <- data + matrix(rnorm(n * J, mean = 0, sd = noise),nrow = n, ncol = J)
  }
  colnames(data) <- sprintf("gene%s", 1:ncol(data))
  positive <- lapply(data_sim, function(dd) colnames(data)[dd$positive])
  data_res <- list(data = data, positive = positive, true.clusters = true.clusters)
  return(data_res)
}


################# normal distribution ####################

simulate_norm <- function(n,## patient
                          j, ## biomark
                          params, ## parameter distribution
                          J ## all biomarkers
){
  m <- params["mean"]
  sd <- params["sd"]
  c= matrix(rnorm(n*j,mean=m,sd=sd),ncol=j,nrow=n)
  c= cbind(c, matrix(0,nrow=n, ncol=J-j))
  idx <- sample(1:ncol(c))
  positive <- sapply(1:j, function (ll) which(ll== idx))
  return(list(C= c[,idx], positive=positive))
}
################# distribution for methylation ####################

simulate_beta <- function(n,## patient
                          j, ## biomark
                          params, ## parameters
                          J ## all biomarkers
){
  m1 <- params["mean1"]
  sd1 <- params["sd1"]
  m2 <- params["mean2"]
  sd2 <- params["sd2"]
  c_1 <- matrix(rnorm(n*3*j/4,mean=m1,sd=sd1),ncol=3*j/4,nrow=n) ## Hypo (neg mean)
  c_2 <- matrix(rnorm(n*j/4,mean=m2,sd=sd2),ncol=j/4,nrow=n) ## Hyper (pos mean)

  c = cbind(c_1, c_2)
  c = cbind(c, matrix(logit(runif(n*(J-j))), ncol=J-j, nrow=n))
  idx <- sample(1:ncol(c)) ## sample columns
  positive <- sapply(1:j, function(ll) which(ll== idx))
  return(list(C= c[,idx], positive=positive))
}
################# binary distribution ####################

simulate_binary <- function(n,## patient
                            j,## biomark
                            params,
                            J ## all biomarkers
){
  p <- params["p"]
  c=matrix(rbinom(n*j,1,p),ncol=j,nrow=n)
  c= cbind(c, matrix(0,nrow=n, ncol=J-j))
  idx <- sample(1:ncol(c))
  positive <- sapply(1:j, function (ll) which(ll== idx))
  return(list(C= c[,idx], positive=positive))
}

logit <- function (x) {
  log(x) - log(1 - x)
}
