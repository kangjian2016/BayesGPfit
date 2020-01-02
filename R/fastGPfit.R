#' Create 256 colors graduately transitioning from
#' Blue to Yellow to Red.
#'
#' @param num  A integer number to specify the number of colors to generate. The default value is 256.
#' @return A vector of RGB colors
#' @author Jian Kang <jiankang@umich.edu>
#' @examples
#' colors = GP.create.cols(101L)
#' require(graphics)
#' filled.contour(volcano,col=colors,nlevels=length(colors)-1,asp=1)
#' filled.contour(volcano,color.palette = GP.create.cols, nlevels=256, asp = 1)
#' @export
GP.create.cols = function(num=256L){
  num = as.integer(num)
  if(num<1L || num>256L){
    stop("The number of colors has to be between 1 and 256!")
  }
  R = c(rep(0,length=96),
        seq(0.015625,0.984375,length=63),
        rep(1,length=65),
        seq(0.984375,0.5,length=32))
  G = c(rep(0,length=32),
        seq(0.015625,0.984375,length=63),
        rep(1,length=65),
        seq(0.984375,0.015625,length=63),
        rep(0,length=33))
  B = c(seq(0.5,0.984375,length=32),
        rep(1,length=64),
        seq(0.984375,0.015625,length=63),
        rep(0,length=97))
  return(rgb(cbind(R[1:num],G[1:num],B[1:num])))
}


#' Create spatial grids.
#'
#' @param d  An integer number for the dimension of the space. The default value is 1.
#' @param num_grids An integer number for the number of grids in each dimension. The default value is 50.
#' @param grids_lim A vector of two real numbers for the range of the grids in each dimension. The default value is c(-1,1).
#' @param random A logical value indicating whether each dimension of the grids is generated from a uniform distribution or fixed as equally-spaced.
#' @return A matrix with d columns and num_grids^d rows.
#' @author Jian Kang <jiankang@umich.edu>
#'
#' @examples
#' x = GP.generate.grids(d=2L)
#' require(lattice)
#' y = sin(abs(x[,1]+x[,2]))
#' levelplot(y~x[,1]+x[,2])
#' @export
GP.generate.grids = function(d = 1L, num_grids = 50L, grids_lim = c(-1,1),random=FALSE){
  if(random){
    base_grids = runif(num_grids,grids_lim[1],grids_lim[2])
  }
  else{
    base_grids = seq(grids_lim[1],grids_lim[2],length=num_grids)
  }

  grids_list = list()
  for(i in 1:d){
    grids_list[[i]] = base_grids
  }
  grids = expand.grid(grids_list)
  names(grids) = paste("x",1:d,sep="")
  return(as.matrix(grids))
}

#' Compute eigen values for the standard modified exponential squared correlation kernel.
#'
#' @param poly_degree A positive integer number specifies the highest degree of Hermite polynomials. The default value is 10L.
#' @param a A positive real number specifying the concentration parameter in the modified exponetial squared kernel. The larger value the more the GP concentrates around the center. The default value is 0.01.
#' @param b A positive real number specifying the smoothness parameter in the modeified exponetial squared kernel. The smaller value the smoother the GP is. The default value is 1.0.
#' @param d A positive integer number specifying the dimension of grid points.
#' @return A matrix represents a set of eigen functions evaluated at grid points.
#' The number of rows is equal to the number of grid points. The number of columns is choose(poly_degree+d,d), where d is the dimnension of the grid points.
#' @details
#' Compute eigen values of the standard modified exponential squared kernel on d-dimensional grids
#'
#' \eqn{cor(X(s_1),X(s_2)) = \exp{-a*(s_1^2+*s_2^2)-b*(s_1-s_2)^2}}
#'
#' where \eqn{a} is the concentration parameter and \eqn{b} is the smoothness parameter. The expected ranges of each coordinate is from -6 to 6.
#'
#'@author Jian Kang <jiankang@umich.edu>
#'
#' @examples
#' library(BayesGPfit)
#' Lambda = GP.eigen.value(poly_degree=10L,a=0.01,b=0.5,d=2)
#' plot(Lambda)
#'
#' @export
GP.eigen.value = function(poly_degree=10,a=1,b=1,d=2){
  cn = sqrt(a^2+2*a*b)
  A = a+b+cn
  B = b/A
  idx = c(0, choose(0:poly_degree+d,d))
  idxlist = sapply(1:(poly_degree+1),function(i) return((idx[i]+1):idx[i+1]))
  k = GP.num.eigen.funs(poly_degree=poly_degree,d=d)
  value = rep(NA,length=k)
  dvalue = (sqrt(pi/A))^d*B^(1:(poly_degree+1))
  for(i in 1:(poly_degree+1))
    value[idxlist[[i]]] = dvalue[i]
  return(value)
}

GP.num.eigen.funs = function(poly_degree=10,d=2){
  return(choose(poly_degree+d,d))
}


sapply.pb <- function(X, FUN, ...) {
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    utils::setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}



GP.fast.gibbs.sampler = function(y_l,lambda,
                            num_results = 500,
                            iters_between_results = 2,
                            burn_in = 500,
                            a_sigma = 0.01,
                            b_sigma = 0.01,
                            a_tau = 0.01,
                            b_tau = 0.01,
                            progress_bar = FALSE){

  L = length(y_l)
  sigma2 = 1/rgamma(1,a_sigma,b_sigma)
  tau2 = 1/rgamma(1,a_tau,b_tau)
  mcmc_sample = matrix(NA,nrow=num_results,ncol=L+2L)
  if(progress_bar){
    pb = txtProgressBar(style=3)
  }
  total_iter = burn_in + num_results*iters_between_results
  k = 1
  for(iter in 1:total_iter){
    mu_theta = y_l/(sigma2/(lambda*tau2)+1)
    sd_theta = sqrt(sigma2*lambda*tau2/(sigma2+lambda*tau2))
    theta_l = rnorm(L, mean = mu_theta, sd=sd_theta)
    sigma2 = 1/rgamma(1,a_sigma+0.5*L,b_sigma+0.5*sum((y_l-theta_l)^2))
    tau2 = 1/rgamma(1,a_tau+0.5*L,b_tau+0.5*sum(theta_l^2/lambda))
    if((iter > burn_in) & ((iter-burn_in)%%iters_between_results==0)& k<=num_results){
      mcmc_sample[k,] = c(theta_l,sigma2,tau2)
      k = k + 1
    }
    if(progress_bar){
      setTxtProgressBar(pb,iter/total_iter)
    }
  }
  if(progress_bar){
    close(pb)
  }
  return(mcmc_sample)
}

#'Compute the standardized grids
#'@param grids A matrix where rows represent grid points and columns are coordinates.
#'@param center A vector of real numbers specifying the centroid parameters in the modified exponential squared kernel. The default value is NULL and set to the center of the grid points: apply(x,2,mean).
#'@param scale A vector of positive numbers specifying the scale parameters in the modified exponential squared kernel. The default value is NULL and set to values such that grid points in a range of (-max_range,max_range) in each dimension.
#'@param max_range A positive real number indicating the maximum range of the grid points to specify the scale parameter. The default value is NULL and set to 6.
#'@return A matrix where rows represent the standardized grids.
#'@author Jian Kang <jiankang@umich.edu>
#'@examples
#'library(BayesGPfit)
#'grids = GP.generate.grids(d=2L)
#'std_grids = GP.std.grids(grids)
#'plot(grids[,1],std_grids[,1],asp=1,type="l")
#'abline(a=0,b=1,lty=2)
#'@export
GP.std.grids = function(grids,center=apply(grids,2,mean),scale=NULL,max_range = 6){
  c_grids = t(grids) - center
  if(is.null(scale)){
    max_grids =pmax(apply(c_grids,1,max),-apply(c_grids,1,min))
    scale=as.numeric(max_grids/max_range)
  }
  grids = t(c_grids/scale)
  return(grids)
}

#' Compute eigen functions for the standard modified exponential squared correlation kernel.
#' @title Compute eigen functions
#' @param grids A matrix where rows represent points and columns are coordinates.
#' @param poly_degree A integer number specifies the highest degree of Hermite polynomials. The default value is 10L.
#' @param a A positive real number specifies the concentration parameter in the modified exponetial squared kernel. The larger value the more the GP concentrates around the center. The default value is 0.01.
#' @param b A positive real number specifies the smoothness parameter in the modeified exponetial squared kernel. The smaller value the smoother the GP is. The default value is 1.0.
#' @return A matrix represents a set of eigen functions evaluated at grid points.
#' The number of rows is equal to the number of grid points. The number of columns is choose(poly_degree+d,d), where d is the dimnension of the grid points.
#' @details
#' Compute eigen values of the standard modified exponential squared kernel on d-dimensional grids
#'
#' \eqn{cor(X(s_1),X(s_2)) = \exp{-a*(s_1^2+*s_2^2)-b*(s_1-s_2)^2}}
#'
#' where \eqn{a} is the concentration parameter and \eqn{b} is the smoothness parameter. The expected ranges of each coordinate is from -6 to 6.
#'
#'
#' @author Jian Kang <jiankang@umich.edu>
#' @examples
#' library(lattice)
#' grids = GP.generate.grids(d=2L)
#' Psi_mat = GP.eigen.funcs.fast(grids)
#' fig = list()
#' for(i in 1:4){
#'    fig[[i]] = levelplot(Psi_mat[,i]~grids[,1]+grids[,2])
#' }
#' plot(fig[[1]],split=c(1,1,2,2),more=TRUE)
#' plot(fig[[2]],split=c(1,2,2,2),more=TRUE)
#' plot(fig[[3]],split=c(2,1,2,2),more=TRUE)
#' plot(fig[[4]],split=c(2,2,2,2))
#' @useDynLib BayesGPfit,.registration = TRUE, Wrapper_R_GP_eigen_funcs
#' @export
GP.eigen.funcs.fast = function(grids,poly_degree=10L,a=0.01,b=1.0){
  num_funcs = GP.num.eigen.funs(poly_degree,ncol(grids))
  eigen_funcs = vector(mode="numeric",length=num_funcs*nrow(grids))
  res = .C('Wrapper_R_GP_eigen_funcs',
           eigen_funcs_R = as.double(eigen_funcs),
           num_funcs_R = as.integer(num_funcs),
           grids_R = as.double(grids),
           grids_size_R = as.integer(nrow(grids)),
           dim_R = as.integer(ncol(grids)),
           poly_degree_R = as.integer(poly_degree),
           a_R = as.double(a),
           b_R = as.double(b))
  return(matrix(res$eigen_funcs_R,nrow=res$grids_size_R,ncol=res$num_funcs_R))
}

#'Fast Bayesian fitting of Gaussian process regression on regular grid points with the modified exponential sqaured kernel.
#'@title Fast Bayesian fitting of Gaussian process
#'@param y A vector of real numbers as the observations for the reponse variable.
#'@param x A matrix of real numbers as grid points where rows are observations and columns are coordinates.
#'@param poly_degree A integer number to specify the highest degree of Hermite polynomials. The default value is 10L.
#'@param a A positive real number to specify the concentration parameter in the standard modified exponential squared kernel. The larger value the more the GP concentrates around the center. The default value is 0.01.
#'@param b A positive real number to specify the smoothness parameter in the standard modified exponential squared kernel. The smaller value the smoother the GP is. The default value is 1.0.
#'@param center A vector of real numbers specifying the centroid parameters in the modified exponential squared kernel. The default value is NULL and set to the center of the grid points: apply(x,2,mean).
#'@param scale A vector of positive numbers specifying the scale parameters in the modified exponential squared kernel. The default value is NULL and set to values such that grid points in a range of (-max_range,max_range) in each dimension.
#'@param max_range A positive real number indicating the maximum range of the grid points to specify the scale parameter. The default value is NULL and set to 6.
#'@param num_results An integer number to specify the number of posterior samples to save over MCMC iterations.
#'@param iters_between_results An integer number to specify the number of iterations to skip between two saved iterations.
#'@param burn_in An integer number to specify the burn-in number. The default value is 500L.
#'@param a_sigma A real number for the shape parameter in the Gamma prior of sigma2. The default value is 0.01.
#'@param b_sigma A real number for the rate parameter in the Gamma prior of sigma2. The default value is 0.01.
#'@param a_tau  A real number for the shape parameter in the Gamma prior of tau2. The default value is 0.01.
#'@param b_tau A real number for the rate parameter in the Gamma prior of tau2. The default value is 0.01.
#'@param progress_bar A logical value to indicate whether a progress bar will be shown.
#'
#'@return A list of variables including the model fitting results
#'\describe{
#'  \item{f}{A vector of real numbers for the posterior mean of the fitted curve.}
#'  \item{x}{A matrix of real numbers for the grid points where rows are observations and columns are coordinates.}
#'  \item{work_x}{A matrix of real numbers for the standardized grid points for the model fitting. It has the same dimension as "x".}
#'  \item{sigma2}{A real number for the posterior mean of the variance parameter of random errors.}
#'  \item{tau2}{A real number for the posterior mean of the variance parameter for the Gaussian process prior.}
#'  \item{theta}{A vector of real numbers for the posterior mean of the basis coefficients for the Gaussian process.}
#'  \item{Xmat}{A matrix real numbers for the basis functions evaluated at the standardized grid points (work_x), where rows are observations and columns are the basis functions}
#'  \item{grid_size}{A real scalar for the grid size}
#'  \item{center}{A vector of real numbers for the centroid parameters in the modified exponential squared kernel.}
#'  \item{scale}{A vector of positive numbers for the scale parameters in the modified exponential squared kernel.}
#'  \item{max_range}{A positive real number indicating the maximum range of the grid points to specify the scale parameter.}
#'  \item{poly_degree}{An integer number to specify the highest degree of Hermite polynomials.}
#'  \item{a}{A positive real number to specify the concentration parameter in the standard modified exponential squared kernel.}
#'  \item{b}{A positive real number to specify the smoothness parameter in the standard modified exponential squared kernel.}
#'  \item{mcmc_sample}{A matrix of real numbers for saved MCMC samples.}
#'  \item{elapsed}{A real number indicating the computing time in second}
#'}
#'
#'@author Jian Kang <jiankang@umich.edu>
#'
#'@examples
#'library(BayesGPfit)
#'library(lattice)
#'set.seed(1224)
#'dat = list()
#'dat$x = GP.generate.grids(d=2,num_grids = 100)
#'curve = GP.simulate.curve.fast(dat$x,a=0.01,b=0.5,poly_degree=20L)
#'dat$f = curve$f + rnorm(length(curve$f),sd=1)
#'fit = GP.fast.Bayes.fit(dat$f,dat$x,a=0.01,b=0.5,poly_degree=20L,progress_bar = TRUE)
#'plot(GP.plot.curve(dat,main="Data"),split=c(1,1,2,2),more=TRUE)
#'plot(GP.plot.curve(curve,main="True curve"),split=c(1,2,2,2),more=TRUE)
#'plot(GP.plot.curve(fit,main="Posterior mean estimates"),split=c(2,2,2,2),more=TRUE)
#'
#'
#'@export
GP.fast.Bayes.fit = function(y,x,poly_degree = 10L,
                             a = 0.01,b = 1.0,center=NULL,
                             scale=NULL,max_range=NULL,
                             num_results = 500L,
                             iters_between_results = 2L,
                             burn_in = 500L,
                             a_sigma = 0.01,
                             b_sigma = 0.01,
                             a_tau = 0.01,
                             b_tau = 0.01,
                             progress_bar = FALSE){
  elapsed = proc.time()[3]
  x = cbind(x)
  d = ncol(x)
  if(is.null(center)){
    center = apply(x,2,mean)
  }
  c_x = t(x) - center
  if(is.null(max_range) & is.null(scale)){
    max_range = 6
  }
  if(is.null(scale)){
    max_grids =pmax(apply(c_x,1,max),-apply(c_x,1,min))
    scale=as.numeric(max_grids/max_range)
  }
  work_x = GP.std.grids(x,center=center,scale=scale,max_range=max_range)
  Xmat = GP.eigen.funcs.fast(work_x,poly_degree=poly_degree,a=a,b=b)
  lambda = GP.eigen.value(poly_degree,a,b,d=ncol(x))

  #grid_size = 1/mean(sapply(1:ncol(Xmat),function(i) sum(Xmat[,i]*Xmat[,i])))
  grid_size = 1/sum(Xmat[,1]*Xmat[,1])
  y_l = t(Xmat)%*%y*grid_size
  mcmc_sample=GP.fast.gibbs.sampler(y_l,lambda,
                                    num_results = num_results,
                                    iters_between_results = iters_between_results,
                                    burn_in = burn_in,
                                    a_sigma = a_sigma,
                                    b_sigma = b_sigma,
                                    a_tau = a_tau,
                                    b_tau = b_tau,
                                    progress_bar = progress_bar)

  theta_l = apply(mcmc_sample[,1:length(y_l)],2,mean)
  sigma2_star = apply(mcmc_sample[,length(y_l)+1,drop=FALSE],2,mean)
  tau2 = apply(mcmc_sample[,length(y_l)+2,drop=FALSE],2,mean)
  f = Xmat%*%theta_l
  return(list(f = f,
              x = x,
              work_x = work_x,
              sigma2 =  sigma2_star/grid_size,
              tau2 = tau2,
              theta = theta_l,
              Xmat = Xmat,
              grid_size = grid_size,
              center = center,
              scale = scale,
              max_range = max_range,
              poly_degree = poly_degree,
              a = a,
              b = b,
              mcmc_sample = mcmc_sample,
              elapsed=proc.time()[3]-elapsed))
}


#'Regular Bayesian fitting of Gaussian process regression on regular grid points with the modified exponential sqaured kernel.
#'
#'@param y A vector of real numbers as the observations for the reponse variable.
#'@param x A matrix of real numbers as grid points where rows are observations and columns are coordinates.
#'@param poly_degree A integer number to specify the highest degree of Hermite polynomials. The default value is 10L.
#'@param a A positive real number to specify the concentration parameter in the standard modified exponential squared kernel. The larger value the more the GP concentrates around the center. The default value is 0.01.
#'@param b A positive real number to specify the smoothness parameter in the standard modified exponential squared kernel. The smaller value the smoother the GP is. The default value is 1.0.
#'@param center A vector of real numbers specifying the centroid parameters in the modified exponential squared kernel. The default value is NULL and set to the center of the grid points: apply(x,2,mean).
#'@param scale A vector of positive numbers specifying the scale parameters in the modified exponential squared kernel. The default value is NULL and set to values such that grid points in a range of (-max_range,max_range) in each dimension.
#'@param max_range A positive real number indicating the maximum range of the grid points to specify the scale parameter. The default value is NULL and set to 6.
#'@param num_results An integer number to specify the number of posterior samples to save over MCMC iterations.
#'@param iters_between_results An integer number to specify the number of iterations to skip between two saved iterations.
#'@param burn_in An integer number to specify the burn-in number. The default value is 500L.
#'@param a_sigma A real number for the shape parameter in the Gamma prior of sigma2. The default value is 0.01.
#'@param b_sigma A real number for the rate parameter in the Gamma prior of sigma2. The default value is 0.01.
#'@param a_zeta  A real number for the shape parameter in the Gamma prior of zeta. The default value is 0.01.
#'@param b_zeta A real number for the rate parameter in the Gamma prior of zeta. The default value is 0.01.
#'@param progress_bar A logical value to indicate whether a progress bar will be shown.
#'
#'@return A list of variables including the model fitting results
#'\describe{
#'  \item{f}{A vector of real numbers for the posterior mean of the fitted curve.}
#'  \item{x}{A matrix of real numbers for the grid points where rows are observations and columns are coordinates.}
#'  \item{work_x}{A matrix of real numbers for the standardized grid points for the model fitting. It has the same dimension as "x".}
#'  \item{sigma2}{A real number for the posterior mean of the variance parameter of random errors.}
#'  \item{tau2}{A real number for the posterior mean of the variance parameter for the Gaussian process prior.}
#'  \item{theta}{A vector of real numbers for the posterior mean of the basis coefficients for the Gaussian process.}
#'  \item{Xmat}{A matrix real numbers for the basis functions evaluated at the standardized grid points (work_x), where rows are observations and columns are the basis functions}
#'  \item{grid_size}{A real scalar for the grid size}
#'  \item{center}{A vector of real numbers for the centroid parameters in the modified exponential squared kernel.}
#'  \item{scale}{A vector of positive numbers for the scale parameters in the modified exponential squared kernel.}
#'  \item{max_range}{A positive real number indicating the maximum range of the grid points to specify the scale parameter.}
#'  \item{poly_degree}{An integer number to specify the highest degree of Hermite polynomials.}
#'  \item{a}{A positive real number to specify the concentration parameter in the standard modified exponential squared kernel.}
#'  \item{b}{A positive real number to specify the smoothness parameter in the standard modified exponential squared kernel.}
#'  \item{mcmc_sample}{A matrix of real numbers for saved MCMC samples.}
#'  \item{elapsed}{A real number indicating the computing time in second.}
#'}
#'
#'@author Jian Kang <jiankang@umich.edu>
#'
#'@examples
#'library(BayesGPfit)
#'library(lattice)
#'set.seed(1227)
#'dat = list()
#'dat$x = GP.generate.grids(d=2,num_grids = 100)
#'curve = GP.simulate.curve.fast(dat$x,a=0.01,b=0.5,poly_degree=20L)
#'dat$f = curve$f + rnorm(length(curve$f),sd=1)
#'fast_fit = GP.fast.Bayes.fit(dat$f,dat$x,a=0.01,b=0.5,poly_degree=20L,progress_bar = TRUE)
#'reg_fit = GP.Bayes.fit(dat$f,dat$x,a=0.01,b=0.5,poly_degree=20L,progress_bar = TRUE)
#'mse = c(reg = mean((reg_fit$f - curve$f)^2),
#'        fast = mean((fast_fit$f - curve$f)^2))
#'print(mse)
#'plot(GP.plot.curve(curve,main="True curve"),split=c(1,2,2,2),more=TRUE)
#'plot(GP.plot.curve(fast_fit,main="Posterior mean estimates (fast)"),split=c(2,2,2,2),more=TRUE)
#'plot(GP.plot.curve(reg_fit,main="Posterior mean estimates (Regular)"),split=c(2,1,2,2))
#'@export
GP.Bayes.fit = function(y,x,poly_degree=60,a = 0.01,b = 20,
                        num_results = 500L,
                        iters_between_results = 2L,
                        burn_in = 500L,
                            a_sigma = 0.01,
                            b_sigma = 0.01,
                            a_zeta = 0.01,
                            b_zeta = 0.01,center=NULL,
                        scale=NULL,max_range=NULL,
                        progress_bar = FALSE){

  elapsed = proc.time()[3]
  x  = cbind(x)
  if(is.null(center)){
    center = apply(x,2,mean)
  }
  c_x = t(x) - center
  if(is.null(max_range) & is.null(scale)){
    max_range = 6
  }
  if(is.null(scale)){
    max_grids =pmax(apply(c_x,1,max),-apply(c_x,1,min))
    scale=as.numeric(max_grids/max_range)
  }

  work_x = GP.std.grids(x,center=center,scale=scale)
  Xmat = GP.eigen.funcs.fast(work_x,poly_degree,a,b)
  lambda = GP.eigen.value(poly_degree,a,b,d=ncol(work_x))

  sres = svd(Xmat)
  U = sres$u
  D = sres$d
  V = sres$v
  y_s = t(U)%*%y
  L = length(y_s)
  inv_sigma2 = 1
  zeta = 1
  total_iter = burn_in + num_results*iters_between_results
  mcmc_sample = matrix(NA,nrow=num_results,ncol=L+2L)
  if(progress_bar){
    pb = txtProgressBar(style=3)
  }
  k = 1
  for(iter in 1:total_iter){
    #update theta
    inv_sqrt_D2_lambda = 1.0/sqrt(D*D+zeta/lambda)
    theta_s = D*y_s*inv_sqrt_D2_lambda+rnorm(L)/sqrt(inv_sigma2)
    theta_s = V%*%(theta_s*inv_sqrt_D2_lambda)
    #update sigma2
    theta_s_2_lambda = sum(theta_s*theta_s/lambda)
    inv_sigma2 = rgamma(1,a_sigma+L,b_sigma+0.5*sum((y_s - D*t(V)%*%theta_s)^2)+0.5*zeta*theta_s_2_lambda)
    #update zeta
    zeta = rgamma(1,a_zeta+0.5*L, b_zeta+0.5*theta_s_2_lambda*inv_sigma2)
    if((iter > burn_in) & ((iter-burn_in)%%iters_between_results==0)& k<=num_results){
      mcmc_sample[k,] = c(theta_s,1/inv_sigma2,1/(zeta*inv_sigma2))
      k = k + 1
    }
    #mcmc_sample[iter,] = c(theta_s,inv_sigma2,zeta)
    if(progress_bar){
      setTxtProgressBar(pb,iter/total_iter)
    }
  }
  if(progress_bar){
    close(pb)
  }

  theta_l = apply(mcmc_sample[,1:L],2,mean)
  sigma2 = apply(mcmc_sample[,L+1,drop=FALSE],2,mean)
  tau2 = apply(mcmc_sample[,L+2,drop=FALSE],2,mean)
  f = Xmat%*%theta_l
  return(list(f = f,
              x = x,
              work_x = work_x,
              sigma2 =  sigma2,
              tau2 = tau2,
              theta = theta_l,
              Xmat = Xmat,
              center = center,
              scale = scale,
              max_range = max_range,
              poly_degree = poly_degree,
              a = a,
              b = b,
              mcmc_sample = mcmc_sample,
              elapsed=proc.time()[3]-elapsed))
}


#'Summary of posterior inference on the Bayesian Gaussian process regression model
#'
#'@param GP_fit  An output object of function \link{GP.Bayes.fit} or \link{GP.fast.Bayes.fit}. Please refer to them for details.
#'
#'@author Jian Kang <jiankang@umich.edu>
#'
#'@return A list object consisting of the following elements:
#'\describe{
#' \item{mean}{A list object for posterior mean of the target function,consisting of two elements (f is a vector for function values; x is
#'   a vector or matrix for points evaluated).}
#'   \item{work_x}{A matrix of real numbers for the standardized grid points for the model fitting.
#'   It has the same dimension as "x".}
#'    \item{uci}{A list object for 95\% upper bound of the creditible interval
#'    (uci) of the taget function,
#'   consisting of two elements (f is a vector for function values;
#'    x is a vector or matrix for points evaluated).}
#'   \item{lci}{A list object for 95\% lower bound of the creditibel interval (lci) of the taget function,
#'   consisting of two elements (f is a vector for function values; x is
#'   a vector or matrix for points evaluated).}
#'   \item{sigma2}{A vector of posteror mean, the 95\% lcl and ucl for
#'   variance of the random error.}
#'  \item{tau2}{A vector of posterior mean, the 95\% lcl and ucl for
#'   variance of the target function (hyperparameters).}
#'}
#'
#'@examples
#'library(BayesGPfit)
#'library(lattice)
#'set.seed(1224)
#'dat = list()
#'dat$x = GP.generate.grids(d=2,num_grids = 30)
#'curve = GP.simulate.curve.fast(dat$x,a=0.01,b=0.5,poly_degree=20L)
#'dat$f = curve$f + rnorm(length(curve$f),sd=1)
#'fast_fit = GP.fast.Bayes.fit(dat$f,dat$x,a=0.01,b=0.5,poly_degree=20L,progress_bar = TRUE)
#'reg_fit = GP.Bayes.fit(dat$f,dat$x,a=0.01,b=0.5,poly_degree=20L,progress_bar = TRUE)
#'sum_fast_fit = GP.summary(fast_fit)
#'sum_reg_fit = GP.summary(reg_fit)
#'curves = list(mean_fast = sum_fast_fit$mean,
#'              mean = sum_reg_fit$mean,
#'              lci_fast = sum_fast_fit$lci,
#'              lci = sum_reg_fit$lci,
#'              uci_fast = sum_fast_fit$uci,
#'              uci = sum_reg_fit$uci)
#'GP.plot.curves(curves,layout=c(2,3))
#'@export
GP.summary = function(GP_fit){
  L = length(GP_fit$theta)
  theta_mcmc = GP_fit$mcmc_sample[,1:L,drop=FALSE]
  M = nrow(theta_mcmc)
  f_mcmc = sapply.pb(1:M, function(i) GP_fit$Xmat%*%theta_mcmc[i,])
  f_CIs= apply(f_mcmc,1,quantile,prob=c(0.025,0.975))
  res = cbind(GP_fit$f,t(f_CIs))
  colnames(res) = c("mean","lci","uci")

  tau2 = mean(GP_fit$mcmc_sample[,L+1L])
  tau2_lci = quantile(GP_fit$mcmc_sample[,L+1L],prob=0.025)
  tau2_uci = quantile(GP_fit$mcmc_sample[,L+1L],prob=0.975)
  sigma2 = mean(GP_fit$mcmc_sample[,L+2L])
  sigma2_lci = quantile(GP_fit$mcmc_sample[,L+2L],prob=0.025)
  sigma2_uci = quantile(GP_fit$mcmc_sample[,L+2L],prob=0.975)

  return(list(mean=list(f=res[,"mean"],x=GP_fit$x),
              work_x=GP_fit$work_x,
              lci=list(f=res[,"lci"],x=GP_fit$x),
              uci=list(f=res[,"uci"],x=GP_fit$x),
              sigma2=c(mean=sigma2,
                       lci=sigma2_lci,
                       uci=sigma2_uci),
              tau2 = c(mean=tau2,
                       lci = tau2_lci,
                       uci = tau2_uci)))
}

#'Gaussian process predictions
#'@param GP_fit  An output object of function \link{GP.Bayes.fit} or \link{GP.fast.Bayes.fit}. Please refer to them for details.
#'@param newx A matrix of real numbers as new grid points for preditions.
#'@param CI A logical value indicating prediction.
#'@author Jian Kang <jiankang@umich.edu>
#'@return A list object
#'When CI is FALSE, the object consists of three elements:
#'\describe{
#'\item{f}{Posterior predictive mean values of the curves.}
#'\item{x}{The grid points for prediction, i.e. "newx".}
#'\item{work_x}{The standardized grid points for prediction.}
#'}
#'When CI is FALSE, the object consists of four elements:
#'\describe{
#' \item{mean}{A list object for posterior predictive mean of the curve,consisting of two elements (f is a vector for the curve values; x is
#'   a vector or matrix for points evaluated).}
#'   \item{work_x}{A matrix of real numbers for the standardized grid points for the model fitting.
#'   It has the same dimension as "newx".}
#'    \item{uci}{A list object for 95\% upper bound of the predictive creditible interval
#'    (uci) of the curve,
#'   consisting of two elements (f is a vector for curve values;
#'    x is a vector or matrix for points evaluated).}
#'   \item{lci}{A list object for 95\% lower bound of the predictive creditibel interval (lci) of the curve,
#'   consisting of two elements (f is a vector for curve value; x is
#'   a vector or matrix for points evaluated).}
#'}
#'@examples
#'
#'set.seed(1224)
#'traindat = list()
#'traindat$x = GP.generate.grids(d=2,num_grids=30,random=TRUE)
#'testdat = list()
#'testdat$x = GP.generate.grids(d=2,num_grids=30,random=FALSE)
#'curve = GP.simulate.curve.fast(rbind(traindat$x,testdat$x),a=0.01,b=0.5,poly_degree=20L)
#'train_curve = list(f=curve$f[1:nrow(traindat$x)],x=traindat$x)
#'test_curve = list(f=curve$f[nrow(traindat$x)+1:nrow(testdat$x)],x=testdat$x)
#'traindat$f = train_curve$f + rnorm(length(train_curve$f),sd=1)
#'testdat$f = test_curve$f + rnorm(length(test_curve$f),sd=1)
#'fast_fit = GP.fast.Bayes.fit(traindat$f,traindat$x,a=0.01,b=0.5,poly_degree=20L,progress_bar = TRUE)
#'reg_fit = GP.Bayes.fit(traindat$f,traindat$x,a=0.01,b=0.5,poly_degree=20L,progress_bar = TRUE)
#'fast_pred = GP.predict(fast_fit,testdat$x,CI=TRUE)
#'reg_pred = GP.predict(reg_fit,testdat$x,CI=TRUE)
#'pmse = c(fast = mean((fast_pred$mean$f-test_curve$f)^2),
#'         reg = mean((reg_pred$mean$f-test_curve$f)^2))
#'print(pmse)
#'curves = list(true = test_curve,
#'              Bayes = reg_pred$mean,
#'              fast = fast_pred$mean)
#'GP.plot.curves(curves,main="Posterior predictive mean")
#'
#'@export
GP.predict = function(GP_fit,newx,CI=TRUE){
  newx = cbind(newx)
  work_newx = GP.std.grids(newx,center = GP_fit$center,scale=GP_fit$scale,max_range=GP_fit$max_range)
  new_Xmat = GP.eigen.funcs.fast(work_newx,GP_fit$poly_degree,GP_fit$a,GP_fit$b)
  f_pred = new_Xmat%*%GP_fit$theta
  if(CI){
    L = length(GP_fit$theta)
    theta_mcmc = GP_fit$mcmc_sample[,1:L,drop=FALSE]
    M = nrow(theta_mcmc)
    f_mcmc = sapply.pb(1:M, function(i) new_Xmat%*%theta_mcmc[i,])
    f_CIs= apply(f_mcmc,1,quantile,prob=c(0.025,0.975))
    res = cbind(f_pred,t(f_CIs))
    colnames(res) = c("pred","lci","uci")
    return(list(mean=list(f=res[,"pred"],x=newx),work_x=work_newx,
                lci=list(f=res[,"lci"],x=newx),
                uci=list(f=res[,"uci"],x=newx)))
  }
  else{
    return(list(f=f_pred,x=newx,work_x=work_newx))
  }
}


#' Simulate curve on d-dimensional Euclidean space based on Gaussian
#' processes via modified exponential squared kernel.
#'
#' @param x A matrix of real numbers as grid points where rows are observations and columns are coordinates.
#' @param poly_degree A integer number to specify the highest degree of Hermite polynomials. The default value is 10L.
#' @param a A positive real number to specify the concentration parameter in the standard modified exponential squared kernel. The larger value the more the GP concentrates around the center. The default value is 0.01.
#' @param b A positive real number to specify the smoothness parameter in the standard modified exponential squared kernel. The smaller value the smoother the GP is. The default value is 1.0.
#' @param center A vector of real numbers specifying the centroid parameters in the modified exponential squared kernel. The default value is NULL and set to the center of the grid points: apply(x,2,mean).
#' @param scale A vector of positive numbers specifying the scale parameters in the modified exponential squared kernel. The default value is NULL and set to values such that grid points in a range of (-max_range,max_range) in each dimension.
#' @param max_range A positive real number indicating the maximum range of the grid points to specify the scale parameter. The default value is NULL and set to 6.
#' @return A list of variables representing the simulated curve:
#' \describe{
#'  \item{f}{A vector of real numbers for the simulated curve.}
#'  \item{x}{A matrix of real numbers for the grid points where rows are observations and columns are coordinates.}
#'  \item{work_x}{A matrix of real numbers for the standardized grid points for the simulated curve. It has the same dimension as "x".}
#' }
#'
#' @author Jian Kang <jiankang@umich.edu>
#'
#' @examples
#'library(BayesGPfit)
#'library(lattice)
#'set.seed(1224)
#'dat = list()
#'dat$x = GP.generate.grids(d=2,num_grids = 100)
#'curve = GP.simulate.curve.fast(dat$x,a=0.01,b=0.5,poly_degree=20L)
#'GP.plot.curve(curve,main="Simulated Curve")
#' @export
GP.simulate.curve.fast = function(x,poly_degree,a,b,
                                  center=NULL,scale=NULL,max_range=6){

  x = cbind(x)
  d = ncol(x)

  if(is.null(center)){
    center = apply(x,2,mean)
  }
  c_grids = t(x) - center
  if(is.null(scale)){
    max_grids =pmax(apply(c_grids,1,max),-apply(c_grids,1,min))
    scale=as.numeric(max_grids/max_range)
  }

  work_x = GP.std.grids(x,center=center,scale=scale,max_range=max_range)
  Xmat = GP.eigen.funcs.fast(grids=work_x,
                             poly_degree =poly_degree,
                             a =a ,b=b)
  lambda = GP.eigen.value(poly_degree=poly_degree,a=a,b=b,d=d)
  betacoef = rnorm(ncol(Xmat),mean=0,sd=sqrt(lambda))
  f = Xmat%*%betacoef
  return(list(f=f,x=x,work_x=work_x))
}

#'Graphical representation of one, two, three-dimensional
#'curves
#'@param curve A list object with two elements:
#'\describe{
#'  \item{f}{A vector of real numbers for the curve.}
#'  \item{x}{A matrix of real numbers for the grid points where rows are observations and columns are coordinates.}
#'}
#'@param xlab A character specifying the label of x-axis for 1D, 2D and 3D case. The default value is NULL and set to "x" for 1D case and "x1" for 2D and 3D cases.
#'@param ylab A character specifying the label of y-axis for 1D curve or coords for 2D and 3D case.  The default value is NULL and set to "x2" for 2D and 3D cases.
#'@param zlab A character specifying the label of z-axis only for 3D case. The default value is NULL and set to "x3".
#'@param xlim A vector of two real numbers specifying the range of x-axis for 1D, 2D and 3D case. The default value is NULL and set to range(curve$x[,1]).
#'@param ylim A vector of two real numbers specifying the range of y-axis only for 2D, 3D case. The default value is NULL and set to range(curve$x[,2]).
#'@param zlim A vector of two real numbers specifying the range of z-axis only for 3D case. The default value is NULL and set to range(curve$x[,3]).
#'@param col.regions A vector of RGB colors for 2D and 3D plots. See \link{GP.create.cols}. The default value is NULL and set to GP.create.cols().
#'@param cut An integer specifying the number of colors in 2D and 3D plots. The default value is NULL and set to length(col.regions)-1.
#'@param num_slices An integer specifying the number of slices cutting through the 3rd dimension to show.
#'@param ... All other parameters for plot (1D case) and levelplot (2D and 3D cases).
#'
#'@return NULL for 1D case. An object of class "trellis" for 2D and 3D cases.
#' @author Jian Kang <jiankang@umich.edu>
#'@examples
#'library(BayesGPfit)
#'library(lattice)
#'set.seed(1224)
#'##plot 1D curve
#'x1d = GP.generate.grids(d=1,num_grids = 1000)
#'curve1d = GP.simulate.curve.fast(x1d,a=0.01,b=0.5,
#'                               poly_degree=10L)
#'GP.plot.curve(curve1d,main="Simulated 1D Curve")
#'
#'##plot 2D curve
#'x2d = GP.generate.grids(d=2L,num_grids = 100)
#'curve2d = GP.simulate.curve.fast(x2d,a=0.01,b=0.5,
#'                               poly_degree=10L)
#'GP.plot.curve(curve2d,main="Simulated 2D Curve")
#'
#'##plot 3D curve
#'x3d = GP.generate.grids(d=3,num_grids = 50)
#'curve3d = GP.simulate.curve.fast(x3d,a=0.01,b=0.5,
#'                               poly_degree=10L)
#'GP.plot.curve(curve3d,main="Simulated 3D Curve",num_slices=10,zlim=c(-0.5,0.5))
#'@importFrom lattice levelplot strip.custom
#'@export
GP.plot.curve = function(curve,xlab=NULL,ylab=NULL,
                         zlab=NULL,
                         xlim=NULL,ylim=NULL,zlim=NULL,
                         col.regions=NULL,cut=NULL,
                         num_slices=NULL,...){
  if(ncol(curve$x)==1L){
    od = order(curve$x)
    if(is.null(xlab))
      xlab = "x"
    if(is.null(ylab))
      ylab = "f(x)"

    plot(curve$x[od],curve$f[od],
         xlab=xlab,ylab=ylab,
         xlim=xlim,
         ylim=ylim,...)
  }

  if(ncol(curve$x)==2L){
    if(is.null(col.regions))
      col.regions = GP.create.cols()
    if(is.null(cut))
      cut = length(col.regions)-1L
    if(is.null(xlim))
      xlim = range(curve$x[,1])
    if(is.null(ylim))
      ylim = range(curve$x[,2])
    if(is.null(xlab))
      xlab = "x1"
    if(is.null(ylab))
      ylab = "x2"


    return(levelplot(curve$f~curve$x[,1]+curve$x[,2],
              xlab=xlab,ylab=ylab,xlim=xlim,
              ylim=ylim,col.regions = col.regions,cut=cut,...))

  }

  if(ncol(curve$x)==3L){
    if(is.null(col.regions))
      col.regions = GP.create.cols()
    if(is.null(cut))
      cut = length(col.regions)-1
    if(is.null(xlim))
      xlim = range(curve$x[,1])
    if(is.null(ylim))
      ylim = range(curve$x[,2])
    if(is.null(zlim))
      zlim = range(curve$x[,3])
    if(is.null(xlab))
      xlab = "x1"
    if(is.null(ylab))
      ylab = "x2"
    if(is.null(zlab))
      zlab = "x3"
    if(is.null(num_slices)){
      num_slices = 4L
    }
    x3 =  curve$x[,3]
    uni_x3 = unique(x3)
    uni_x3 = uni_x3[uni_x3>=zlim[1] & uni_x3<=zlim[2]]
    if(length(uni_x3)<num_slices){
      num_slices = length(uni_x3)
    }

    uni_idx = round(seq(1,length(uni_x3),length=num_slices))
    slice_idx = which(is.element(x3,uni_x3[uni_idx]))
    #od = order(x3[slice_idx],decreasing = FALSE)
    #slice_idx = slice_idx[od]
    #x3_print = paste(zlab,"=", round(x3[slice_idx],digits=2L))
    x3 = round(x3[slice_idx],digits=2L)

    return(levelplot(curve$f[slice_idx]~curve$x[slice_idx,1] + curve$x[slice_idx,2] | x3,
              xlab=xlab,ylab=ylab,xlim=xlim,
              ylim=ylim,
              col.regions = col.regions,cut=cut,
              strip = strip.custom(strip.names = TRUE,
                                   var.name = zlab,
                                   strip.levels = TRUE),...))
  }
}

#'Graphical representation of multiple curves in one and two-dimensional
#'curves
#'@param curves A list object of multiple curves and each curve is a list with two elements:
#'\describe{
#'  \item{f}{A vector of real numbers for the curve.}
#'  \item{x}{A matrix of real numbers for the grid points where rows are observations and columns are coordinates.}
#'}
#'@param xlab A character specifying the label of x-axis for 1D, 2D and 3D case. The default value is NULL and set to "x" for 1D case and "x1" for 2D case.
#'@param ylab A character specifying the label of y-axis for 1D curve or coords for 2D and 3D case.  The default value is NULL and set to "x2" for 2D case.
#'@param cols A vector of integer numbers or characters to specify the plot colors for 1D curve. The default value is NULL and set to 1:length(curves).
#'@param lwd A positive number to specify the width of lines for 1D curve.
#'@param type A character specifying what type of plot should be drawn for 1D curve. Possible types are the same as \link{plot}.
#'@param leg_pos A character spaecifying the position of legend for multiple 1D curves. Possible valeus are "topleft", "topright","bottemright","bottemleft".
#'@param xlim A vector of two real numbers specifying the range of x-axis for 1D, 2D and 3D case. The default value is NULL and set to range(curve$x[,1]).
#'@param ylim A vector of two real numbers specifying the range of y-axis only for 2D, 3D case. The default value is NULL and set to range(curve$x[,2]).
#'@param col.regions A vector of RGB colors for 2D and 3D plots. See \link{GP.create.cols}. The default value is NULL and set to GP.create.cols().
#'@param cut An integer specifying the number of colors in 2D  plots. The default value is NULL and set to length(col.regions)-1.
#'@param ... All other parameters for plot (1D case) and levelplot (2D case).
#'
#'@return NULL for 1D case. An object of class "trellis" for 2D and 3D cases.
#'@author Jian Kang <jiankang@umich.edu>
#'@importFrom lattice levelplot
#'@importFrom grDevices rgb
#'@importFrom graphics legend lines plot
#'@importFrom stats quantile rgamma rnorm runif
#'@importFrom utils setTxtProgressBar txtProgressBar
#'@examples
#'library(BayesGPfit)
#'library(lattice)
#'set.seed(1227)
#'dat = list()
#'dat$x = GP.generate.grids(d=2L,num_grids = 100)
#'curve = GP.simulate.curve.fast(dat$x,a=0.01,b=0.5,poly_degree=20L)
#'dat$f = curve$f + rnorm(length(curve$f),sd=1)
#'fast_fit = GP.fast.Bayes.fit(dat$f,dat$x,a=0.01,b=0.5,poly_degree=20L,progress_bar = TRUE)
#'reg_fit = GP.Bayes.fit(dat$f,dat$x,a=0.01,b=0.5,poly_degree=20L,progress_bar = TRUE)
#'curves = list(True = curve,
#'Bayes_fast = fast_fit,
#'Bayes = reg_fit)
#'GP.plot.curves(curves,
#'               main="Comparisons of Bayesian model fitting")
#'
#'@export
GP.plot.curves = function(curves,xlab=NULL,ylab=NULL,
                         cols = NULL, lwd=NULL, type=NULL, leg_pos=NULL,
                         xlim=NULL,ylim=NULL,
                         col.regions=NULL,cut=NULL,...){

  if(ncol(curves[[1]]$x)==1L){

    od = order(curves[[1]]$x)
    if(is.null(cols))
      cols=1:length(curves)
    if(is.null(xlab))
      xlab = "x"
    if(is.null(ylab))
      ylab = "f(x)"
    if(is.null(lwd))
      lwd = 2
    if(is.null(type))
      type = "l"
    if(is.null(leg_pos))
      leg_pos = "topleft"

    plot(curves[[1]]$x[od],curves[[1]]$f[od],
         xlab=xlab,ylab=ylab,
         xlim=xlim,
         ylim=ylim,type=type,col=cols[1],lwd=lwd,...)

    if(length(curves)>1){
      for(i in 2:length(curves)){
        lines(curves[[i]]$x[od],curves[[i]]$f[od],col=cols[i],lwd=lwd)
      }
    }

    nms = names(curves)
    if(is.null(nms)){
      nms = 1:length(curves)
    }
    legend(leg_pos,nms,col=cols,lwd=lwd,lty=1)
  }

  if(ncol(curves[[1]]$x)==2L){
    if(is.null(col.regions))
      col.regions = GP.create.cols()
    if(is.null(cut))
      cut = length(col.regions)-1L
    if(is.null(xlim))
      xlim = range(curves[[1]]$x[,1])
    if(is.null(ylim))
      ylim = range(curves[[1]]$x[,2])
    if(is.null(xlab))
      xlab = "x1"
    if(is.null(ylab))
      ylab = "x2"

    nms = names(curves)
    if(is.null(nms)){
      nms = 1:length(curves)
    }

    x1 = NULL
    x2 = NULL
    group = NULL
    f = NULL
    for(i in 1:length(curves)){
      x1 = c(x1, curves[[i]]$x[,1])
      x2 = c(x2, curves[[i]]$x[,2])
      group = c(group,  rep(nms[i],length=length(curves[[i]]$x[,1])))
      f = c(f,  curves[[i]]$f)
    }

    group = factor(group,levels=nms)

    return(levelplot(f~x1+x2 | group,
                     xlab=xlab,ylab=ylab,xlim=xlim,
                     ylim=ylim,col.regions = col.regions,cut=cut,...))

  }

}


