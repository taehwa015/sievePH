#' @importFrom stats as.formula lm predict sd quantile optim pnorm
#' @importFrom splines2 mSpline
#' @importFrom MASS ginv
NULL
#' Fit the full likelihood proportional hazards model
#'
#' Fit the proportional hazards model with maximum full likelihood estimation. Sieve estimation is used for estimating the baseline hazard function.
#'
#'
#' @param y survival time (> 0).
#' @param d right-censoring indicator, \code{1}: observed; \code{0}: right-censored.
#' @param x p-dimensional covariates matrix.
#'
#' @return \code{smle_ph} returns a list containing the following components:
#' \itemize{
#'   \item \code{Coef}: regression estimator and its inferential results.
#'   \item \code{Cum.hazard}: baseline cumulative hazard function estimates.
#' }
#'
#' @details
#' see Halabi et al., (2024+) for detailed method explanation.
#'
#' @references
#' Halabi et al., (2024+) Sieve maximum full likelihood estimation for the proportional hazards model
#'
#'
#' @examples
#' library(smlePH)
#' set.seed(111)
#' n = 200
#' beta = c(1, -1, 0.5, -0.5, 1)
#' p = length(beta)
#' beta = matrix(beta, ncol = 1)
#' R = matrix(c(rep(0, p^2)), ncol = p)
#' diag(R) = 1
#' mu = rep(0, p)
#' SD = rep(1, p)
#' S = R * (SD %*% t(SD))
#' x = MASS::mvrnorm(n, mu, S)
#' T = (-log(runif(n)) / (2 * exp(x %*% beta)))^(1/2)
#' C = runif(n, min = 0, max = 2.9)
#' y = apply(cbind(T,C), 1, min)
#' d = (T <= C)+0
#' ord = order(y)
#' y = y[ord]; x = x[ord,]; d = d[ord]
#' smle_ph(y = y, d = d, x = x)
#' @export

smle_ph = function(y,
                   d,
                   x)
{
  phfunc_sieve = function(y,
                          d,
                          x,
                          degree,
                          nknots)
  {
    ord = order(y)
    ut = y = y[ord]
    d = d[ord]
    x = as.matrix(as.matrix(x)[ord,])
    if (nknots == 0) {
      knots = numeric(0)
    } else {
      qq = seq(0,1,by=1/(nknots+1))[-c(1,nknots+2)]
      knots = quantile(ut,qq)
    }
    msp = splines2::mSpline(ut,degree=degree,knots=knots,intercept = FALSE)
    n = nrow(x)
    p = ncol(x)
    q = ncol(msp)
    slike = function(theta) {
      beta = theta[1:p]
      gamma = theta[-(1:p)]
      haz = pmax(as.vector(msp%*%gamma),1e-5)
      Haz = cumsum(haz*diff(c(0,ut)))
      xbeta=drop(x%*%beta)
      sum((d*(log(haz)+xbeta) - Haz*exp(xbeta)))
    }

    theta0 = c(rep(0,p), cumsum(rep(1,q)))
    theta = optim(theta0, slike,control = list(fnscale=-1))$par
    coef = theta[1:p]
    gamma = theta[-(1:p)]
    haz = pmax(as.vector(msp%*%gamma),1e-5)
    Haz = cumsum(haz*diff(c(0,ut)))
    list(coef = coef, like = slike(theta), msp = msp, scoef = gamma,
         Haz = Haz, haz = haz, time = y)
  }

  var_func = function(y,
                      d,
                      x,
                      fit,
                      hn)
  {
    slike = function(theta)
    {
      beta = theta[1:p]
      gamma = theta[-(1:p)]
      haz = pmax(as.vector(msp%*%gamma),1e-5)
      Haz = cumsum(haz*diff(c(0,y)))
      xbeta=drop(x%*%beta)
      ((d*(log(haz)+xbeta) - Haz*exp(xbeta)))
    }
    theta = c(fit$coef,fit$scoef)
    msp = fit$msp
    ord = order(y)
    ut = y = y[ord]
    d = d[ord]
    x = as.matrix(as.matrix(x)[ord,])
    n = nrow(x); p = ncol(x); q = length(theta) - p
    cvec = diag(1, p+q)
    A = matrix(0, n, p+q)
    for (i in 1:(p+q)) {
      A[,i] = slike(theta + hn*cvec[i,])
    }
    B = slike(theta)
    Cmat = crossprod((A - B), (A - B))/hn^2
    sqrt(diag(MASS::ginv( Cmat + diag(1e-4, p+q) )))[1:p]
  }

  ord = order(y)
  ut = y = y[ord]
  d = d[ord]
  x = as.matrix(as.matrix(x)[ord,])
  grids = expand.grid(3, 0:(floor(nrow(x)^(1/3))))
  like_val = apply(grids, 1, function(a) phfunc_sieve(y, d, x, a[1], a[2])$like)
  val = -2*like_val + 2*rowSums(grids)*log(log(nrow(x)))
  opt = as.numeric(grids[which.min(val),])
  tmp = phfunc_sieve(y, d, x, opt[1], opt[2])
  tmp$opt = opt
  tmp$se = var_func(y, d, x, tmp,hn=0.1*nrow(x)^(-1/2))

  coef.smr = data.frame("Coefficients" = tmp$coef,
                        "Std.Err"= tmp$se,
                        "T0" = abs(tmp$coef/tmp$se),
                        "p-value" = pnorm(abs(tmp$coef/tmp$se),
                                          lower.tail = FALSE))
  c.haz = data.frame("chaz" = tmp$Haz, "time" = tmp$time)
  res = list("Coef" = round(coef.smr, 3), "Cum.hazard" = c.haz)
  res
}



# roxygen2::roxygenize()
# devtools::build_manual()
# devtools::check(cran=TRUE)
# https://win-builder.r-project.org/upload.aspx