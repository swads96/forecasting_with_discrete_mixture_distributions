library(scoringutils)

md1 <- distr::UnivarMixingDistribution(distr::Norm(1.5,1), 
                                       distr::Norm(2,2),
                                       mixCoeff = c(.3,.7))

md2 <- distr::UnivarMixingDistribution(distr::Lnorm(2,1),
                                       distr::Norm(2.1,1),
                                       mixCoeff = c(.4,.6))

md2 <- distr::Truncate(md2,lower=0,upper=4)

pmd1 <- function(x) {p(md1)(x)}
pmd2 <- function(x) {p(md2)(x)}


crps_integrand <- function(x,dist,y) {(dist(x) - as.numeric(y <= x))^2}
CRPS <- function(y,dist) {
  int <- integrate(crps_integrand,-Inf,Inf,y,dist=dist)
  return(int$value)
}

CRPS(pmd1,y=4)
CRPS(pmd2,y=4)
ens(lw1,x=4)
CRPS(pmd1+.6*pmd2,y=2.2)

ens <- function(lw1,x) {
  w1 <- boot::inv.logit(lw1)
  w2 <- 1-w1
  p <- function(y) {w1*pmd1(y) + w2*pmd2(y)}
  s <- CRPS(x,p)
  return(s)
}

j <- integrate(dnorm,-Inf,Inf)
test <- optim(fn=ens,par=c(12),x=4)
optimize(f=CRPS,interval=c(-3,3),dist=pmd1)

weight_ensemble <- function(w1,w2,x) {
  ens <- w1*pmd1(x) + w2*pmd2(x)
  return(ens)
}

optim()

crps(y=2.2,family='normal',mean=0,sd=1)

install.packages('scoringRules')

crps(y=2.1,family='mixnorm',m=c(1,2),s=c(1,1),w=c(.5,.5))












bruh(fn=ens,par=0,x=.5,lower=-30,upper=30,method='Brent')


bruh <- function (par, fn, gr = NULL, ..., method = c("Nelder-Mead", 
                                              "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), lower = -Inf, 
          upper = Inf, control = list(), hessian = FALSE) 
{
  fn1 <- function(par) fn(par, ...)
  gr1 <- if (!is.null(gr)) 
    function(par) gr(par, ...)
  method <- match.arg(method)
  if ((any(lower > -Inf) || any(upper < Inf)) && !any(method == 
                                                      c("L-BFGS-B", "Brent"))) {
    warning("bounds can only be used with method L-BFGS-B (or Brent)")
    method <- "L-BFGS-B"
  }
  npar <- length(par)
  con <- list(trace = 0, fnscale = 1, parscale = rep.int(1, 
                                                         npar), ndeps = rep.int(0.001, npar), maxit = 100L, abstol = -Inf, 
              reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, 
              gamma = 2, REPORT = 10, warn.1d.NelderMead = TRUE, type = 1, 
              lmm = 5, factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
  nmsC <- names(con)
  if (method == "Nelder-Mead") 
    con$maxit <- 500
  if (method == "SANN") {
    con$maxit <- 10000
    con$REPORT <- 100
  }
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  if (con$trace < 0) 
    warning("read the documentation for 'trace' more carefully")
  else if (method == "SANN" && con$trace && as.integer(con$REPORT) == 
           0) 
    stop("'trace != 0' needs 'REPORT >= 1'")
  if (method == "L-BFGS-B" && any(!is.na(match(c("reltol", 
                                                 "abstol"), namc)))) 
    warning("method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'")
  if (npar == 1 && method == "Nelder-Mead" && isTRUE(con$warn.1d.NelderMead)) 
    warning("one-dimensional optimization by Nelder-Mead is unreliable:\nuse \"Brent\" or optimize() directly")
  if (npar > 1 && method == "Brent") 
    stop("method = \"Brent\" is only available for one-dimensional optimization")
  lower <- as.double(rep_len(lower, npar))
  upper <- as.double(rep_len(upper, npar))
  res <- if (method == "Brent") {
    if (any(!is.finite(c(upper, lower)))) 
      stop("'lower' and 'upper' must be finite values")
    res <- optimize(function(par) fn(par, ...)/con$fnscale, 
                    lower = lower, upper = upper, tol = con$reltol)
    names(res)[names(res) == c("minimum", "objective")] <- c("par", 
                                                             "value")
    res$value <- res$value * con$fnscale
    c(res, list(counts = c(`function` = NA, gradient = NA), 
                convergence = 0L, message = NULL))
  }
  else .External2(C_optim, par, fn1, gr1, method, con, lower, 
                  upper)
  if (hessian) 
    res$hessian <- .External2(C_optimhess, res$par, fn1, 
                              gr1, con)
  res
}