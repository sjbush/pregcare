# LT 2/12/2021
# integrated likelihood for PREGCARE study

estimatep = function (c1, n1, c2, n2) {
# c1: number of ALT in the control
# n1: total number of reads in the control
# c2: number of ALT in the case
# n2: total number of reads in the case
  
  # joint likelihood  
  f = function(x, p, c1,n1,c2,n2) {
    # x is a vector
    # joint likelihood: sum on log-scale for numerical accuracy
    exp(dbinom(c1, n1, x, log=TRUE) + dbinom(c2, n2, p + (1-p) * x, log=TRUE))
  }
  
  # integrated log-likelihood
  ll = function(p, c1, n1, c2, n2) {
    # naively integrating over (0,1) does not work well
    # find smart upper bound of integration interval to avoid numerical instability
    ub = uniroot(function(x) dbinom(c1,n1, x, log=TRUE) + 100, c(c1/n1, 1))
    #cat('integration upper bound: ', ub$root, '\n')
    log(integrate(f, 0, ub$root, p=p, c1=c1, n1=n1, c2=c2, n2=n2, rel.tol=.Machine$double.eps^0.5)$value)
    #integrate(f, 0, 1, p=p, c1=c1, n1=n1, c2=c2, n2=n2, abs.tol=0)
  }
  
  
  # if no ALTs observed mle is 0 and upper bound is calculated with a mixture of chi square
  
  if (c2 == 0) {
    mle = 0
    lb = 0
    
    # likelihood at mle
    ll.mle = ll(mle, c1, n1, c2, n2)
    # penalty function
    h2 = function(...) {
      (ll.mle - ll(...) - qchisq(0.90, df=1)/2)^1
    }
    fit.ub =uniroot(h2, c(mle, min(1.5*binom.test(c2,n2)$conf.int[2],1)), tol=.Machine$double.eps, c1, n1, c2,n2)
    ub = fit.ub$root
    
  } else {
  # work out good parameter space bounds
  # lower bound 0 is OK
  # use c2/N2 as upper bound. if c2 == 0 use arbitrary small upper bound of 0.01
  
  # find mle
  fit = optimize(ll, c(0,ifelse(c2==0, 0.01, c2/n2)), maximum=TRUE, tol=.Machine$double.eps, c1=c1, n1=n1, c2=c2, n2=n2)
  
  if (!is.finite(fit$objective)) {
    cat('NA: c1=', c1, ', n1=', n1, ', c2=', c2, ', n2=', n2, '\n')
    return(c(vaf=NA, lower=NA, upper=NA))
  }
  
  mle=fit$maximum
  
  # find upper bound of CI
  h = function(...) {
    (fit$objective - ll(...) - qchisq(0.95, df=1)/2)^2
  }
  
  # penalty function
  h2 = function(...) {
    (fit$objective - ll(...) - qchisq(0.95, df=1)/2)^1
  }
  
  # initial try with quadratic penalty: numerically unstable
  #fit.ub = optimize(h, c(mle, min(1.5*binom.test(c2,n2)$conf.int[2],1)), maximum=FALSE, tol=.Machine$double.eps, c1=c1, n1=n1, c2=c2, n2=n2)
  #ub = fit.ub$minimum
  
  # switch to uniroot
  fit.ub =uniroot(h2, c(mle, min(1.5*binom.test(c2,n2)$conf.int[2],1)), tol=.Machine$double.eps, c1, n1, c2,n2)
  ub = fit.ub$root
  
  # find lower bound
  # check whether likelihood at zero is above critical value
  if (ll(mle, c1, n1, c2, n2) - ll(0, c1, n1, c2, n2)> qchisq(0.95, df=1)/2) {
    fit.lb = optimize(h, c(0, fit$maximum), maximum=FALSE, tol=.Machine$double.eps, c1=c1, n1=n1, c2=c2, n2=n2)
    lb = fit.lb$minimum
  } else lb=0
  
  }
  return(c(vaf=mle, lower=lb, upper=ub))
}

# example FAM 57
estimatep(c1=25, n1=64474, c2=4, n2=7762)

# example FAM 30  
estimatep(c1=44, n1=110308, c2=13, n2=41851)
