#########################################################################################################################################
####### Non-critical b-d-process 
#########################################################################################################################################

#' Non-critical clone size distribution (exact).
#'
#' @description Exact probability to grow from a clone of size "a" to a clone of size "b" within time "t" according to a non-critical birth-death process.
#'
#' @param lambda proliferation rate 
#' @param delta loss rate
#' @param t time
#' @param a clone size at t=0
#' @param b clone size at t=t
#' @return The probability that a clone of size a grows to size "b" within "t".
#' @references  Bailey, NTJ (1964). The elements of stochastic processes with applications to the natural sciences, Wiley (New York).
#' @examples
#' density.a.b.exact(1, 0, 10, 1, 2)
#' @export

density.a.b.exact <- function(lambda, delta, t, a, b){
  if(a==1){
    if(b==0){
      .alpha(lambda, delta, t)
    }else{
      (1-.alpha(lambda, delta, t))*(1-.beta(lambda, delta, t))*
        .beta(lambda, delta, t)^(b-1)
    }
  }else{
    if(b==0){
      .alpha(lambda, delta, t)^a
    }else{
      sum(sapply(seq(0,min(a,b)), function(j){
        choose(a, j)*choose(a+b-j-1, a-1)*.alpha(lambda, delta, t)^(a-j)*.beta(lambda, delta, t)^(b-j)*(1-.alpha(lambda, delta, t)-
                                                                                                          .beta(lambda, delta, t))^j
      }))
    }
  }
}

.alpha <- function(lambda, delta, t){
  (delta*exp((lambda - delta)*t) - delta)/
    (lambda*exp((lambda - delta)*t) - delta)
}

.beta <- function(lambda, delta, t){
  (lambda*exp((lambda - delta)*t) - lambda)/
    (lambda*exp((lambda-delta)*t) - delta)
}

#' Clone size distribution in a noncritical birth-death process (approximate).
#' @description Probability of a clone of size "a" to grow to size "b" within "t" according to a noncritical linear birth-death process. 
#' @param lambda proliferation rate 
#' @param delta loss rate
#' @param t time
#' @param a clone size at t=0
#' @param b clone size at t=t
#' @param mode "density" if density distribution is to be returned , "cumulative" if cumulative distribution is to be returned. Defaults to "cumulative"
#' @param approx Approximation to be used. Defaults to "highnumbers"; i.e. the distribution is approximated with a gamma distribution if `a` and `b` are large.
#' @details
#' If `approx="highnumbers"`, the function is approximated with a \eqn{\Gamma}-distribution if `a+b>100` and `mode="density"` or if `a+b>10` and `mode="cumulative"`. The \eqn{\Gamma}-distribution is parametrized with 
#' \eqn{shape = \mu^2/\sigma, scale = \sigma/\mu}, where
#' \eqn{\mu = a e^{(\lambda - \delta)t}, \sigma = a \frac{\lambda + \delta}{\lambda - \delta}e^{(\lambda - \delta)t}(e^{(\lambda - \delta)t}-1)}
#' @references Bailey, NTJ (1964). The elements of stochastic processes with applications to the natural sciences, Wiley (New York).
#' @return The probability of growing from size a to size b within t. The Function switches between the exact solution and an approximate solution according to a parametrized gamma distribution.
#' @export
#' 
p.a.b <- function(lambda, delta, t, a, b, mode="cumulative", approx="highnumbers"){
  if(a==0){
    if(mode=="density"){
      if(b==0){
        return(1)
      }else{
        return(0)
      }
    }else{
      return(1)
    }
  }
  if(mode=="density"){
    if(a+b<=100 & approx=="highnumbers"){
      density.a.b.exact(lambda, delta, t, a, b)
    }else{
      ## Approximate with gamma-distribution, parametrized by mean and variance
      mean.g <- a*exp((lambda-delta)*t)
      var.g <- a*(lambda+delta)/(lambda-delta)*exp((lambda-delta)*t)*(exp((lambda-delta)*t)-1)
      dgamma(b, shape=mean.g^2/var.g, scale=var.g/mean.g)
    }
    
  }else if(mode=="cumulative"){
    if(a+b<=10 & approx=="highnumbers"){
      ## cumulative distribution
      sum(sapply(seq(0, b), function(b){
        ## Pab according to Bailey, 1964
        density.a.b.exact(lambda, delta, t, a, b)
      }))

    }else{
      ## Approximate with gamma-distribution, parametrized by mean and variance
      mean.g <- a*exp((lambda-delta)*t)
      var.g <- a*(lambda+delta)/(lambda-delta)*exp((lambda-delta)*t)*(exp((lambda-delta)*t)-1)
      if(mean.g==0 & var.g==0){
        return(1)
      }
      pgamma(b, shape=mean.g^2/var.g, scale=var.g/mean.g)
    }
  }
}


#' Mutation accumulation in a growing tissue
#' 
#' @description Expected number of neutral mutations that are present in at least `n.min` cells at `t.end` in an exponentially growing or contracting tissue.
#' @param lambda proliferation rate 
#' @param delta loss rate
#' @param t.end time
#' @param mu mutation rate per cell division
#' @param n.min minimal clone size at `t.end`; can be a value or a vector
#' @param N0 initial population size
#' @param N final population size
#' @param mode if "approx" the sum is approximated by integration. If "exact" the sum is exactly computed for clone sizes between 1 and 10 but beyond that also approximated.
#' @details The expected number of mutations present in at least `n.min` cells is computed as \cr
#' \eqn{M(n_\mathrm{min})=\sum_{n_\mathrm{min}}^N \mu \lambda \int_0^t e^{(\lambda - \delta)(t-t')} P(1,n_{min},t-t')dt'}, which is approximated to \cr
#' \eqn{M(n_\mathrm{min})\approx \mu \lambda \int_0^t e^{(\lambda - \delta)(t-t')} \frac{ P(1,N,t-t') -  P(1,n_\mathrm{min},t-t')}{\log y(t-t')}dt'}, \cr where
#' \eqn{y(t) = \frac{\lambda e^{(\lambda - \delta)t} - \lambda}{\lambda e^{(\lambda - \delta)t}-\delta}}, if `mode=="approx"` or if \eqn{n_\mathrm{min} \le 10}
#' @references  Bailey, NTJ (1964). The elements of stochastic processes with applications to the natural sciences, Wiley (New York).
#' @return The expected number of mutations present in at least `n.min` cells at `t.end` in an exponentially growing tissue. The function assumes that mutations are continuously acquired at a constant rate.
#' @export

mutations.noncritical.bd <- function(lambda, delta, t.end, mu, n.min, N0=1, N=N, mode="approx"){
  
  if(mode=="exact"){
    
    res <- sapply(n.min, function(n.min){
      if(n.min <=10){
        total <- .approx.count(mu, lambda, delta, n.min=11, n.max=100*max(N,N0), t.end, N0) + .exact.count(t, mu, lambda, delta, 10, t.end, N0)
        if(n.min > 1){
          res <- total - .exact.count(t, mu, lambda, delta, n.min-1, t.end, N0)
        }else{
          res <- total 
        }
        return(res)
      }else{
        .approx.count(mu, lambda, delta, n.min=1, n.max=Inf, t.end, N0) - .approx.count(mu, lambda, delta, n.min=1, n.max=n.min, t.end, N0)
      }
      
    })
    return(res)
    
  }
  ## Compute the number of mutations that are present in at least 1 cell and at most n.min cells
  ## The sum necessary in order to compute the cumulative distribution, is here replaced by integration.
  sapply(n.min, function(n.min){
    
    res <- .approx.count(mu, lambda, delta, n.min=1, n.max=100*max(N0, N), t.end, N0) - .approx.count(mu, lambda, delta, n.min=1, n.max=n.min, t.end, N0) 
    return(res)
  })
  
}

## Exact number of mutations present in at least n cells at time t
.exact.count <- function(t, mu, lambda, delta, n, t.end, N0){
  sum(sapply(seq(1,n), function(n){
    integrand <- function(t, mu, lambda, delta, n){
      mu*lambda*N0*exp((lambda - delta)*t)*(density.a.b.exact(lambda, delta, t.end-t, 1, n)) 
      
    }
    ## total number of mutations acquired during exponential growth that survived:
    res <- integrate(integrand, lower=0, upper=t.end, mu=mu, lambda=lambda, delta=delta, n=n)$value
    res
  }))
}

## Approximate number of mutations present in at least n.min cells at time t

.approx.count <- function(mu, lambda, delta, n.min, n.max, t.end, N0){
  integrand <- function(t, mu, lambda, delta, n.min, n.max){
    mu*lambda*N0*exp((lambda - delta)*t)/log(.beta(lambda, delta, t.end-t))*(density.a.b.exact(lambda, delta, t.end-t, 1, n.max) - 
                                                                               density.a.b.exact(lambda, delta, t.end-t, 1, n.min))
    
  }
  total <- integrate(integrand, lower=0, upper=t.end, mu=mu, lambda=lambda, delta=delta, n.min=n.min, n.max=n.max, rel.tol = .Machine$double.eps^0.1)$value
  total
}

#########################################################################################################################################
####### Critical b-d-process 
#########################################################################################################################################

#' Exact solution to grow from a clone of size a to a clone of size b within t in a critical birth-death process
#'
#' @param lambda proliferation rate 
#' @param a clone size at t=0
#' @param b clone size at t=t
#' @param t time
#' @return The probability to growth from a to b within t.
#' @export

p.ss.exact <- function(lambda, a, b, t){
  
  if(a==0 & b==0){return(1)}
  if(a==0 & b!=0){return(0)}
  p <- lambda*t/(1+lambda*t)
  
  res <- sapply(1:min(a,b), function(k){
    k/b*dbinom(x = k, size = a, p=1-p)*dbinom(x = k, size = b, p=1-p)
  })
  
  if(!is.matrix(res)){
    res <- sum(res)
  }else{
    res <- rowSums(res)
  }
  
  res[b==0] <- p^{a}
  
  res
}


#' Approximate solution to grow from a clone of size a to a clone of size b within t in a critical birth-death process using gamma distribution
#'
#' @param lambda proliferation rate 
#' @param a clone size at t=0
#' @param b clone size at t=t
#' @param t time
#' @param mode either 'density' if density distribution is to be returned or 'cumulative'.
#' @return The approximate probability to grow from a to b within t.
#' @export

p.ss.approx <- function(lambda, a, b, t, mode="density"){
  
  if(a==0 & b==0){return(1)}
  if(a==0 & b!=0){return(0)}
  ## mean and variance to parametrize gamma distribution
  mean=a
  var=2*a*lambda*t
  
  if(mode=="density"){
    return(dgamma(b, shape=mean^2/var, scale=var/mean))
  }else{
    return(pgamma(b, shape=mean^2/var, scale=var/mean))
  }
}

#' Clone size distribution in a critical b-d process. 
#' @description Function to compute the probability to grow from `a` to `b` in a critical b-d process using automatic switching between exact and approximate solution. 
#' @param lambda proliferation rate 
#' @param a clone size at `t=0`
#' @param b clone size at `t=t`; single value or vector
#' @param t time
#' @details The function automatically switches between the exact solution and an approximation with a parametrized \eqn{\Gamma}-distribution at a cutoff criterion of \eqn{a*p*(1-p)>=9 \& b*p*(1-p)>=9}
#' @return The probability to grow from `a` to `b` within `t`
#' @export

p.ss <- function(lambda, a, b, t){
  
  p <- lambda*t/(1+lambda*t)
  res <- c()
  for(j in 1:length(b)){
    if(a*p*(1-p)>=9 & b[j]*p*(1-p)>=9){
      res[j] <- p.ss.approx(lambda, a, b[j], t)
    }else{
      a <- round(a)
      b[j] <- round(b[j])
      
      if(a==1){
        res[j] <- p^(b[j]-1)*(1-p)^2
      }else{
        res[j] <- p.ss.exact(lambda, a, b[j], t)
        
      }
    }
    if(b[j]==0){
      res[j] <- p^(a)
    }
  }
  res
}

#' Neutral mutation accumulation during steady state. 
#' 
#' @param lambda proliferation rate 
#' @param N population size
#' @param mu mutation rate per cell division
#' @param n.min minimal clone size
#' @param t.end time at end point
#' @return The number of mutations that were acquired during steady state and are present in at least `n.min` cells.
#' @export

mutations.during.steady.state <- function(lambda, N, mu, n.min, t.end){
  
  ## Compute the number of mutations that are present in at least 1 cell and at most n cells
  
  if(t.end==0){
    return(0)
  }
  ## mutations are acquired at rate mu*lambda*N and drift according to a critical b-d-process starting with size 1. 
  integrand <- function(t, mu, lambda, n, N, t.end){
    res <- c()
    for(j in 1:length(t)){
      p <- lambda*(t.end-t[j])/(1+lambda*(t.end-t[j]))
      ## for small n sum up
      if(n < 100){
        res[j] <- mu*lambda*N*(1-p - sum(sapply(1:n, function(x){p.ss(lambda, 1, x, t.end-t[j])})))
      }else{
        ## for large n approximate with integration
        res[j] <- mu*lambda*N*(1/log(p)*(p^(N*2-1) - p^(n))*(1-p)^2)
      }
    }
    res
    
  }
  
  ## total number of mutations acquired during exponential growth that survived:
  res <- integrate(integrand, lower=0, upper=t.end, mu=mu, lambda=lambda, n=n.min, N=N, t.end=t.end)$value
  res
}

#########################################################################################################################################
####### Supercritical b-d-process followed by a second phase of either acritical b-d-process or a noncritical b-d-process
#########################################################################################################################################

#' Mutation accumulation during exponential expansion followed by homeostasis. 
#' 
#' @param mu mutation rate per cell division
#' @param N population size
#' @param lambda.exp proliferation rate during expansion
#' @param delta.exp loss rate during expansion
#' @param lambda.ss proliferation and loss rate during homeostasis
#' @param t.end end point (starting from homeostasis)
#' @param b minimal clone size of interest. Number or vector. 
#' @param accuracy.a step size in which mutations accumulated during expansion are evaluated (evaluation runs between 5 and 100\%); defaults to 5\%
#' @param phase return variants from "both" phases, or from "expansion" or "homeostasis" only
#' @return This function returns the approximate number of mutations in clones of at least `b` cells, by first computing the distribution at the transition time within intervals, then averaging the fate of each interval during homeostasis and adding newly acquired mutations. 
#' @export

mutational.burden <- function(mu, N, lambda.exp, delta.exp, lambda.ss, t.end, b, accuracy.a = 0.05, phase = "both"){
  
  ## Compute the VAF distribution after expansion in discretized intervals (a). Sample more densely for large as. 
  a <- 10^seq(0, log10(N)-1, 0.05)
  a <- unique(round(c(a, seq(0.105, 1, accuracy.a)*N)))
  
  ## compute for each interval the expected number of mutations at Tss, the transition point between expansion and homeostasis
  t.ss <- log(N)/(lambda.exp-delta.exp)
  ## cumulative distribution
  m.tss <- sapply(a, function(x){mutations.noncritical.bd(lambda = lambda.exp, delta = delta.exp, t.end = t.ss, mu = mu, n.min = x, N=N)})
  ## subtract mutation count from higher interval to get mutations per bin
  m.tss <- m.tss - c(m.tss[-1],0)
  
  ## no homeostatic phase if measuring directly after exponential expansion
  if(t.end==0){
    return(mutations.noncritical.bd(lambda = lambda.exp, delta = delta.exp, t.end = t.ss, mu = mu, n.min = b, N=N))
  }
  
  ## compute for each interval the expected number of mutations that will have grown larger or equal than b cells at t.end
  ## do this by integration if a is large and by sums else
  
  mutations.from.expansion <- sapply(b, function(b){
    
    res <- histogram.drift(lower.bins.1 = a, n.muts = m.tss, lower.bins.2 = b, N = N,
                           bin.p1 = 1, bin.p2 = 1,
                           lambda = lambda.ss, delta = lambda.ss, t = t.end)
    res
    
  })
  
  ## in addition, need to add new mutations that were acquired during steady state phase
  
  mutations.from.steady.state <- sapply(b, function(b){mutations.during.steady.state(lambda = lambda.ss, N = N, mu = mu, n.min = b, t.end = t.end)})
  
  if(phase=="both"){
    return(mutations.from.steady.state + mutations.from.expansion)
  }else if(phase=="early"){
    return(mutations.from.expansion)
  }else if(phase=="homeostasis"){
    return(mutations.from.steady.state)
  }
  
}


#' Mutation accumulation during exponential expansion followed by second phase of either expansion, homeostasis or decline 
#' 
#' @param mu mutation rate per cell division
#' @param N.1 population size after first phase
#' @param lambda.1 proliferation rate during initial expansion
#' @param delta.1 loss rate during initial expansion
#' @param lambda.2 proliferation rate during second phase
#' @param delta.2 loss rate during second phase
#' @param t.end end point (starting after initial expansion)
#' @param b minimal clone size of interest. Number or vector. 
#' @param accuracy.a step size in which mutations accumulated during expansion are evaluated (evaluation runs between 5 and 100\%); defaults to 5\%
#' @param phase return variants from "both" phases, or from "first" or "second" only
#' @return This function returns the approximate number of mutations in clones of at least `b` cells, by first computing the distribution at the transition time within intervals, then averaging the fate of each interval during homeostasis and adding newly acquired mutations. 
#' @export


mutational.burden.general <- function(mu, N.1, lambda.1, delta.1, lambda.2, delta.2=NULL, t.end, b, accuracy.a = 0.05, phase = "both"){
  
  if(is.null(delta.2)){
    delta.2 <- lambda.2 
    warning("Assuming homeostatic tissue.")
  }
  
  # Initial expansion
  ## Compute the VAF distribution after expansion in discretized intervals (a). Sample more densely for large as. 
  a <- 10^seq(0, log10(N.1)-1, 0.05)
  a <- unique(round(c(a, seq(0.105, 1, accuracy.a)*N.1)))
  
  ## compute for each interval the expected number of mutations at Tss, the transition point between expansion and homeostasis
  t.ss <- log(N.1)/(lambda.1-delta.1)
  ## cumulative distribution
  m.tss <- sapply(a, function(x){mutations.noncritical.bd(lambda = lambda.1, delta = delta.1, t.end = t.ss, mu = mu, n.min = x, N=N.1)})
  ## subtract mutation count from higher interval to get mutations per bin
  m.tss <- m.tss - c(m.tss[-1],0)
  # ## how much more mutations lie on the left as compared to the right border of the bin?
  # bin.factor <- m.tss/c(m.tss[-1],0)
  
  ## no homeostatic phase if measuring directly after exponential expansion
  if(t.end==0){
    return(mutations.noncritical.bd(lambda = lambda.1, delta = delta.1, t.end = t.ss, mu = mu, n.min = b, N=N.1))
  }
  
  
  # II. Second phase
  
  N.final <- N.1*exp((lambda.2 - delta.2)*t.end)
  ## compute for each interval the expected number of mutations that will have grown larger or equal than b cells at t.end
  ## do this by integration if a is large and by sums else
  
  mutations.from.first.phase <- sapply(b, function(b){
    
    res <- histogram.drift(lower.bins.1 = a, n.muts = m.tss, lower.bins.2 = b, N = N.1,
                           bin.p1 = 1, bin.p2 = 1,
                           lambda = lambda.2, delta = delta.2, t = t.end)
    res
    
  })
  
  ## in addition, need to add new mutations that were acquired during second phase
  
  if(lambda.2 == delta.2){ # steady state
    mutations.from.second.phase <- sapply(b, function(b){mutations.during.steady.state(lambda = lambda.2, N = N.1, mu = mu, n.min = b, t.end = t.end)})
  }else{
    mutations.from.second.phase <- sapply(b, function(b){mutations.noncritical.bd(lambda = lambda.2, delta = delta.2, t.end = t.end, mu = mu, n.min = b, N0=N.1, N=N.final, mode="approx")})
  }
  
  if(phase=="both"){
    return(mutations.from.second.phase + mutations.from.first.phase)
  }else if(phase=="first"){
    return(mutations.from.first.phase)
  }else if(phase=="second"){
    return(mutations.from.second.phase)
  }
  
}

#########################################################################################################################################
####### Supercritical birth-death process with selection of an advantageous subclone
#########################################################################################################################################

#' Mutation accumulation during exponential expansion with clonal selection. 
#' 
#' @param mu mutation rate per cell division
#' @param N population size
#' @param lambda proliferation rate 
#' @param delta loss rate 
#' @param t.end end point 
#' @param t.s time point at which selective advantage is acquired.
#' @param s selective advantage
#' @param b minimal clone size of interest. Number or vector. 
#' @return This function returns an approximation by first computing the distribution at the transition time within intervals, then averaging the fate of each interval during homeostasis and adding newly acquired mutations in a scenario where a subpopulation is under positive selection. Returns the number of mutations present in at least `b` cells.
#' @export

mutational.burden.selection.expansion=function(mu,lambda,delta,s,t.s,t.end, b){
  if (s<0){
    message("No selective advantage (s<0)")
  }
  
  ## final tissue size
  N <- exp((lambda - delta)*t.end) - exp((lambda - delta)*(t.end - t.s)) + exp((lambda - delta*s)*(t.end - t.s))
  
  ## size of the selected clone at t.end
  sel.size <- exp((lambda - s*delta)*(t.end-t.s))
  
  
  ## Compute the number of mutations that are present in at least 1 cell and at most n.min cells
  ## The sum necessary in order to compute the cumulative distribution, is here replaced by integration.
  mutations.in.selected.clone.prior.t.s <- sapply(b, function(n.min){
    integrand <- function(t, mu, lambda, delta, n){
      p.mut.in.sel <- exp((lambda - delta)*(t.s - t))/exp((lambda - delta)*t.s)
      if(n==0){
        p.mut.in.sel*mu*lambda*exp((lambda - delta)*t)
      }else{
        p.mut.in.sel*mu*lambda*exp((lambda - delta)*t)*((density.a.b.exact(lambda, delta, t.end-t, 1, N*100) -                                                                                                   
                                                           density.a.b.exact(lambda, delta, t.end-t, 1, n))/log(.beta(lambda, delta, t.end-t)))
        
      }
    }
    ## total number of mutations acquired during exponential growth that survived:
    res <- integrate(integrand, lower=0, upper=t.s, mu=mu, lambda=lambda, delta=delta, n=max(0,n.min-sel.size))$value
    return(res)
  })
  
  mutations.not.in.selected.clone.prior.t.s <- sapply(b, function(n.min){
    integrand <- function(t, mu, lambda, delta, n){
      p.mut.in.sel <- exp((lambda - delta)*(t.s - t))/exp((lambda - delta)*t.s)
      (1-p.mut.in.sel)*mu*lambda*exp((lambda - delta)*t)/log(.beta(lambda, delta, t.end-t))*(density.a.b.exact(lambda, delta, t.end-t, 1, 100*N) - 
                                                                                               density.a.b.exact(lambda, delta, t.end-t, 1, n))
      
    }
    
    ## total number of mutations acquired during exponential growth that survived:
    res <- integrate(integrand, lower=0, upper=t.s, mu=mu, lambda=lambda, delta=delta, n=n.min)$value
    return(res)
  })
  
  
  
  mutations.before.t.s <- mutations.in.selected.clone.prior.t.s + mutations.not.in.selected.clone.prior.t.s
  
  
  mutations.from.selected.clone.after.t.s <- sapply(b, function(n.min){
    mutations.noncritical.bd(lambda = lambda, delta = delta *s, t.end = t.end - t.s, 
                             mu = mu, n.min = n.min, N = sel.size)})
  mutations.from.founder.clone.after.t.s <- exp((lambda - delta)*t.s)*sapply(b, function(n.min){
    mutations.noncritical.bd(lambda = lambda, delta = delta, t.end = t.end - t.s, 
                             mu = mu, n.min = n.min, N = exp((lambda - delta)*(t.end - t.s)), mode="exact")})
  
  mutations.before.t.s + mutations.from.selected.clone.after.t.s + mutations.from.founder.clone.after.t.s
  
}

#########################################################################################################################################
####### Suprecritical b-d-process followed by critical b-d-process with selection of an advantageous subclone; subclonal mutation acquired during homeostasis
#########################################################################################################################################

#' Mutation accumulation during exponential expansion followed by homeostasis with clonal selection. 
#' 
#' @param mu mutation rate per cell division
#' @param N population size
#' @param lambda.exp proliferation rate during expansion
#' @param delta.exp loss rate during expansion
#' @param lambda.ss proliferation and loss rate during homeostasis
#' @param t.end end point (starting from homeostasis)
#' @param t.s time point at which selective advantage is acquired.
#' @param s selective advantage
#' @param b minimal clone size of interest. Number or vector. 
#' @param accuracy.a step size in which mutations accumulated during expansion are evaluated between 5 and 100\%; defaults to 5\%
#' @param min.clone.size the lower detection limit for selected clones
#' @return This function returns an approximation by first computing the distribution at the transition time within intervals, then averaging the fate of each interval during homeostasis and adding newly acquired mutations in a scenario where a subpopulation is under positive selection. Returns the number of mutations present in at least `b` cells.
#' @export

mutational.burden.with.selection <- function(mu, N, lambda.exp, delta.exp, lambda.ss, t.end, t.s, s, b, accuracy.a = 0.05, min.clone.size = 0.05){
  
  ## initialize mutation count
  mutations.at.t.end <- 0
  
  ## compute the relative size of the selected clone at t.end
  f.sel <- exp((lambda.ss - s*lambda.ss)*(t.end-t.s))/N
  
  ## in case the selected clone took over, add an offset of the number of mutations in the founder cell to the solution of mutations during steady state (it's again a neutrally expanding clone)
  if(f.sel>1){
    t.ss <- log(N*0.5)/(lambda.ss*(1-s))
    ## the mutation rate is per division, i.e., each daughter cell receives mu/2 mutations. However, the founder cell is a surviving lineage and hence divided with rate 2lambda during homeostasis, to compensate for the loss of the dying lineages
    mutations.in.founder.cell <- (log(N)/(lambda.exp - delta.exp)*mu/2) + t.s*lambda.ss*mu
    mutations.at.t.end <- mutational.burden(mu, N, 1, s, lambda.ss, t.end-t.s-t.ss, b) + mutations.in.founder.cell
    return(mutations.at.t.end)
  }
  
  ## compute when the selected clone reaches 5% as we cannot resolve smaller subclones
  t.min.clone.size <- log(min.clone.size*N)/(lambda.ss - s*lambda.ss)
  ## if 5% is reached after t end, just take the predicted output at t.end according to homeostatic turnover and neglect expansion of the selected clone
  
  if(t.min.clone.size > (t.end - t.s)){
    mutations.at.t.end <- mutational.burden(mu, N, lambda.exp, delta.exp, lambda.ss, t.end, b, accuracy.a = accuracy.a)
    return(mutations.at.t.end)
  }
  
  ## else if 5% is reached before, take the distribution at t.s as the input for selection and neglect mutations acquired in the founder cell population after t.s. We thus have two contributions:
  ## selected clone + subclonal mutations acquired during expansion of the latter
  ## founder population:drift of mutations acquired prior to t.s
  
  ## The drift of the founder cell population is non-trivial, due to the non-exponential decay in a competing system.
  ## To assess it, we approximate the loss with exponential decay, requiring the same number of death events in the founder cell population 
  
  cell.states <- .forward_dynamics(N = N, init.cells = c(1, 0), lambda.ss = lambda.ss, delta.ss = lambda.ss, lambda.exp = lambda.exp,
                                   delta.exp = delta.exp, s = c(1,s), t.s = c(0,t.s), mother.daughter = matrix(c(1,2),nrow=1), t = seq(0, t.end, length.out = 1000))
  
  delta.founder <- .approximate.delta(lambda.ss, N, t.end - t.s, cell.states[which(t.end <= cell.states[,1])[1],4] -
                                        cell.states[which(cell.states[,1] >= t.s)[1], 4])
  
  ## Next, compute the mutation distribution at t.s, i.e. the time point at which the driver mutation is acquired. Do this in a discretized way, as before. 
  
  if((log10(N)-1)>2){
    a <- c(1, 10, 25, 50, 75, 10^seq(2, log10(N)-1, 0.05))
  }else{
    a <- c(1, 10, 25, 50, 75, 100)
  }
  a <- sort(unique(round(c(a, seq(0.105, 1, accuracy.a)*N))))
  
  ## compute the mutational burden at t.s
  mutations.at.t.s <- mutational.burden(mu=mu, N=N, lambda.exp=lambda.exp, delta.exp=delta.exp, lambda.ss=lambda.ss, t.end=t.s, b=a, accuracy.a = accuracy.a)
  ## mutations per bin (cumulative --> discrete)
  mutations.at.t.s <- mutations.at.t.s - c(mutations.at.t.s[-1],0)
  
  ## Now, every mutation can be either present in the founder population, the selected population or both. Its fate is determined by these cases.
  ## For each bin at ts, the probability that the mutation is present in the selected clone, reads a/N. In this case, it will be present in all selected cells + putatively in a subset of founder cells
  
  ## take the average from the lower and upper border
  mutations.at.t.end <- sapply(b, function(b){
    
    ## if b is bigger than the selected daughter clones, it may contain mutations present in both the selected clone and the founder population (prob.this.combination)
    
    if(b >= f.sel*N){ 
      
      prob <- probability.this.combination(bin.size = a, clone.size.mother = N, n.daughters.present = 1, n.daughters.absent = 0)
      prob.upper <- probability.this.combination(bin.size = c(a[-1], 2*N), clone.size.mother = N, n.daughters.present = 1, n.daughters.absent = 0)
      
      # if the mutation ends up in the daughter clones, drift accounts only for the remaining size (b - total.size.of all.daughters.in.comb)
      res <- histogram.drift(lower.bins.1 = a - 1, n.muts = mutations.at.t.s, lower.bins.2 = round(b - f.sel*N), N = N,
                             bin.p1 = prob, bin.p2 = prob.upper,
                             lambda = lambda.ss, delta = delta.founder, t = t.end-t.s) +
        histogram.drift(lower.bins.1 = a, n.muts = mutations.at.t.s, lower.bins.2 = b, N = N,
                        bin.p1 = 1-prob, bin.p2 = 1-prob.upper,
                        lambda = lambda.ss, delta = delta.founder, t = t.end-t.s)
      
      res
    }else{ ## mutations present in at least b cells, where b <= f.sel. --> cumulative distribution from founder cell population +
      
      prob <- probability.this.combination(bin.size = a, clone.size.mother = N, n.daughters.present = 1, n.daughters.absent = 0)
      prob.upper <- probability.this.combination(bin.size = c(a[-1], 2*N), clone.size.mother = N, n.daughters.present = 1, n.daughters.absent = 0)
      
      # if the mutation ends up in the daughter clones, drift accounts only for the remaining size (b - total.size.of all.daughters.in.comb)
      res <- histogram.drift(lower.bins.1 = a - 1, n.muts = mutations.at.t.s, lower.bins.2 = 0, N = N,
                             bin.p1 = prob, bin.p2 = prob.upper,
                             lambda = lambda.ss, delta = delta.founder, t = t.end-t.s) +
        histogram.drift(lower.bins.1 = a, n.muts = mutations.at.t.s, lower.bins.2 = b, N = N,
                        bin.p1 = 1 - prob, bin.p2 = 1 - prob.upper,
                        lambda = lambda.ss, delta = delta.founder, t = t.end-t.s)
      
      res
    }
  })
  
  ## Finally, we need to add up new mutations acquired during the expansion of the selected clone
  
  mutations.from.selected.clone <- mutations.noncritical.bd(lambda.ss, lambda.ss*s, t.end-t.s, mu, b, N=N)
  
  if(lambda.ss == delta.founder){
    mutations.from.founder.clone <- sapply(b, function(b){mutations.during.steady.state(lambda.ss, N, mu, b, t.end - t.s)})
  }else{
    mutations.from.founder.clone <- mutations.noncritical.bd(lambda.ss, delta.founder, t.end - t.s, mu, b, N0=N, N=(1-f.sel)*N)
  }
  ## truncate at b= f.sel*N
  mutations.from.selected.clone[b>f.sel*N] <- 0
  mutations.at.t.end <- mutations.at.t.end + mutations.from.selected.clone + mutations.from.founder.clone
  
  mutations.at.t.end
}



#' Mutation accumulation during exponential expansion followed by a second phase with clonal selection. 
#' 
#' @param mu mutation rate per cell division
#' @param N population size during homeostasis in absence of CH
#' @param lambda.exp proliferation rate during expansion
#' @param delta.exp loss rate during expansion
#' @param lambda.ss proliferation and loss rate during homeostasis
#' @param t.end end point (starting from homeostasis)
#' @param t.s time point at which selective advantage is acquired.
#' @param s selective advantage
#' @param b minimal clone size of interest. Number or vector. 
#' @param accuracy.a step size in which mutations accumulated during expansion are evaluated between 5 and 100\%; defaults to 5\%
#' @param min.clone.size the lower detection limit for selected clones
#' @return This function returns an approximation by first computing the distribution at the transition time within intervals, then averaging the fate of each interval during homeostasis and adding newly acquired mutations in a scenario where a subpopulation is under positive selection. Returns the number of mutations present in at least `b` cells.
#' @export

mutational.burden.with.selection.no.size.compensation <- function(mu, N, lambda.exp, delta.exp, lambda.ss, t.end, t.s, s, b, accuracy.a = 0.05, min.clone.size = 0.05){
  
  ## initialize mutation count
  mutations.at.t.end <- 0
  
  ## compute the relative size of the selected clone at t.end, assuming that the mutant clone adds on top of the normal cells
  N.final <- N + exp((lambda.ss - s*lambda.ss)*(t.end-t.s))
  f.sel <- exp((lambda.ss - s*lambda.ss)*(t.end-t.s))/N.final
  ## if the clone does not exceed the minimal clone size, just take the predicted output at t.end according to homeostatic turnover and neglect expansion of the selected clone
  if(f.sel < min.clone.size){
    mutations.at.t.end <- mutational.burden(mu, N, lambda.exp, delta.exp, lambda.ss, t.end, b, accuracy.a = accuracy.a)
    return(mutations.at.t.end)
  }
  
  ## else if the minimal clone size is reached before, take the distribution at t.s as the input for selection. We thus have two contributions:
  ## selected clone + subclonal mutations acquired during expansion of the latter
  ## founder population:drift of mutations acquired prior to t.s and after t.s
  
  ## Compute the mutation distribution at t.s, i.e. the time point at which the driver mutation is acquired. Do this in a discretized way, as before. 
  
  
  if((log10(N)-1)>2){
    a <- c(1, 10, 25, 50, 75, 10^seq(2, log10(N)-1, 0.05))
  }else{
    a <- c(1, 10, 25, 50, 75, 100)
  }
  a <- sort(unique(round(c(a, seq(0.105, 1, accuracy.a)*N))))
  
  ## compute the mutational burden at t.s
  mutations.at.t.s <- mutational.burden(mu=mu, N=N, lambda.exp=lambda.exp, delta.exp=delta.exp, lambda.ss=lambda.ss, t.end=t.s, b=a, accuracy.a = accuracy.a)
  ## mutations per bin (cumulative --> discrete)
  mutations.at.t.s <- mutations.at.t.s - c(mutations.at.t.s[-1],0)
  
  ## Now, every mutation can be either present in the founder population, the selected population or both. Its fate is determined by these cases.
  ## For each bin at ts, the probability that the mutation is present in the selected clone, reads a/N. In this case, it will be present in all selected cells + putatively in a subset of founder cells
  
  mutations.at.t.end <- sapply(b, function(b){
    
    ## if b is bigger than the selected daughter clones, it may contain mutations present in both the selected clone and the founder population (prob.this.combination)
    
    if(b >= f.sel*N.final){ 
      
      prob <- probability.this.combination(bin.size = a, clone.size.mother = N, n.daughters.present = 1, n.daughters.absent = 0)
      prob.upper <- probability.this.combination(bin.size = c(a[-1], 2*N), clone.size.mother = N, n.daughters.present = 1, n.daughters.absent = 0)
      
      # if the mutation ends up in the daughter clones, drift accounts only for the remaining size (b - total.size.of all.daughters.in.comb)
      res <- histogram.drift(lower.bins.1 = a - 1, n.muts = mutations.at.t.s, lower.bins.2 = round(b - f.sel*N.final), N = N,
                             bin.p1 = prob, bin.p2 = prob.upper,
                             lambda = lambda.ss, delta = lambda.ss, t = t.end-t.s) +
        histogram.drift(lower.bins.1 = a, n.muts = mutations.at.t.s, lower.bins.2 = b, N = N,
                        bin.p1 = 1-prob, bin.p2 = 1-prob.upper,
                        lambda = lambda.ss, delta = lambda.ss, t = t.end-t.s)
      
      res
    }else{ ## mutations present in at least b cells, where b <= f.sel. --> cumulative distribution from founder cell population +
      
      prob <- probability.this.combination(bin.size = a, clone.size.mother = N, n.daughters.present = 1, n.daughters.absent = 0)
      prob.upper <- probability.this.combination(bin.size = c(a[-1], 2*N), clone.size.mother = N, n.daughters.present = 1, n.daughters.absent = 0)
      
      # if the mutation ends up in the daughter clones, drift accounts only for the remaining size (b - total.size.of all.daughters.in.comb)
      res <- histogram.drift(lower.bins.1 = a - 1, n.muts = mutations.at.t.s, lower.bins.2 = 0, N = N,
                             bin.p1 = prob, bin.p2 = prob.upper,
                             lambda = lambda.ss, delta = lambda.ss, t = t.end-t.s) +
        histogram.drift(lower.bins.1 = a, n.muts = mutations.at.t.s, lower.bins.2 = b, N = N,
                        bin.p1 = 1 - prob, bin.p2 = 1 - prob.upper,
                        lambda = lambda.ss, delta = lambda.ss, t = t.end-t.s)
      
      res
    }
  })
  
  
  ## Finally, we need to add up new mutations acquired during the expansion of the selected clone
  
  mutations.from.selected.clone <- mutations.noncritical.bd(lambda.ss, lambda.ss*s, t.end-t.s, mu, b, N=f.sel*N.final)
  
  mutations.from.founder.clone <- sapply(b, function(b){mutations.during.steady.state(lambda.ss, N, mu, b, t.end - t.s)})
  ## truncate at b= f.sel*N
  mutations.from.selected.clone[b>f.sel*N.final] <- 0
  
  mutations.at.t.end <- mutations.at.t.end + mutations.from.selected.clone + mutations.from.founder.clone
  
  mutations.at.t.end
}

#########################################################################################################################################
####### Histogram evolution from t1 to t2
#########################################################################################################################################

#' Computes the expected histogram of variants in a clone of interest at time t, given a VAF histogram at time zero
#' 
#' @param lower.bins.1 vector of bin sizes of the histogram at t0
#' @param n.muts vector of mutation counts per bin at t0
#' @param bin.p1 vector of probabilities that the variants in each bin will remain in the clone of interest
#' @param bin.p2 as `bin.p1` but for the upper border of the bin
#' @param lower.bins.2 vector of clone sizes at t 
#' @param N the number of cells in the system
#' @param lambda the division rate
#' @param delta the loss rate
#' @param t the time point of evaluation
#' @return The cumulative VAF histogram at time t
#' @export

histogram.drift <- function(lower.bins.1, n.muts, bin.p1=1, bin.p2, lower.bins.2, N, lambda, delta, t){
  
  ## the clone size at the upper border of the bin
  upper.bins <- c(lower.bins.1[-1], 2*N)
  
  # get the weights for the two bin borders
  bin.factor <- n.muts/c(n.muts[-1],0)
  bin.factor[is.infinite(bin.factor)] <- 10^8
  bin.factor[is.na(bin.factor)] <- 1
  
  if(lambda == delta){
    res <- .critical.drift(hist.1 = rbind(lower.bins.1, upper.bins, bin.factor, bin.p1, bin.p2), n.muts, lower.bins.2 = lower.bins.2, lambda =  lambda, t = t)
  }else{
    res <- .noncritical.drift(hist.1 = rbind(lower.bins.1, upper.bins, bin.factor, bin.p1, bin.p2), n.muts, lower.bins.2 = lower.bins.2, lambda =  lambda, delta = delta, t = t)
  }
  
  res
}

.critical.drift <- function(hist.1, n.muts, lower.bins.2,  lambda, t){
  res <- sapply(lower.bins.2, function(b){
    ## weighted average of the probability to drift from the lower and upper boundary of the interval to a size of at least b
    weighted.average <- apply(hist.1, 2, function(x){
      ## from critical b-d-process
      p <- lambda*t/(1+lambda*t)
      ## mean and variance from critical b-d- process to parametrize gamme distribution
      ## mean = a, variance = 2a lambda*t
      ## translates to shape = mean^2/var, scale=var/mean
      if(b>20){ # approximate with gamma distribution if b > 20
        ## lower border of bin
        mean=x[1]
        var=2*x[1]*lambda*t
        ## upper border of the bin
        upper.mean=x[2]
        upper.var=2*x[2]*lambda*t
        
        if(x[1]==0){
          mean.prob.b <- sum(c(0,
                               x[5]*(1 - pgamma(b-1, shape=upper.mean^2/upper.var, scale=upper.var/upper.mean))))/(x[3]+1)
          
        }else{
          ## weighted average of cumulative distribution
          mean.prob.b <- sum(c(x[3]*x[4]*(1-  pgamma(b-1, shape=mean^2/var, scale=var/mean)),
                               x[5]*(1 - pgamma(b-1, shape=upper.mean^2/upper.var, scale=upper.var/upper.mean))))/(x[3]+1)
        }
      }else{
        if(b == 0){
          ## weighted average of cumulative distribution
          mean.prob.b <- sum(c(x[3]*x[4]*1),
                             x[5]*1)/(x[3]+1)
          
        }else{
          ## weighted average of cumulative distribution
          mean.prob.b <- sum(c(x[3]*x[4]*(1-sum(sapply(0:(b-1), function(b){p.ss.exact(lambda, x[1], b, t)}))),
                               x[5]*(1-sum(sapply(0:(b-1), function(b){p.ss.exact(lambda, x[2], b, t)})))))/(x[3]+1)
          
        }
      }
      mean.prob.b
      
    })
    
    sum(n.muts*weighted.average)
  })
  res
}

.noncritical.drift <- function(hist.1, n.muts, lower.bins.2,  lambda, delta, t){
  res <- sapply(lower.bins.2, function(b){
    
    weighted.average <- apply(hist.1, 2, function(x){
      
      if(b==0){
        return( (x[3]*x[4] + x[5])/(x[3] + 1))
      }
      
      (x[3]*x[4]*(1-p.a.b(a=x[1], b=b-1, lambda=lambda, delta=delta, t=t)) +
         x[5]*(1-p.a.b(a=x[2], b=b-1, lambda=lambda, delta=delta, t=t)) 
      )/(x[3] + 1)
      
    })
    sum(n.muts*weighted.average)
  })
  res
}


#' Probability of a variant present in a given bin size ends up in a particular combination of selected daughters if a driver is acquired in a random cell of the mother clone 
#' 
#' @param bin.size vector of bin sizes the variant are present in
#' @param clone.size.mother the size of the mother clone
#' @param n.daughters.present the number of daughters the variant ends up in
#' @param n.daughters.absent the number of daughters the variant does not end up in
#' @return Computes the probability that a given variant that is present in `bin.size` cells of the mother clone ends up in `n.daughters.present` daughter clones, but not in the remaining `n.daughters.absent` if the cells giving rise to the selected daughters are randomly sampled.
#' @export

probability.this.combination <- function(bin.size, clone.size.mother, n.daughters.present, n.daughters.absent){
  
  p <- sapply(bin.size, function(a){
    min(1,a/clone.size.mother)^(n.daughters.present)*(1 - min(1, a/clone.size.mother))^(n.daughters.absent)
  })
  p 
  
}
#########################################################################################################################################
####### Multi-clone problem
#########################################################################################################################################

#' Mutation accumulation during exponential expansion followed by homeostasis. 
#' 
#' @param mu mutation rate per cell division
#' @param N population size
#' @param lambda.exp proliferation rate during expansion
#' @param delta.exp loss rate during expansion
#' @param lambda.ss proliferation and loss rate during homeostasis
#' @param t.end end point (starting from homeostasis)
#' @param t.s vector of time points at which selective advantages are acquired.
#' @param s vector of selective advantages associated with driver mutations
#' @param b minimal clone size of interest. Number or vector. 
#' @param mother.daughter a matrix containing the mother (1st column) - daughter (2nd column) relationships between the subclones 
#' @return This function returns an approximation by first computing the distribution at the transition time within intervals, then averaging the fate of each interval during homeostasis and adding newly acquired mutations in a scenario with 2 nested clonal selections (clone starts growing at t.s >t.ss and grows with a selective advantage s; clone 2 starts. Returns the number of mutations present in at least b cells
#' @export

mutational.burden.multiclone <- function(mu, N, lambda.exp, delta.exp, lambda.ss, t.end, t.s, s, mother.daughter, b, min.clone.size = 0.05, accuracy.a = 0.05){
  
  # add column names to the mother-daughter relationship matrix
  colnames(mother.daughter) <- c("M", "D")
  
  # order clones by time point of appearance
  clone.order <- data.frame(new.id = 1:length(t.s), 
                            old.id = order(t.s))
  
  mother.daughter <- rbind(c(1,1), mother.daughter)
  mother.daughter <- apply(mother.daughter, 2, function(x){
    sapply(x, function(y){
      clone.order[clone.order$old.id == y,"new.id"]
    })
  })
  mother.daughter <- mother.daughter[order(mother.daughter[,"D"]),,drop=F]

  s <- s[clone.order$new.id]
  t.s <- sort(t.s)
  
  # analyze the population dynamics and identify the time points at which individual clones peak
  # initiate the system with N normal cells and 0 mutant cells. The system also returns the total number of death events, which will be needed below to compute drift during contraction
  cell.states <- .forward_dynamics(N = N, init.cells = c(1, rep(0, length(s)-1)), lambda.ss = lambda.ss, delta.ss = lambda.ss, lambda.exp = lambda.exp,
                                   delta.exp = delta.exp, s = s, t.s = t.s, mother.daughter = mother.daughter, t = seq(0, t.end, length.out = 1000))
  final.sizes <- cell.states[nrow(cell.states),1+(1:length(s))]
  
  # neglect clones smaller than the minimal clone size; only filter daughter cells; do iteratively!
  
  to.remove <- c()
  for(clone in mother.daughter[,"D"]){
    daughters.this.clone <- setdiff(.get.progeny(mother.daughter, clone), clone)
    total.size <- sum(cell.states[nrow(cell.states),c(clone, daughters.this.clone) + 1])
    # if the total size of this clone is smaller than the minimal clone size, remove the clone and all it's daughters
    if(total.size/N < min.clone.size){
      to.remove <- unique(c(to.remove, clone, daughters.this.clone))
    }
  }

  if(length(to.remove)>0){
    s <- s[-to.remove]
    t.s <- t.s[-to.remove]
    final.sizes <- final.sizes[-to.remove]
    cell.states <- cell.states[,-c(to.remove + 1, (length(s)+2+to.remove)),drop=F]
    mother.daughter <- mother.daughter[-which(mother.daughter[,"M"] %in% to.remove |
                                                mother.daughter[,"D"] %in% to.remove),,drop=F]
  }
  if(nrow(mother.daughter)==1){
    return(mutational.burden(mu, N, lambda.exp, delta.exp, lambda.ss, t.end, b, accuracy.a = accuracy.a))
  }else{

    # rename daughters in order of appearance
    id.conversion <- data.frame(new.id = 1:length(t.s), 
                                old.id = mother.daughter[,"D"])
    
    mother.daughter <- apply(mother.daughter, 2, function(x){
      sapply(x, function(y){
        id.conversion[id.conversion$old.id == y,"new.id"]
      })
    })
    
    cell.states <- .forward_dynamics(N = N, init.cells = c(1, rep(0, length(s)-1)), lambda.ss = lambda.ss, delta.ss = lambda.ss, lambda.exp = lambda.exp,
                                     delta.exp = delta.exp, s = s, t.s = t.s, mother.daughter = mother.daughter, t = seq(0, t.end, length.out = 1000))
    final.sizes <- cell.states[nrow(cell.states),1+(1:length(s))]
  }

  t.peaks <- cell.states[,1][apply(cell.states, 2, function(x){which.max(round(x/10))})[seq(2,1+length(s))]]
  
  t.peaks <- t.peaks[t.peaks>0]
  # we evaluate the problem step-wise, pausing at the birth dates of each clone and at the time point at which they peak
  time.points.of.interest <- sort(c(t.s, t.peaks))
  time.points.of.interest <- setdiff(time.points.of.interest, c(0, t.end))
  # store the site frequency spectra of mutations generated in each clone in a matrix
  muts.tend <- matrix(0, nrow=length(t.s), ncol=length(b))
  # implement a function that returns the evolution of mutations already present in a cell population if one daughter cell becomes positively selected
  
  # discretize the site frequency spectrum
  if((log10(N)-1)>2){
    a <- c(1, 10, 25, 50, 75, 10^seq(2, log10(N)-1, 0.05))
  }else{
    a <- c(1, 10, 25, 50, 75, 100)
  }
  a <- sort(unique(round(c(a, seq(0.105, 1, accuracy.a)*N))))
  #a <- sort(unique(c(a, b)))
  
  clones <- 1:length(s)
  
  ## initialize new.muts
  new.muts <- matrix(0, nrow=length(t.s), ncol=length(a)) # define a earlier, fine-grained histogram for computation
  new.muts[1,] <- mutational.burden(mu, N, lambda.exp, delta.exp, lambda.ss, time.points.of.interest[1], a)
  for(timepoint in 1:(length(time.points.of.interest))){

    lower.t <- time.points.of.interest[timepoint]
    upper.t <- c(time.points.of.interest[-1], t.end)[timepoint]
    
    # new mutations acquired between lower.t and upper.t
    # iterate through all clones
    for(clone in clones){
      # only proceed if clone exists
      if(t.s[clone] > lower.t){next}
      
      ## subtract mutation count from higher interval to get mutations per bin
      new.muts[clone,] <- new.muts[clone,] - c(new.muts[clone,-1],0)
      new.muts[clone,] <- replace(new.muts[clone,], new.muts[clone,] < 0, 0)
      # I. fate of mutations already present in the system
      clone.size.now <- cell.states[which(cell.states[,1] >= lower.t)[1], 1 + clone]
      clone.size.upper <- cell.states[which(cell.states[,1] >= upper.t)[1], 1 + clone]
      clone.size.tend <- final.sizes[clone]
      # if the clone expands, set its delta according to the selective advantage; otherwise approximate by linear b-d process
      if(clone.size.upper > clone.size.now){
        delta.this.clone.until.upper <- lambda.ss * s[clone]
        # our computation does an exponential growth approximation, but we model the population overall with a competition model. Hence, if the clone is already in the competition phase, we would overshoot if not correcting the time span that accounts for exponential expansion only.
        time.span.upper <- log(clone.size.upper/clone.size.now)/(lambda.ss - delta.this.clone.until.upper) 
        
      }else if(clone.size.upper == clone.size.now){
        delta.this.clone.until.upper <- lambda.ss
        time.span.upper <- upper.t - lower.t
      }else{
        delta.this.clone.until.upper <- .approximate.delta(lambda.ss, clone.size.now, upper.t - lower.t, cell.states[which(upper.t <= cell.states[,1])[1],1 + length(s) + clone] -
                                                             cell.states[which(cell.states[,1] >= lower.t)[1], 1 + length(s) + clone])
        time.span.upper <- upper.t - lower.t
      }
      if(clone.size.tend > clone.size.now){
        delta.this.clone.until.tend <- lambda.ss * s[clone]
        # our computation does an exponential growth approximation, but we model the population overall with a competition model. Hence, if the clone is already in the competition phase, we would overshoot if not correcting the time span that accounts for exponential expansion only.
        time.span.tend <- log(clone.size.tend/clone.size.now)/(lambda.ss - delta.this.clone.until.tend) 
        
      }else if(clone.size.tend == clone.size.now){
        delta.this.clone.until.tend <- lambda.ss
        time.span.tend <- t.end - lower.t
      }else{
        delta.this.clone.until.tend <- .approximate.delta(lambda.ss, clone.size.now, t.end - lower.t, cell.states[nrow(cell.states),1 + length(s) + clone] -
                                                            cell.states[which(cell.states[,1] >= lower.t)[1], 1 + length(s) + clone])
        time.span.tend <- t.end - lower.t
      }
      
      total.size.of.all.daughters.this.clone <- sum(final.sizes[.get.progeny(mother.daughter, clone)])
      # for each bin size, compute the probability that the mutation ends up in the newly founded daughter cell 
      daughters.this.clone <- setdiff(mother.daughter[mother.daughter[,"M"]==clone,"D"], clone)
      # take only the daughters that are born after the current lower.t
      #daughters.this.clone <- setdiff((1:length(t.s))[t.s[daughters.this.clone] >= lower.t], clone)
      daughters.this.clone <- setdiff(intersect(daughters.this.clone, (1:length(t.s))[t.s >= lower.t]), clone)
      combinations.of.daughters <- .clonal.combinations(daughters.this.clone)
      if(length(daughters.this.clone) > 0){
        muts.tend[clone,] <- muts.tend[clone,] + sapply(b, function(b){
         # if(b >= total.size.of.all.daughters.this.clone){return(0)}
          # expected number of mutations in clone size b if present in the combination of daughters
          res.daughters <- sum(unlist(lapply(combinations.of.daughters,function(comb){
            # total size of the daughter subclones
            
            total.size.of.all.daughters.in.comb <- sum(final.sizes[.get.progeny(mother.daughter, comb)])
            
              ## if b is bigger than the selected daughter clones, it may contain mutations present in both the selected clone and the founder population (prob.this.combination)
              
              if(b >= total.size.of.all.daughters.in.comb){ 
                
                prob.this.combination <- probability.this.combination(bin.size = a, clone.size.mother = clone.size.now, n.daughters.present = length(comb), n.daughters.absent = length(daughters.this.clone) - length(comb))
                prob.this.combination.upper <- probability.this.combination(bin.size = c(a[-1], 2*N), clone.size.mother = clone.size.now, n.daughters.present = length(comb), n.daughters.absent = length(daughters.this.clone) - length(comb))
                
                # if the mutation ends up in the daughter clones, drift accounts only for the remaining size (b - total.size.of all.daughters.in.comb)
                res.this.combination <- histogram.drift(lower.bins.1 = a - 1, n.muts = new.muts[clone,], 
                                                        lower.bins.2 = round(b - total.size.of.all.daughters.in.comb), N = N,
                                       bin.p1 = prob.this.combination, bin.p2 = prob.this.combination.upper,
                                       lambda = lambda.ss, delta = delta.this.clone.until.tend, t = time.span.tend)
                
                
              }else{ ## mutations present in at least b cells, where b <= f.sel. --> cumulative distribution from founder cell population +
                
                prob.this.combination <- probability.this.combination(bin.size = a, clone.size.mother = clone.size.now, n.daughters.present = length(comb), n.daughters.absent = length(daughters.this.clone) - length(comb))
                prob.this.combination.upper <- probability.this.combination(bin.size = c(a[-1], 2*N), clone.size.mother = clone.size.now, n.daughters.present = length(comb), n.daughters.absent = length(daughters.this.clone) - length(comb))
                
                # if the mutation ends up in the daughter clones, drift accounts only for the remaining size (b - total.size.of all.daughters.in.comb)
                res.this.combination <- histogram.drift(lower.bins.1 = a - 1, n.muts = new.muts[clone,], lower.bins.2 = 0, N = N,
                                       bin.p1 = prob.this.combination, bin.p2 = prob.this.combination.upper,
                                       lambda = lambda.ss, delta = delta.this.clone.until.tend, t = time.span.tend)
                
                
              }
            
            res.this.combination
          })))
          # expected number of mutations in clone size b if not present in any daughter
          
          prob.no.daughter <- probability.this.combination(bin.size = a, clone.size.mother = clone.size.now, n.daughters.present = 0, n.daughters.absent = length(daughters.this.clone))
          prob.no.daughter.upper <- probability.this.combination(bin.size = c(a[-1], 2*N), clone.size.mother = clone.size.now, n.daughters.present = 0, n.daughters.absent = length(daughters.this.clone))
          
          res.no.daughter <- histogram.drift(lower.bins.1 = a, n.muts = new.muts[clone,], lower.bins.2 = b, N = N,
                                             bin.p1 = prob.no.daughter, bin.p2 = prob.no.daughter.upper,
                                             lambda = lambda.ss, delta = delta.this.clone.until.tend, t = time.span.tend)
          
          res.daughters + res.no.daughter
          
        })
      }else{ # if b is smaller than the selected daughter clones, all mutations must stem from the current mother clone 
        muts.tend[clone,] <- muts.tend[clone,] + sapply(b, function(b){
          
        #  if(b >= total.size.of.all.daughters.this.clone){return(0)}
          res <- histogram.drift(lower.bins.1 = a, n.muts = new.muts[clone,], lower.bins.2 = b, N = N,
                                 bin.p1 = 1, bin.p2 = 1,
                                 lambda = lambda.ss, delta = delta.this.clone.until.tend, t = time.span.tend)
          
          res
          
        })
      }
      
      
    }
    
    # II. new mutations introduced in this time.interval
    
    if(upper.t == t.end){
      new.muts <- matrix(0, nrow=length(t.s), ncol=length(b)) # define a earlier, fine-grained histogram for computation
    }
    
    for(clone in clones){
      
      if(t.s[clone] > lower.t){next}
      
      # I. fate of mutations already present in the system
      clone.size.now <- cell.states[which(cell.states[,1] >= lower.t)[1], 1 + clone]
      clone.size.upper <- cell.states[which(cell.states[,1] >= upper.t)[1], 1 + clone]
      clone.size.tend <- final.sizes[clone]
      # if the clone expands, set its delta according to the selective advantage; otherwise approximate by linear b-d process
      if(clone.size.upper > clone.size.now){
        delta.this.clone.until.upper <- lambda.ss * s[clone]
        # our computation does an exponential growth approximation, but we model the population overall with a competition model. Hence, if the clone is already in the competition phase, we would overshoot if not correcting the time span that accounts for exponential expansion only.
        time.span.upper <- log(clone.size.upper/clone.size.now)/(lambda.ss - delta.this.clone.until.upper) 
      }else if(clone.size.upper == clone.size.now){
        delta.this.clone.until.upper <- lambda.ss
        time.span.upper <- upper.t - lower.t
      }else{
        delta.this.clone.until.upper <- .approximate.delta(lambda.ss, clone.size.now, upper.t - lower.t, cell.states[which(upper.t <= cell.states[,1])[1],1 + length(s) + clone] -
                                                             cell.states[which(cell.states[,1] >= lower.t)[1], 1 + length(s) + clone])
        time.span.upper <- upper.t - lower.t
      }
      if(clone.size.tend > clone.size.now){
        delta.this.clone.until.tend <- lambda.ss * s[clone]
        # our computation does an exponential growth approximation, but we model the population overall with a competition model. Hence, if the clone is already in the competition phase, we would overshoot if not correcting the time span that accounts for exponential expansion only.
        time.span.tend <- log(clone.size.tend/clone.size.now)/(lambda.ss - delta.this.clone.until.tend) 
      }else if(clone.size.tend == clone.size.now){
        delta.this.clone.until.tend <- lambda.ss
        time.span.tend <- t.end - lower.t
      }else{
        delta.this.clone.until.tend <- .approximate.delta(lambda.ss, clone.size.now, t.end - lower.t, cell.states[nrow(cell.states),1 + length(s) + clone] -
                                                            cell.states[which(cell.states[,1] >= lower.t)[1], 1 + length(s) + clone])
        time.span.tend <- t.end - lower.t
      }
      if(round(lambda.ss, digits=10) != round(delta.this.clone.until.upper, digits=10)){
        if(upper.t == t.end){
          new.muts[clone,] <- sapply(b, function(b){
            mutations.noncritical.bd(lambda = lambda.ss, delta = delta.this.clone.until.upper, n.min = b, t.end = time.span.tend, mu = mu, N0 = clone.size.now,
                                     N = final.sizes[clone])
          })
          new.muts[clone,b>=final.sizes[clone]] <- 0
        }else{
          new.muts[clone,] <- sapply(a, function(b){
            mutations.noncritical.bd(lambda = lambda.ss, delta = delta.this.clone.until.upper, n.min = b, t.end = time.span.upper, mu = mu, N0 = clone.size.now,
                                     N = clone.size.upper)
          })
          new.muts[clone,a>=clone.size.upper] <- 0
        }
        
        new.muts[clone,][new.muts[clone,]<0] <- 0 # for very small expansions, our approximation can yield negative results; take them out
      }else{
        time.span.upper <- upper.t - lower.t
        if(upper.t == t.end){
          new.muts[clone,] <- sapply(b, function(b){
            if(b > 2*clone.size.now){return(0)}
            mutations.during.steady.state(lambda = lambda.ss, n.min = b, t.end = time.span.tend, N = clone.size.now, mu = mu)
          })
          new.muts[clone,b>=final.sizes[clone]] <- 0
        }else{
          new.muts[clone,] <- sapply(a, function(b){
            if(b > 2*clone.size.now){return(0)}
            mutations.during.steady.state(lambda = lambda.ss, n.min = b, t.end = time.span.upper, N = clone.size.now, mu = mu)
          })
          new.muts[clone,a>=clone.size.upper] <- 0
        }
        
        new.muts[clone,][new.muts[clone,]<0] <- 0 # for very small expansions, our approximation can yield negative results; take them out
      }
      
    }
    
  }
  
  muts.tend <- colSums(muts.tend + new.muts)
  muts.tend
}



#' Approximate deltas to approximate non-linear decline with a linear b-d-process parametrized by the number of death events
#' 
#' @param N the clone size at the start of contraction
#' @param lambda the division rate
#' @param t the time span
#' @param D the number of death events
#' @return the death rate yielding the same number of death events if modeling exponential decay

.approximate.delta <- function(lambda, N, t, D){
  ## function to determine the loss rate in case of exponential decay when requiring a fixed number of death events D.
  fun <- function(lambda, delta, N, t, D){
    delta*N/(lambda-delta)*(exp((lambda-delta)*t) - 1) - D
  }
  ## re-scale lambda to 1 to make things easier
  upper <- fun(lambda=1, N=N, delta=10, t=t*lambda, D=D)
  lower <- fun(lambda=1, N=N, delta=1+10^-10, t=t*lambda, D=D)
  if((upper>0 & lower>0) | (upper < 0 & lower < 0)){
    return(lambda)
  }
  ## approximate death rate
  delta.founder <- uniroot(fun, interval=c((1+10^-10), 10), lambda=1, N=N, t=t*lambda, D=D)
  ## scale back
  delta.founder <- delta.founder$root*lambda
  delta.founder
}


#' Forward dynamics of j clones
#' 
#' @param N the compartment size
#' @param init.cells vector with the initial condition of the system (number of cells per clone)
#' @param lambda.ss cell division rate during homeostasis
#' @param delta.ss differentiation rate during homeostasis
#' @param lambda.exp division rate during initial expansion
#' @param delta.exp loss rate during initial expansion
#' @param s vector with selective advantages associated with the j-th driver. Selection is modeled as a reduction of the differentiation rate, so 0 <= s <= 1
#' @param t.s vector of length j with the time points at which the selected advantages were introduced
#' @param mother.daughter mother-daughter relationships. Matrix with 2 columns encoding mother and daughter for each pair
#' @param t the time point of evaluation
#' @param resolution the time resolution of the simulation
#' @return A data.frame reporting the system state in the following order: time, cell count for each clone
.forward_dynamics <- function(N, init.cells, lambda.ss, delta.ss, lambda.exp, delta.exp, s, t.s, mother.daughter, t, resolution = 0.01){
  
  colnames(mother.daughter) <- c("M", "D")
  
  if(any(s>1)){
    stop("Values of ", s, " must be <= 1.")
  }
  
  output <-c() # collect states here
  
  cell.states <- init.cells
  death.states <- rep(0, length(s)) # collect the number of death events per clone
  live.index <- which(init.cells!=0) # model only clones that actually exist; better as it prevents numerical instability
  zero.index <- which(init.cells==0)
  # bundle the parameters in a list
  parms <- list(N=N, lambda=lambda.exp, delta=delta.exp, s=s[live.index])
  
  # if we start at t=0, model initial expansion
  if(t[1]==0){
    pre.init <- c(cell.states[live.index], death.states[live.index])
    
    # initial expansion:
    duration <- log(N)/(lambda.exp[1] - delta.exp[1])
    state <- deSolve::ode(as.vector(pre.init), seq(0, duration, resolution), .clonal_dynamics_odes, parms = parms)
    state[,"time"] <- state[,"time"] - duration
    # get cell state and death event state trajectories
    cell.state.traj <- state[,1+(1:length(s[live.index]))]
    death.state.traj <- state[,(length(s[live.index]) + 2): (2*length(s[live.index])+1)]
    # append zeros for non-existing clones and order correctly
    cell.state.traj <- cbind(cell.state.traj, matrix(0, nrow(state), length(zero.index)))
    cell.state.traj <- cell.state.traj[,order(c(live.index, zero.index))]
    death.state.traj <- cbind(death.state.traj, matrix(0, nrow(state), length(zero.index)))
    death.state.traj <- death.state.traj[,order(c(live.index, zero.index))]
    state <- cbind(state[,1], cell.state.traj, death.state.traj)
       
    output <- rbind(output, state)
    cell.states <- state[nrow(state),1+(1:length(s))]
    death.states <- state[nrow(state),(length(s) + 2): (2*length(s)+1)]
  }
  # homeostasis:
  init <- c(cell.states[live.index], death.states[live.index])
  parms <- list(N=N, lambda=lambda.ss, delta=delta.ss, s=s[live.index])
  # only consider clones born after the start of the simulation
  t.s <- t.s[t.s >= min(t),drop=F]
  
  for(k in 1:(length(t.s)+1)){
    
    end.interval <- c(t.s, max(t))[k]
    start.interval <- c(t[1], t.s)[k]
    time.points.this.interval <- sort(unique(c(t[t >= start.interval & t <= end.interval], t.s[t.s >= start.interval & t.s <= end.interval])))
    if(length(time.points.this.interval)<=1){next}
    state <- deSolve::ode(as.vector(init), time.points.this.interval, .clonal_dynamics_odes, parms = parms)
    # get cell state and death event state trajectories
    cell.state.traj <- state[,1+(1:length(s[live.index]))]
    death.state.traj <- state[,(length(s[live.index]) + 2): (2*length(s[live.index])+1)]
    # append zeros for non-existing clones and order correctly
    cell.state.traj <- cbind(cell.state.traj, matrix(0, nrow(state), length(zero.index)))
    cell.state.traj <- cell.state.traj[,order(c(live.index, zero.index))]
    cell.state.traj <- replace(cell.state.traj, cell.state.traj < 0, 0)
    death.state.traj <- cbind(death.state.traj, matrix(death.states[zero.index], nrow(state), length(zero.index)))
    death.state.traj <- death.state.traj[,order(c(live.index, zero.index))]
    state <- cbind(state[,1], cell.state.traj, death.state.traj)
    
    output <- rbind(output, state)
    cell.states <- state[nrow(state),1+(1:length(s))]
    death.states <- state[nrow(state),(length(s)+2):(2*length(s)+1)]
    
    if(c(t.s, max(t))[k]==max(t)){
      break
    }
    # introduce the new clone
    output <- output[-nrow(output),]
    mdc <- mother.daughter[mother.daughter[,"D"]==(k),]
    cell.states[mdc[2]] <- 1
    cell.states[mdc[1]] <- cell.states[mdc[1]] -1
    if(cell.states[mdc[1]] < 0){
      cell.states[mdc[1]] <- 0
    }
    live.index <- which(cell.states!=0) # model only clones that actually exist; better as it prevents numerical instability
    zero.index <- which(cell.states==0)
    # bundle the parameters in a list
    parms <- list(N=N, lambda=lambda.ss, delta=delta.ss, s=s[live.index])
    # initial condition for next interval
    init <- c(cell.states[live.index], death.states[live.index])
  }
  if(any(duplicated(output[,1]))){
    output <- output[-duplicated(output[,1]),]
  }
  
  return(output)
}


#' Dynamics of normal cells and j selected clones 
#' 
#' @param N the carrying capacity
#' @param init vector of length j with the initial condition of the system (number of cells per clone)
#' @param lambda cell division rate
#' @param delta differentiation rate
#' @param s vector with selective advantage associated with the j-th driver. Selection is modeled as a reduction of the differentiation rate, so 0 <= s <= 1
#' @param t the time point of evaluation
#' @return The system state at time t
#' @export
.clonal_dynamics <- function(N, init, lambda, delta, s, t){
  parms <- list(N=N, lambda=lambda, delta=delta, s=s)
  ## ode requires state as vector and parameters as list
  state <- ode(as.vector(init), c(0,t), .clonal_dynamics_odes, parms = parms)
  state <- state[2,-1]
  return(state)
}

.clonal_dynamics_odes <- function(time, state, parms){
  with(as.list(c(state, parms)), {
    loss.rates <- .loss_rates(delta, s, N, state[1:length(s)])
    dn <- lambda*state[1:length(s)] - loss.rates*state[1:length(s)]
    dndeath <- loss.rates*state[1:length(s)]
    return(list(c(dn, dndeath)))
  })
}


#' Compute the final size of the clones from the input parameters
#' 
#' @param N the compartment size
#' @param lambda.ss cell division rate during homeostasis
#' @param delta.ss differentiation rate during homeostasis
#' @param lambda.exp division rate during initial expansion
#' @param delta.exp loss rate during initial expansion
#' @param size vector of length j with the input size of the selcted clones. Using exponential growth approximation, the selective advantage, `s` will be computed from `size`. The actual clone sizes will then be computed using the values of `s` and a clonal competition model.
#' @param t.s vector of length j with the time points at which the selected advantages were introduced
#' @param mother.daughter mother-daughter relationships. Matrix with 2 columns encoding mother and daughter for each pair
#' @param t.end the time point of evaluation
#' @return A vector reporting the system state at the final time point for each clone


.compute_actual_size <- function(t.s, mother.daughter, N, lambda.ss, delta.ss, lambda.exp, 
                                 delta.exp, size, t.end){
  
  s <- c(1, 1 - log(size*N)/(lambda.ss*(t.end - t.s[-1])))
  s[s<0] <- 0
  s[s>1] <- 1
  
  colnames(mother.daughter) <- c("M", "D")
  clone.order <- data.frame(new.id = 1:length(t.s), old.id = order(t.s))
  
  mother.daughter <- rbind(c(1, 1), mother.daughter)
  mother.daughter <- apply(mother.daughter, 2, function(x) {
    sapply(x, function(y) {
      clone.order[clone.order$old.id == y, "new.id"]
    })
  })
  mother.daughter <- mother.daughter[order(mother.daughter[, 
                                                           "D"]), , drop = F]
  s <- s[clone.order$new.id]
  t.s <- sort(t.s)
  
  cell.states <- SCIFER:::.forward_dynamics(N = N, init.cells = c(1, 
                                                                  rep(0, length(s) - 1)), lambda.ss = lambda.ss, delta.ss = lambda.ss, 
                                            lambda.exp = lambda.exp, delta.exp = delta.exp, s = s, 
                                            t.s = t.s, mother.daughter = mother.daughter, t = seq(0, 
                                                                                                  t.end, length.out = 1000))
  final.sizes <- cell.states[nrow(cell.states), 1 + (1:length(s))]
  to.remove <- c()
  for (clone in mother.daughter[, "D"]) {
    daughters.this.clone <- setdiff(.get.progeny(mother.daughter, 
                                                 clone), clone)
    total.size <- sum(cell.states[nrow(cell.states), c(clone, 
                                                       daughters.this.clone) + 1])
    if (total.size/N < min.clone.size) {
      to.remove <- unique(c(to.remove, clone, daughters.this.clone))
    }
  }
  to.keep <- setdiff(1:length(s), to.remove)
  if (length(to.remove) > 0) {
    s <- s[-to.remove]
    t.s <- t.s[-to.remove]
    final.sizes <- final.sizes[-to.remove]
    cell.states <- cell.states[, -c(to.remove + 1, (length(s) + 
                                                      2 + to.remove)), drop = F]
    mother.daughter <- mother.daughter[-which(mother.daughter[, 
                                                              "M"] %in% to.remove | mother.daughter[, "D"] %in% 
                                                to.remove), , drop = F]
  }
  if (nrow(mother.daughter) == 1) {
    return(rep(0, length(t.s) -1))
  }
  else {
    id.conversion <- data.frame(new.id = 1:length(t.s), old.id = mother.daughter[, 
                                                                                 "D"])
    mother.daughter <- apply(mother.daughter, 2, function(x) {
      sapply(x, function(y) {
        id.conversion[id.conversion$old.id == y, "new.id"]
      })
    })
    cell.states <- SCIFER:::.forward_dynamics(N = N, init.cells = c(1, 
                                                                    rep(0, length(s) - 1)), lambda.ss = lambda.ss, delta.ss = lambda.ss, 
                                              lambda.exp = lambda.exp, delta.exp = delta.exp, s = s, 
                                              t.s = t.s, mother.daughter = mother.daughter, t = seq(0, 
                                                                                                    t.end, length.out = 1000))
    final.sizes <- cell.states[nrow(cell.states), 1 + (1:length(s))]
  }
  final.sizes.all.clones <- rep(NA, length(s))
  final.sizes.all.clones[to.keep] <- final.sizes
  final.sizes.all.clones[to.remove] <- 0
  return(final.sizes.all.clones)
}



#' Determine loss rates from logistic growth
#' @description
#' We assume that different clones compete with each other and that competition acts on the loss rate. This function computes the net loss rates for a given state in dependence of the selective advantages, according to
#' \enumerate{
#' \item rho_D0D0=rho_D1D1=rho_D2D2=...=rho_DiDi=delta/K
#' \item rho_D1D0=r1/K, rhoD2D0=r2/K, rhoD2D1=r2/r1/K, ...
#' \item rho_D0D1=2-r1, rho_D0D2=2-r2, rhoD1D2=2-r2/r1, ...
#' },
#' where `D`stands for `D`river
#' @param delta a vector of length i containing the loss rate per compartment
#' @param s a vector of selective advantages associated with clone j
#' @param N the compartment size
#' @param cell.state vector of length j with the number of cells of the j-th clone
#' @return an vector of length j with the current loss rates
#' @keywords internal
#' @export
#' 

.loss_rates <- function(delta, s, N, cell.state){
  if(any(s>1)){
    stop("Values of ", s, " must be <= 1.")
  }
  
  ## competition matrix; can be computed when requiring steady cell.state a + b + c=N
  rho <- matrix(delta/N, length(s), length(s))
  for(k in 1:nrow(rho)){
    for(l in 1:ncol(rho)){
      if(k==l){
        next
      }else if(k > l){
        rho[k,l] <- rho[k,l]*s[k]/s[l]
      }else{
        rho[k,l] <- rho[k,l]*(2 - s[l]/s[k])
      }
    }
  }
  if(N==0){
    rho <- matrix(0, length(s), length(s))
  }
  res <- colSums(t(rho)*cell.state)
  
  return(res)
}


#' Get all subclonal descendants from a given clone
#' 
#' @param mother.daughter mother-daughter relationships. Matrix with 3 columns encoding mother, daughter and birth-compartment for each pair
#' @param id IDs of the clone of interest
#' @return A vector with the progeny IDs.

.get.progeny <- function(mother.daughter, id){
  res <- id
  while(length(id)>0){
    daughter.ids <- mother.daughter[mother.daughter[,"M"] %in% id,"D"]
    id <- setdiff(daughter.ids, res)
    res <- c(res, daughter.ids)
  }
  return(unique(res))
}


#' Clonal combinations
#' 
#' @param clone.ids IDs of daughter clones from the same mother
#' @return All clonal combinations in which a mutation acquired in the mother can end up.
#' @export
.clonal.combinations <- function(clone.ids){
  n.clones <- length(clone.ids)
  # get all combinations of the mutation being in any subset of the clones
  combinations <- sapply(1:n.clones, function(x){apply(combn(n.clones,x), 2, list)})
  combinations <- unlist(combinations, recursive = F)
  combinations <- unlist(combinations, recursive = F)
  # rename using clone ids
  combinations <- lapply(combinations, function(x){
    clone.ids[x]
  })
  return(combinations)
}


#########################################################################################################################################
####### Simulate whole genome sequencing
#########################################################################################################################################

#' Simulate read sampling in WGS by Binomial sampling
#' 
#' @param clone.sizes vector of clone sizes at which cumulative mutation counts were measured
#' @param expected.mutations expected number of mutations at each clone size
#' @param depth sequencing depth
#' @param sensitivity logical, if sensitivity of sequencing method should be taken into account in addition to binomial noise. Requires a specification for `false.negative.per.vaf`.
#' @param false.negative.per.vaf optional, a matrix with columns corresponding to the measured VAFs and rows corresponding to individual measurements of the false negative rate at this VAF in addition to binomial noise. Must be provided if `sensitivity=T`
#' @return A vector of simulated VAFs

.simulated.wgs.data <- function(clone.sizes, expected.mutations, depth=90, sensitivity=T, false.negative.per.vaf, min.vaf=0.05){
  
  simulated.vafs <- c()
  sim <- apply(rbind(clone.sizes, expected.mutations), 2, function(x){
    if(x[2] < 1){
      x[2][x[2]<0] <- 0
      x[2] <- sample(size = 1, c(0,1), prob = c(1-x[2], x[2]))
    }
    simulated.depth <- rpois(n = round(x[2]), lambda = depth)
    res <- rbinom(n=round(x[2]), size = simulated.depth, prob = x[1]/2)
    res[res < 3] <- 0
    res <- res/simulated.depth
  })
  sim <- unlist(sim)
  
  sim <- sim[sim!=0]
  
  ## simulate accuracy based on blood-bone-marrow-comparison
  if(sensitivity){
    false.negative.per.vaf <- false.negative.per.vaf[,seq(0.05, 1, 0.01)>=min.vaf]
    ## additional unexplained error
    false.negative.per.vaf. <- t(t(false.negative.per.vaf) - sapply(seq(0.05, 1, 0.01), function(x){pbinom(0.05*depth, size=depth, prob=x)} ))
    false.negative.per.vaf.[false.negative.per.vaf. < 0 ] <- 0
    sim <- sapply(sim, function(x){
      ## round to 2 digits
      y <- round(x, digits=2)
      which.fn <- which.min((vafs.of.interest - y)^2)[1]
      ## sample a false negative rate from mean and variance at this vaf
      fn.rate <- rnorm(1, mean= mean(false.negative.per.vaf.[,which.fn]),
                       sd = sd(false.negative.per.vaf.[,which.fn]))
      
      
      
      if(fn.rate < 0){
        fn.rate <- 0
      }else if(fn.rate>1){
        fn.rate <- 1
      }
      
      detected <- sample(x=c(1,0), 1, prob = c(1-fn.rate, fn.rate))
      if(detected==1){
        x
      }else{
        0
      }
    })
  }
  sim
}


#' Simulate sampling as in single-cell sequencing by simulating the actual sampling of cells (i.e. it's not the expected value but a sampling instance)
#' 
#' @param clone.sizes vector of clone sizes at which cumulative mutation counts were measured
#' @param expected.mutations expected number of mutations at each clone size
#' @param ncells sequenced cells
#' @return A vector of simulated VAFs

.simulated.scWGS.data <- function(clone.sizes, expected.mutations, ncells=100, min.vaf=0.05){
  
  simulated.vafs <- c()
  sim <- apply(rbind(clone.sizes, expected.mutations), 2, function(x){
    if(x[2] < 1){
      x[2][x[2]<0] <- 0
      x[2] <- sample(size = 1, c(0,1), prob = c(1-x[2], x[2]))
    }
    res <- rbinom(n=round(x[2]), size = ncells, prob = x[1])/(2*ncells) ## transform into VAF
    res 
  })
  sim <- unlist(sim)
  
  sim <- sim[sim!=0]
  
  sim
  
}

#' Wrapper function to simulate sequencing either by bulk WGS or by scWGS
#' 
#' @param seqtype string, specifying the sequencing method. Must be either "bulk" or "sc"
#' @param clone.sizes vector of clone sizes at which cumulative mutation counts were measured
#' @param expected.mutations expected number of mutations at each clone size
#' @param depth sequencing depth, only specify for bulk WGS
#' @param false.negative.per.vaf optional, a matrix with columns corresponding to the measured VAFs and rows corresponding to individual measurements of the false negative rate at this VAF in addition to binomial noise.
#' @param sensitivity logical, if sensitivity of sequencing method should be taken into account in addition to binomial noise. Requires a specification for `false.negative.per.vaf`.
#' @param min.vaf the minimal VAF to return
#' @param ncells the number of sequenced cells, only specify if `seqtype=="sc"`.
#' @return A vector of simulated VAFs
#' @export

simulated.data <- function(seqtype, clone.sizes, expected.mutations, depth=90, ncells=100, sensitivity=T, false.negative.per.vaf, min.vaf=0.05){
  
  if(seqtype=="bulk"){
    sim <- .simulated.wgs.data(clone.sizes, expected.mutations, depth, sensitivity, false.negative.per.vaf, min.vaf)
  }else if(seqtype=="sc"){
    sim <- .simulated.scWGS.data(clone.sizes, expected.mutations, ncells, min.vaf)
  }
  sim
  
}
