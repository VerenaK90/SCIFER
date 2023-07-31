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
        total <- .approx.count(t, mu, lambda, delta, n.min=11, n.max=100*max(N,N0), t.end) + .exact.count(t, mu, lambda, delta, 10, t.end)
        if(n.min > 1){
          res <- total - .exact.count(t, mu, lambda, delta, n.min-1, t.end)
        }else{
          res <- total 
        }
        return(res)
      }else{
        .approx.count(t, mu, lambda, delta, n.min=1, n.max=Inf, t.end) - .approx.count(t, mu, lambda, delta, n.min=1, n.max=n.min, t.end)
      }
      
    })
    return(res)
    
  }
  ## Compute the number of mutations that are present in at least 1 cell and at most n.min cells
  ## The sum necessary in order to compute the cumulative distribution, is here replaced by integration.
  sapply(n.min, function(n.min){
    
    res <- .approx.count(t, mu, lambda, delta, n.min=1, n.max=100*max(N0, N), t.end) - .approx.count(t, mu, lambda, delta, n.min=1, n.max=n.min, t.end) 
    return(res)
  })
  
}

## Exact number of mutations present in at least n cells at time t
.exact.count <- function(t, mu, lambda, delta, n, t.end){
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

.approx.count <- function(t, mu, lambda, delta, n.min, n.max, t.end){
  integrand <- function(t, mu, lambda, delta, n.min, n.max){
    mu*lambda*N0*exp((lambda - delta)*t)/log(.beta(lambda, delta, t.end-t))*(density.a.b.exact(lambda, delta, t.end-t, 1, n.max) - 
                                                                               density.a.b.exact(lambda, delta, t.end-t, 1, n.min))
    
  }
  total <- integrate(integrand, lower=0, upper=t.end, mu=mu, lambda=lambda, delta=delta, n.min=n.min, n.max=n.max)$value
  total
}

#########################################################################################################################################
####### Critical b-d-process 
#########################################################################################################################################

#' Clone size distribution of a critical birth-death process (exact).
#' 
#' @description Exact solution to grow from a clone of size a to a clone of size b within t in a critical birth-death process.
#' @param lambda proliferation rate 
#' @param a clone size at `t=0`
#' @param b clone size at `t=t`
#' @param t time
#' @references  Bailey, NTJ (1964). The elements of stochastic processes with applications to the natural sciences, Wiley (New York).
#' @return The probability to grow from `a` to `b` within `t`.
#' @export

p.ss.exact <- function(lambda, a, b, t){
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
  res[a==0] <- 0
  
  res
}

#' Clone size distribution of a critical birth-death process (approximate).
#' 
#' @description Probabilitly to grow from a clone of size `a` to a clone of size `b` within `t` in a critical birth-death process using a parametrized \eqn{\Gamma}-distribution for approximation.
#' @param lambda proliferation rate 
#' @param a clone size at `t=0`
#' @param b clone size at `t=t`
#' @param t time
#' @param mode either "density" if density distribution is to be returned or "cumulative", defaults to "density".
#' @references  Bailey, NTJ (1964). The elements of stochastic processes with applications to the natural sciences, Wiley (New York).
#' The function is approximated with a \eqn{\Gamma}-distribution that is parametrized with 
#' \eqn{shape = \mu^2/\sigma, scale = \sigma/\mu}, where
#' \eqn{\mu = a, \sigma = 2a \lambda t}
#' @return The approximate probability to grow from `a` to `b` within `t`.
#' @export

p.ss.approx <- function(lambda, a, b, t, mode="density"){
  
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
####### Supercritical b-d-process followed by critical b-d-process
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
  ## how much more mutations lie on the left as compared to the right border of the bin?
  bin.factor <- m.tss/c(m.tss[-1],0)
  
  ## no homeostatic phase if measuring directly after exponential expansion
  if(t.end==0){
    return(mutations.noncritical.bd(lambda = lambda.exp, delta = delta.exp, t.end = t.ss, mu = mu, n.min = b, N=N))
  }
  
  ## compute for each interval the expected number of mutations that will have grown larger or equal than b cells at t.end
  ## do this by integration if a is large and by sums else
  
  
  mutations.from.expansion <- sapply(b, function(b){
    
    ## Probability that a mutation present in bin a at t.ss will drift to a clone size larger than b within t.end
    ## Takes the weighted average between left and right side of the bin. Weighted according to the number of mutations present in the bin and the next bin (bin.factor, see above)
    mean.prob.b <- c()
    for(j in 1:length(a)){
      
      bin.factor[is.infinite(bin.factor)] <- 10^8
      bin.factor[is.na(bin.factor)] <- 1
      
      ## the clone size at the upper border of the bin. Take 2*N for the largest bin. 
      
      upper.a <- round(c(a[-1],2*N)[j])
      ## from critical b-d-process
      p <- lambda.ss*t.end/(1+lambda.ss*t.end)
      ## mean and variance from critical b-d- process to parametrize gamme distribution
      ## mean = a, variance = 2a lambda*t
      ## translates to shape = mean^2/var, scale=var/mean
      
      #if(b>100){## use gamma distribution if b >20
      if(b>20){
        ## lower border of bin
        mean=a[j]
        var=2*a[j]*lambda.ss*t.end
        ## upper border of the bin
        upper.mean=upper.a
        upper.var=2*upper.a*lambda.ss*t.end
        ## weighted average of cumulative distribution
        mean.prob.b[j] <- sum(c(bin.factor[j]*(1-  pgamma(b-1, shape=mean^2/var, scale=var/mean)),
                                1 - pgamma(b-1, shape=upper.mean^2/upper.var, scale=upper.var/upper.mean)))/(bin.factor[j]+1)
      }else{
        ## weighted average of cumulative distribution
        mean.prob.b[j] <- sum(c(bin.factor[j]*(1-sum(sapply(0:(b-1), function(b){p.ss.exact(lambda.ss, a[j], b, t.end)}))),
                                (1-sum(sapply(0:(b-1), function(b){p.ss.exact(lambda.ss, upper.a, b, t.end)})))))/(bin.factor[j]+1)
        
      }
    }
    
    mutations.from.expansion <- sum(m.tss*mean.prob.b)
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
        p.mut.in.sel*mu*lambda*exp((lambda - delta)*t)*(density.a.b.exact(lambda, delta, t.end-t, 1, 0) )
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
  
  ode.system.competition <- function(t, y, parms){
    with(as.list(c(parms, y)), {
      ## number of selected cells
      dm <- lambda*m*(1-m/N-s*(N-m)/N)
      ## number of death events in founder cell population
      dDeath <- lambda*((N-m)/N + (2-s)*m/N)*(N-m)
      ## number of divisions in the expanding population
      dDiv <- lambda*m
      list(c(dm, dDeath, dDiv))
    })
  }
  
  number.of.events <- deSolve::ode(func=ode.system.competition, times=c(0,t.end-(t.s)), y=c(m=1, Death=0, Div=0), parms=c(lambda=lambda.ss, s=s, N=N))
  number.of.deaths.in.founder <- number.of.events[2,"Death"]
  ## function to determine the loss rate in case of exponential decay when requiring a fixed number of death events D.
  fun <- function(lambda, delta, N, t, D){
    delta*N/(lambda-delta)*(exp((lambda-delta)*t) - 1) - D
  }
  ## re-scale lambda to 1 to make things easier
  upper <- fun(lambda=1, N=N, delta=10, t=(t.end-t.s)*lambda.ss, D=number.of.deaths.in.founder)
  lower <- fun(lambda=1, N=N, delta=1+10^-10, t=(t.end-t.s)*lambda.ss, D=number.of.deaths.in.founder)
  if((upper>0 & lower>0) | (upper < 0 & lower < 0)){
    return(1000)
  }
  ## approximate death rate
  delta.founder <- uniroot(fun, interval=c((1+10^-10), 10), lambda=1, N=N, t=(t.end-(t.s ))*lambda.ss, D=number.of.deaths.in.founder)
  ## scale back
  delta.founder <- delta.founder$root*lambda.ss
  
  
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
  
  bin.factor <- mutations.at.t.s/c(mutations.at.t.s[-1],0)
  bin.factor[is.infinite(bin.factor)] <- 10^8
  bin.factor[is.na(bin.factor)] <- 1
  
  ## the clone size at the upper border of the bin
  upper.a <- c(a[-1], 2*N)
  
  ## take the average from the lower and upper border
  mutations.at.t.end <- sapply(b, function(b){
    ## if b is bigger than the selected clone, it may contain mutations present in the founder cell population only (probability 1-a/N)
    ## or mutations present in both the selected clone and the founder population (probability a/N)
    
    if(b >= f.sel*N){ ## weighted average of the probability to drift from the lower and upper boundary of the interval to a size of at least b
      weighted.average <- apply(rbind(a, upper.a, bin.factor), 2, function(x){
        ## probability that the mutation is present in the selected clone computed at both borders of the bin
        prob.sel <- x[1]/N
        prob.sel.upper <- x[2]/N
        (x[3]* (prob.sel*(1-p.a.b(a=x[1]-1, b=round(b-f.sel*N), lambda=lambda.ss, delta=delta.founder, t=t.end-t.s)) + ## mutation in both populations (-1 because 1 cell is the selected clone)
                  (1-prob.sel)*(1-p.a.b(a=x[1], b=b, lambda=lambda.ss, delta=delta.founder, t=t.end-t.s))) +## mutation only in founder cell 
            
            prob.sel.upper*(1-p.a.b(a=x[2]-1, b=round(b-f.sel*N), lambda=lambda.ss, delta=delta.founder, t=t.end-t.s)) + ## mutation in both populations
            (1-prob.sel.upper)*(1-p.a.b(a=x[2], b=b, lambda=lambda.ss, delta=delta.founder, t=t.end-t.s)) ## mutation only in founder cell 
          
        )/(x[3] + 1)
        
      })
      
      res <- sum(mutations.at.t.s*weighted.average)
    }else{ ## mutations present in at least b cells, where b <= f.sel. --> cumulative distribution from founder cell population +
      ## cumulative distribution from selected clone, but in this case start with f.sel at the lower end as this smaller sizes are not compatible with the selected clone.
      weighted.average <- apply(rbind(a, upper.a, bin.factor), 2, function(x){
        
        prob.sel <- x[1]/N
        prob.sel.upper <- x[2]/N
        
        
        p <- (x[3]* (prob.sel*(1-p.a.b(a=x[1]-1, b=0, lambda=lambda.ss, delta=delta.founder, t=t.end-t.s)) + ## mutation in both populations (-1 because 1 cel is the selected clone)
                       (1-prob.sel)*(1-p.a.b(a=x[1], b=b, lambda=lambda.ss, delta=delta.founder, t=t.end-t.s))) +## mutation only in founder cell 
                
                prob.sel.upper*(1-p.a.b(a=x[2]-1, b=0, lambda=lambda.ss, delta=delta.founder, t=t.end-t.s)) + ## mutation in both populations
                (1-prob.sel.upper)*(1-p.a.b(a=x[2], b=b, lambda=lambda.ss, delta=delta.founder, t=t.end-t.s)) ## mutation only in founder cell 
              
        )/(x[3] + 1)
        
        return(p)
        
      })
      
      res <- sum(mutations.at.t.s*weighted.average)
    }
    res
  })
  
  ## Finally, we need to add up new mutations acquired during the expansion of the selected clone
  
  mutations.from.selected.clone <- mutations.noncritical.bd(lambda.ss, lambda.ss*s, t.end-t.s, mu, b, N=N)
  
  mutations.from.founder.clone <- mutations.noncritical.bd(lambda.ss, delta.founder, t.end - t.s, mu, b, N0=N, N=(1-f.sel)*N)
  ## truncate at b= f.sel*N
  mutations.from.selected.clone[b>f.sel*N] <- 0
  mutations.at.t.end <- mutations.at.t.end + mutations.from.selected.clone + mutations.from.founder.clone
  
  mutations.at.t.end
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
