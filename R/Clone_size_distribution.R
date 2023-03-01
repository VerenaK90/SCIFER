#########################################################################################################################################
####### Non-critical b-d-process 
#########################################################################################################################################

#' Exact probability to grow from a to be within time t according to a supercritical birth-death process. Exact solution according to Bailey, 1964
#'
#' @param lambda proliferation rate 
#' @param delta loss rate
#' @param t time
#' @param a clone size at t=0
#' @param b clone size at t=t
#' @return The probability that a clone of size a grows to size b within t.
#' @examples
#' density.a.b.exact(1, 0, 10, 1, 2)
#' @export

.alpha <- function(lambda, delta, t){
  (delta*exp((lambda - delta)*t) - delta)/
    (lambda*exp((lambda - delta)*t) - delta)
}

.beta <- function(lambda, delta, t){
  (lambda*exp((lambda - delta)*t) - lambda)/
    (lambda*exp((lambda-delta)*t) - delta)
}

density.a.b.exact <- function(lambda, delta, t, a, b){
  if(a==1){
    (1-.alpha(lambda, delta, t))*(1-.beta(lambda, delta, t))*
      .beta(lambda, delta, t)^(b-1)
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

#' Probability of a clone of size a to grow to size b within t according to a supercritical linear birth-death process. 
#'
#' @param lambda proliferation rate 
#' @param delta loss rate
#' @param t time
#' @param a clone size at t=0
#' @param b clone size at t=t
#' @param mode 'density' if density distribution is to be returned , 'cumulative' if cumulative distribution is to be returned. Defaults to 'cumulative'
#' @param approx Approximation to be used. Defaults to 'highnumbers'; i.e. we want to approximate with a gamma distribution if a and b are high.
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


#' Expected number of neutral mutations that are present in at least n.min cells at t.ss in an exponentially growing or contracting tissue.
#'
#' @param lambda proliferation rate 
#' @param delta loss rate
#' @param t.end time
#' @param mu mutation rate per cell division
#' @param n.min minimal clone size at t.ss; can be a value or a vector
#' @param N0 initial population size
#' @param N final population size
#' @param mode if "approx" the sum is approximated by integration. If "exact" the sum is exactly computed for clone sizes between 1 and 10 but beyond that also approximated.
#' @return The expected number of mutations present in at least n.min cells at t.ss in an exponentially growing tissue. The function assumes that mutations are continuously acquired at a constant rate.
#' @export

mutations.noncritical.bd <- function(lambda, delta, t.end, mu, n.min, N0=1, N=N, mode="approx"){
  
  ## Exact number of mutations present in at least n cells at time t
  exact.count <- function(t, mu, lambda, delta, n, t.end){
    sum(sapply(seq(1,n), function(n){
      integrand <- function(t, mu, lambda, delta, n){
        mu*lambda*N0*exp((lambda - delta)*t)*(density.a.b.exact(lambda, delta, t.end-t, 1, n)) 
        
      }
      ## total number of mutations acquired during exponential growth that survived:
      res <- integrate(integrand, lower=0, upper=t.end, mu=mu, lambda=lambda, delta=delta, n=n)$value
      res
    }))
  }
  
  approx.count <- function(t, mu, lambda, delta, n.min, n.max, t.end){
    integrand <- function(t, mu, lambda, delta, n.min, n.max){
      mu*lambda*N0*exp((lambda - delta)*t)/log(.beta(lambda, delta, t.end-t))*(density.a.b.exact(lambda, delta, t.end-t, 1, n.max) - 
                                                                                 density.a.b.exact(lambda, delta, t.end-t, 1, n.min))
      
    }
    total <- integrate(integrand, lower=0, upper=t.end, mu=mu, lambda=lambda, delta=delta, n.min=n.min, n.max=n.max)$value
    total
  }
  
  if(mode=="exact"){
    

    res <- sapply(n.min, function(n.min){
      if(n.min <=10){
        total <- approx.count(t, mu, lambda, delta, n.min=11, n.max=100*max(N,N0), t.end) + exact.count(t, mu, lambda, delta, 10, t.end)
        if(n.min > 1){
          res <- total - exact.count(t, mu, lambda, delta, n.min-1, t.end)
        }else{
          res <- total 
        }
        return(res)
      }else{
        approx.count(t, mu, lambda, delta, n.min=1, n.max=Inf, t.end) - approx.count(t, mu, lambda, delta, n.min=1, n.max=n.min, t.end)
      }
      
    })
    return(res)
    
  }
  ## Compute the number of mutations that are present in at least 1 cell and at most n.min cells
  ## The sum necessary in order to compute the cumulative distribution, is here replaced by integration.
  sapply(n.min, function(n.min){
    
    res <- approx.count(t, mu, lambda, delta, n.min=1, n.max=100*max(N0, N), t.end) - approx.count(t, mu, lambda, delta, n.min=1, n.max=n.min, t.end) 
    return(res)
  })
  
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
  
  ## mean and variance to parametrize gamma distribution
  mean=a
  var=2*a*lambda*t
  
  if(mode=="density"){
    return(dgamma(b, shape=mean^2/var, scale=var/mean))
  }else{
    return(pgamma(b, shape=mean^2/var, scale=var/mean))
  }
}

#' Function to compute the probabilty to grow from a to b in a critical b-d process using automatic switching between exact and approximate solution.
#' Cutoff criterion: a*p*(1-p)>=9 and b*p*(1-p)>=9
#' @param lambda proliferation rate 
#' @param a clone size at t=0
#' @param b clone size at t=t; single value or vector
#' @param t time
#' @param mode either 'density' if density distribution is to be returned or 'cumulative'.
#' @return The probability to grow from a to b within t
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
#' @return The number of mutations that were acquired during steady state and are present in at least n.min cells.
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

#' Mutation accumulation during epxonential expansion followed by homeostasis. 
#' 
#' @param mu mutation rate per cell division
#' @param N population size
#' @param lambda.exp proliferation rate during expansion
#' @param delta.exp loss rate during expansion
#' @param lambda.ss proliferation and loss rate during homeostasis
#' @param t.end end point (starting from homeostasis)
#' @param b minimal clone size of interest. Number or vector. 
#' @param accuracy.a step size in which mutations accumulated during expansion are evaluated between 5 and 100%; defaults to 5%
#' @param return This function returns an approximation by first computing the distribution at the transition time within intervals, then averaging the fate of each interval during homeostasis and adding newly acquired mutations. Returns the number of mutations present in at least b cells.
#' @export


mutational.burden <- function(mu, N, lambda.exp, delta.exp, lambda.ss, t.end, b, accuracy.a = 0.05){

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
  
  mutations.from.steady.state + mutations.from.expansion
}



#########################################################################################################################################
####### Supercritical birth-death process with selection of an advantageous subclone
#########################################################################################################################################

#' Mutation accumulation during exponential expansion followed by homeostasis. 
#' 
#' @param mu mutation rate per cell division
#' @param N population size
#' @param lambda proliferation rate 
#' @param delta loss rate 
#' @param t.end end point 
#' @param t.s time point at which selective advantage is acquired.
#' @param s selective advantage
#' @param b minimal clone size of interest. Number or vector. 
#' @return This function returns an approximation by first computing the distribution at the transition time within intervals, then averaging the fate of each interval during homeostasis and adding newly acquired mutations in a scenario where a subpopulation is under positive selection. Returns the number of mutations present in at least b cells
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
      p.mut.in.sel*mu*lambda*exp((lambda - delta)*t)/log(.beta(lambda, delta, t.end-t))*(density.a.b.exact(lambda, delta, t.end-t, 1, n) - 
                                                                             density.a.b.exact(lambda, delta, t.end-t, 1, 1))
      
    }
    
    ## total number of mutations acquired during exponential growth that survived:
    total <- integrate(integrand, lower=0, upper=t.s, mu=mu, lambda=lambda, delta=delta, n=100*N)$value
    res <- total - integrate(integrand, lower=0, upper=t.s, mu=mu, lambda=lambda, delta=delta, n=max(0,n.min-sel.size))$value
    return(res)
  })
  
  mutations.not.in.selected.clone.prior.t.s <- sapply(b, function(n.min){
    integrand <- function(t, mu, lambda, delta, n){
      p.mut.in.sel <- exp((lambda - delta)*(t.s - t))/exp((lambda - delta)*t.s)
      (1-p.mut.in.sel)*mu*lambda*exp((lambda - delta)*t)/log(.beta(lambda, delta, t.end-t))*(density.a.b.exact(lambda, delta, t.end-t, 1, n) - 
                                                                                           density.a.b.exact(lambda, delta, t.end-t, 1, 1))
      
    }
    
    ## total number of mutations acquired during exponential growth that survived:
    total <- integrate(integrand, lower=0, upper=t.s, mu=mu, lambda=lambda, delta=delta, n=100*N)$value
    res <- total - integrate(integrand, lower=0, upper=t.s, mu=mu, lambda=lambda, delta=delta, n=n.min)$value
    return(res)
  })
  
  mutations.before.t.s <- mutations.in.selected.clone.prior.t.s + mutations.not.in.selected.clone.prior.t.s
  
  
  mutations.from.selected.clone.after.t.s <- sapply(b, function(n.min){
    mutations.after.expansion(lambda = lambda, delta = delta *s, t.ss = t.end - t.s, 
                                                                       mu = mu, n.min = n.min, N = sel.size)})
  mutations.from.founder.clone.after.t.s <- exp((lambda - delta)*t.s)*sapply(b, function(n.min){
    mutations.after.expansion(lambda = lambda, delta = delta, t.ss = t.end - t.s, 
                              mu = mu, n.min = n.min, N = exp((lambda - delta)*(t.end - t.s)), mode="exact")})
  
  
  mutations.before.t.s + mutations.from.selected.clone.after.t.s + mutations.from.founder.clone.after.t.s
  
}

#########################################################################################################################################
####### Suprecritical b-d-process followed by critical b-d-process with selection of an advantageous subclone; subclonal mutation acquired during homeostasis
#########################################################################################################################################

#' Mutation accumulation during exponential expansion followed by homeostasis. 
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
#' @param accuracy.a step size in which mutations accumulated during expansion are evaluated between 5 and 100%; defaults to 5%
#' @param min.clone.size the lower detection limit for selected clones
#' @return This function returns an approximation by first computing the distribution at the transition time within intervals, then averaging the fate of each interval during homeostasis and adding newly acquired mutations in a scenario where a subpopulation is under positive selection. Returns the number of mutations present in at least b cells
#' @export

mutational.burden.with.selection <- function(mu, N, lambda.exp, delta.exp, lambda.ss, t.end, t.s, s, b, accuracy.a = 0.05, min.clone.size = 0.05){
  
  ## initialize mutation count
  mutations.at.t.end <- 0
  
  ## compute the relative size of the selected clone at t.end
  f.sel <- exp((lambda.ss - s*lambda.ss)*(t.end-t.s))/N

  ## in case the selected clone took over, add an offset of the number of mutations in the founder cell to the solution of mutations during steady state (it's again a neutrally expanding clone)
  if(f.sel>1){
    t.ss <- log(N*0.5)/(lambda.ss*(1-s))
    mutations.in.founder.cell <- (log(N)/(lambda.exp - delta.exp)*mu + t.s*lambda.ss*mu)/2
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
  
  number.of.events <- ode(func=ode.system.competition, times=c(0,t.end-(t.s)), y=c(m=1, Death=0, Div=0), parms=c(lambda=lambda.ss, s=s, N=N))
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
####### Suprecritical b-d-process followed by critical b-d-process with selection of two consecutive advantageous subclone
#########################################################################################################################################

#' Mutation accumulation during exponential expansion followed by homeostasis. 
#' 
#' @param mu mutation rate per cell division
#' @param N population size
#' @param lambda.exp proliferation rate during expansion
#' @param delta.exp loss rate during expansion
#' @param lambda.ss proliferation and loss rate during homeostasis
#' @param t.end end point (starting from homeostasis)
#' @param t.s time point at which selective advantage is acquired.
#' @param t.s2 time point of second mutation
#' @param s1 selective advantage associated with first mutation
#' @param s2 selective advantage associated with second mutation
#' @param b minimal clone size of interest. Number or vector. 
#' @return This function returns an approximation by first computing the distribution at the transition time within intervals, then averaging the fate of each interval during homeostasis and adding newly acquired mutations in a scenario with 2 nested clonal selections (clone starts growing at t.s >t.ss and grows with a selective advantage s; clone 2 starts. Returns the number of mutations present in at least b cells
#' @export


mutational.burden.with.nested.selection <- function(mu, N, lambda.exp, delta.exp, lambda.ss, t.end, t.s, t.s2, s1, s2, b){
  

  ## initialize mutation count
  mutations.at.t.end <- 0
  
  ## compute the relative size of the selected clone at t.end
  f.sel <- (exp((lambda.ss - s1*lambda.ss)*(t.end - t.s)) - exp((lambda.ss - s1*lambda.ss)*(t.end - t.s2)) + 
              exp((lambda.ss - s1*s2*lambda.ss)*(t.end - t.s2)))/N

  ## in case the selected clone took over, add an offset of the number of mutations in the founder cell to the solution of mutations during steady state (it's again a pure clone)
  if(f.sel>1){
    ## when did the first clone reach 50% of N?
    func <- function(t.ss){
      if(t.s2 < t.ss){
        (exp((lambda.ss - s1*lambda.ss)*(t.ss - t.s)) - exp((lambda.ss - s1*lambda.ss)*(t.ss - t.s2)) + 
           exp((lambda.ss - s1*s2*lambda.ss)*(t.ss - t.s2)))/N - 0.5
      }else{
        (exp((lambda.ss - s1*lambda.ss)*(t.ss - t.s)))/N - 0.5
      }

    }
    
    t.ss <- uniroot(func, lower = 0, upper = t.end)$root - t.s
    #t.ss <- log(N*0.5)/(lambda.ss*(1-s1))
    t.end.1st.clone <- t.end - t.ss - t.s
    t.s2.rel.to.1st.clone <- t.s2 - t.ss - t.s
    if(t.s2.rel.to.1st.clone < 0){
      t.s2.rel.to.1st.clone <- 0
    }
    mutations.in.founder.cell <- (log(N)/(lambda.exp - delta.exp)*mu + t.s*lambda.ss*mu)/2
    mutations.at.t.end <- mutational.burden.with.selection(mu, N, 1, delta.exp, lambda.ss, t.end.1st.clone, t.s2.rel.to.1st.clone, s2, b) + mutations.in.founder.cell
    return(mutations.at.t.end)
  }else if(f.sel <= 0.05){ ##just take the predicted output at t.end according to homeostatic turnover and neglect expansion of the selected clone
    mutations.at.t.end <- mutational.burden(mu, N, lambda.exp, delta.exp, lambda.ss, t.end, b)
  }
  
  ## else if 5% is reached before, take the distribution at t.s as the input for selection and neglect mutations acquired in the founder cell population after t.s. We thus have two contributions:
  ## selected clone + subclonal mutations acquired during expansion of the latter
  ## founder population:drift of mutations acquired prior to t.s
  
  ## The drift of the founder cell population is non-trivial, due to the non-exponential decay in a competing system.
  ## To assess it, we approximate the loss with exponential decay, requiring the same number of death events in the founder cell population 
  
  ## First, death events while there's only one clone
  ode.system.competition.one.mutant <- function(t, y, parms){
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
  ## count the number of death events until t.s2
  number.of.events <- ode(func=ode.system.competition.one.mutant, times=c(0,t.s2-(t.s)), y=c(m=1, Death=0, Div=0), parms=c(lambda=lambda.ss, s=s1, N=N))

  ## Second, death events with two clones
  ode.system.competition.two.mutants <- function(t, y, parms){
    with(as.list(c(parms, y)), {
      rhom1m2 <- 2-s2
      rhonm1 <- 2-s1
      rhonm2 <- 2-s1*s2
      ## number of selected m1 cells
      dm1 <- lambda*m1*(1-m1/N - s1*(N-m1-m2)/N - rhom1m2*m2/N)
      ## number of selected m2 cells
      dm2 <- lambda*m2*(1-m2/N-s1*s2*(N-m1-m2)/N - s2*m1/N)
      ## number of death events in founder cell population
      dDeath <- lambda*((N-m1-m2)/N + rhonm1*m1/N + rhonm2*m2/N)*(N-m1-m2)

      
      ## number of divisions in the expanding population
      dDiv <- lambda*(m1+m2)
      list(c(dm1, dm2, dDeath, dDiv))
    })
  }
  
  number.of.events <- ode(func=ode.system.competition.two.mutants, times=c(0,t.end-(t.s2)), 
                          y=c(m1 = unname(number.of.events[2,"m"]), m2=1, number.of.events[2,"Death"], number.of.events[2,"Div"]), 
                          parms=c(lambda=lambda.ss, s1=s1, s2=s2, N=N))
  
  number.of.deaths.in.founder <- number.of.events[2,"Death"]
  
  ## solve for delta if modeling with exponential decay
  fun <- function(lambda, delta, N, t, D){
    delta*N/(lambda-delta)*(exp((lambda-delta)*t) - 1) - D
  }
  ## re-scale to 1 to make things easier
  upper <- fun(lambda=1, N=N, delta=10, t=(t.end-t.s)*lambda.ss, D=number.of.deaths.in.founder)
  lower <- fun(lambda=1, N=N, delta=1+10^-10, t=(t.end-t.s)*lambda.ss, D=number.of.deaths.in.founder)
  if((upper>0 & lower>0) | (upper < 0 & lower < 0)){
    return(1000)
  }
  delta.founder <- uniroot(fun, interval=c((1+10^-10), 10), lambda=1, N=N, t=(t.end-(t.s ))*lambda.ss, D=number.of.deaths.in.founder)
  ## scale back
  delta.founder <- delta.founder$root*lambda.ss
  
  
  ## Next, compute the mutation distribution at t.s, i.e. the time point at which the driver mutation is acquired. Do this in a discretized way, as before. 

  if((log10(N)-1)>2){
    a <- c(1, 10, 25, 50, 75, 10^seq(2, log10(N)-1, 0.05))
    
  }else{
    a <- c(1, 10, 25, 50, 75, 100)
  }
  a <- sort(unique(round(c(a, seq(0.105, 1, 0.05)*N))))
  
  ## compute the mutational burden at t.s
  mutations.at.t.s <- mutational.burden(mu=mu, N=N, lambda.exp=lambda.exp, delta.exp=delta.exp, lambda.ss=lambda.ss, t.end=t.s, b=a)

  ## mutations per bin (cumulative --> discrete)
  mutations.at.t.s <- mutations.at.t.s - c(mutations.at.t.s[-1],0)
  
  ## Now, every mutation can be either present in the founder population, the selected population or both. Its fate is determined by these cases.
  ## For each bin at ts, the probability that the mutation is present in the selected clone, reads a/N. In this case, it will be present in all selected cells + putatively in a subset of founder cells
  
  bin.factor <- mutations.at.t.s/c(mutations.at.t.s[-1],0)
  bin.factor[is.infinite(bin.factor)] <- 10^8
  bin.factor[is.na(bin.factor)] <- 1
  
  upper.a <- c(a[-1], 2*N)
  
  ## next, compute the evolution of these mutations by taking the average from the lower and upper border of each bin
  
  mutations.at.t.end <- sapply(b, function(b){
    ## if b is bigger than the selected clone, it may contain mutations present in the founder cell population only (probability 1-a/N)
    ## or mutations present in both the selected clone and the founder population (probability a/N)
    
    if(b >= f.sel*N){  ## weighted average of the probability to drift from the lower and upper boundary of the interval to a size of at least b
      weighted.average <- apply(rbind(a, upper.a, bin.factor), 2, function(x){
        ## probability that the mutation is present in the selected clone computed at both borders of the bin
        prob.sel <- x[1]/N
        prob.sel.upper <- x[2]/N
        
        
        p <- (x[3]* (prob.sel*(1-p.a.b(a=x[1]-1, b=b-round(f.sel*N), lambda=lambda.ss, delta=delta.founder, t=t.end-t.s)) + ## mutation in both populations (-1 because 1 cel is the selected clone)
                       (1-prob.sel)*(1-p.a.b(a=x[1], b=b, lambda=lambda.ss, delta=delta.founder, t=t.end-t.s))) +## mutation only in founder cell 
                
                prob.sel.upper*(1-p.a.b(a=x[2]-1, b=b-round(f.sel*N), lambda=lambda.ss, delta=delta.founder, t=t.end-t.s)) + ## mutation in both populations
                (1-prob.sel.upper)*(1-p.a.b(a=x[2], b=b, lambda=lambda.ss, delta=delta.founder, t=t.end-t.s)) ## mutation only in founder cell 
              
        )/(x[3] + 1)
        
             
        return(p)
        
        
        
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
  
  ## Now, consider mutations acquired in the selected clone until t.s2
  
  size.first.clone.at.t.s2 <- exp((lambda.ss - s1*lambda.ss)*(t.s2 - t.s))
  if(size.first.clone.at.t.s2>N){
    size.first.clone.at.t.s2 <- N
  }
  f.sel.2 <- exp((lambda.ss - s1*s2*lambda.ss)*(t.end - t.s2))/N
  
  a <- a[a <= size.first.clone.at.t.s2]
  upper.a <- c(a[-1], 2*max(a))
  mutations.from.selected.clone <- mutations.noncritical.bd(lambda.ss, lambda.ss*s1, t.s2-t.s, mu, a, N=size.first.clone.at.t.s2)

  ## At t.s2 again, mutations could be present in the second subclone, in the first subclone or in both. Thus, consider the evolution of these mutations 
  ## mutations per bin 
  mutations.from.selected.clone <- mutations.from.selected.clone - c(mutations.from.selected.clone[-1],0)
  
  bin.factor <- mutations.from.selected.clone/c(mutations.from.selected.clone[-1],0)
  bin.factor[is.infinite(bin.factor)] <- 10^8
  bin.factor[is.na(bin.factor)] <- 1

  
  mutations.at.t.end.within.selected.clone <- sapply(b, function(b){
    
    if(b >= f.sel.2*N){ ## weighted average between lower and upper boundary of the interval
      
      weighted.average <- apply(rbind(a, upper.a, bin.factor), 2, function(x){
        ## probability that the mutation is present in the selected clone computed at both borders of the bin
        prob.sel <- x[1]/size.first.clone.at.t.s2
        prob.sel.upper <- x[2]/size.first.clone.at.t.s2
        (x[3]* (prob.sel*(1-p.a.b(a=x[1]-1, b=round(b-f.sel.2*N), lambda=lambda.ss, delta=s1*lambda.ss, t=t.end-t.s2)) + ## mutation in both populations (-1 because 1 cell is the selected clone)
                  (1-prob.sel)*(1-p.a.b(a=x[1], b=b, lambda=lambda.ss, delta=s1*lambda.ss, t=t.end-t.s2))) +## mutation only in founder cell 
            
            prob.sel.upper*(1-p.a.b(a=x[2]-1, b=round(b-f.sel.2*N), lambda=lambda.ss, delta=s1*lambda.ss, t=t.end-t.s2)) + ## mutation in both populations
            (1-prob.sel.upper)*(1-p.a.b(a=x[2], b=b, lambda=lambda.ss, delta=s1*lambda.ss, t=t.end-t.s2)) ## mutation only in founder cell 
          
        )/(x[3] + 1)
        
      })
      
      res <- sum(mutations.from.selected.clone*weighted.average)
      
    }else{ ## mutations present in at least b cells, where b <= mutations.from.selected.clone --> cumulative distribution from first selected clone +
      ## cumulative distribution from second selected clone, but in this case start with mutations.from.selected.clone at the lower end
      
      ## cumulative distribution from selected clone, but in this case start with f.sel.2 at the lower end as this smaller sizes are not compatible with the selected clone.
      
      
      weighted.average <- apply(rbind(a, upper.a, bin.factor), 2, function(x){
        
        prob.sel <- x[1]/size.first.clone.at.t.s2
        prob.sel.upper <- x[2]/size.first.clone.at.t.s2
        
        
        p <- (x[3]* (prob.sel*(1-p.a.b(a=x[1]-1, b=0, lambda=lambda.ss, delta=s1*lambda.ss, t=t.end-t.s2)) + ## mutation in both populations (-1 because 1 cel is the selected clone)
                       (1-prob.sel)*(1-p.a.b(a=x[1], b=b, lambda=lambda.ss, delta=s1*lambda.ss, t=t.end-t.s2))) +## mutation only in founder cell 
                
                prob.sel.upper*(1-p.a.b(a=x[2]-1, b=0, lambda=lambda.ss, delta=s1*lambda.ss, t=t.end-t.s2)) + ## mutation in both populations
                (1-prob.sel.upper)*(1-p.a.b(a=x[2], b=b, lambda=lambda.ss, delta=s1*lambda.ss, t=t.end-t.s2)) ## mutation only in founder cell 
              
        )/(x[3] + 1)
        
        return(p)
        
        
     })
      
      res <- sum(mutations.from.selected.clone*weighted.average)

    }
    res
  })
  
  
  ## Finally, we need to add up new mutations acquired during the expansion of both selected clones
  
  mutations.from.first.selected.clone <- mutations.noncritical.bd(lambda.ss, lambda.ss*s1, t.end-t.s2, mu, b, N=N)*(size.first.clone.at.t.s2-1)

  mutations.from.second.selected.clone <- mutations.noncritical.bd(lambda.ss, lambda.ss*s1*s2, t.end-t.s2, mu, b, N=N)
  
  mutations.from.founder.clone <- mutations.noncritical.bd(lambda.ss, delta.founder, t.end - t.s, mu, b, N0=N, N=(1-f.sel)*N)
  
  mutations.at.t.end <- mutations.at.t.end + mutations.at.t.end.within.selected.clone + mutations.from.first.selected.clone + mutations.from.second.selected.clone +
    mutations.from.founder.clone
  
  mutations.at.t.end
}


#########################################################################################################################################
####### Simulate whole genome sequencing
#########################################################################################################################################

#' Simulate sampling as in WGS by simulating the actual sampling (i.e. it's not the expected value but a sampling instance)
#' 
#' @param vafs vector of VAF at which cumulative mutation counts were measured
#' @param expected.mutations expected number of mutations at each VAF
#' @param depth sequencing depth
#' @param false.negative.per.vaf matrix with columns corresponding to the measured VAFs and rows corresponding to individual measurements of the false negative rate at this VAF.
#' @return A vector of simulated VAFs
#' @export

simulated.data <- function(vafs, expected.mutations, depth=90, sensitivity=T, false.negative.per.vaf, min.vaf=0.05){

  simulated.vafs <- c()
  sim <- apply(rbind(vafs, expected.mutations), 2, function(x){
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


#' Compute the expected number of mutations after WGS (note: runs slowly!)
#' 
#' @param M vector of cumulative mutation counts with VAF as given in 'vaf' 
#' @param vaf vector of variant allele frequencies for which M is reported
#' @param depth sequencing depth
#' @param false.negative.per.vaf estimated fraction of false negatives per vaf; columns correspond to VAFs as in vaf, rows to individual measurements
#' @return A vector of simulated VAFs
#' @export

ExpectedNumberAfterWGS <- function(M, vafs, depth, false.negative.per.vaf, min.vaf=0.05){
  

  ## counts per interval
  M <- M - c(M[-1], 0)
  
  ## additional unexplained error
  false.negative.per.vaf. <- t(t(false.negative.per.vaf) - sapply(seq(0.05, 1, 0.01), function(x){pbinom(0.05*depth, size=depth, prob=x)} ))
  false.negative.per.vaf.[false.negative.per.vaf. < 0 ] <- 0
  
  ## expected fraction to be recovered per bin at lower end: conditional probability on coverage (Poisson distributed)
  
  expected.number.to.measure <- function(expected.number, vafs, depth, vaf.of.interest){
    sum(sapply(seq(1, length(expected.number)), function(i){
      expected.number[i]*
        sum(sapply(seq(0, 1000), function(coverage){
          dpois(coverage, depth)*pbinom(q=round(vaf.of.interest*coverage), size=coverage, prob=vafs[i]/2, lower.tail = F)
        }))*(1-mean(false.negative.per.vaf.[,which.min((vaf.of.interest - seq(0.05, 1, 0.01))^2):ncol(false.negative.per.vaf.)]))
    }))
    
  }
  
  sim.vafs <- sapply( vafs, function(vaf){expected.number.to.measure(M, vafs, depth, vaf)})
  
  sim.vafs
  
}
