#' Propensity function for "S"tem cells and 2 different mature cells (Type 1 and Type 2)
#' 
#' This function computes the propensities for a linear, hierarchical system of cell division, death and differentiation
#' @param cell.count a vector with cell numbers (S, T1, T2)
#' @param parms a named parameter vector containing the rates of cell division (lambda.s/lambda.t1/lambda.t2) and differentiation (alpha.s1/alpha.s2) or loss (delta.s/delta.t1/delta.t2), where "s" and "t" indicate stem and type (1 or 2) cells, respectively
#' @return the propensities of the possible reactions
#' @export

.props.het <- function(cell.count, parms) {
  with(as.list(c(cell.count, parms)), {
    
    # division of S
    div.s <- lambda.s*cell.count[1]
    # differentiation of S into T1
    loss.s1 <- alpha.s1*cell.count[1]
    # differentiation of S into T2
    loss.s2 <- alpha.s2*cell.count[1]
    # division of Type 1 cells
    div.t1 <- lambda.t1*cell.count[2]
    # division of Type 2 cells
    div.t2 <- lambda.t2*cell.count[3]
    # differentiation of Type 1 cells
    loss.t1 <- delta.t1*cell.count[2]
    # differentiation of Type 2 cells
    loss.t2 <- delta.t2*cell.count[3]
    # death of S
    death.s <- delta.s*cell.count[1]

    res <- c(div.s,  loss.s1, loss.s2, div.t1, loss.t1, div.t2, loss.t2, death.s)
    return(res)
    
  })
}



#' Tree simulation
#' 
#' Simulate the phylogenetic tree of a physiological population that grows in 2 regimes: initial exponential expansion and subsequent homeostasis. The function can simuate stem cells only or stem and progenitor cells
#' @param parms.exp expansion parameters, vector that must contain lambda.s (stem cell division rate), lambda.t1 (type 1 cell division rate), lambda.t2 (type 2 cell division rate), alpha.s1 (stem cell differentiation rate into type 1 cells), alpha.s2 (stem cell division rate into type 2 cells), delta.s (stem cell loss rate), delta.t1 (type 1 cell differentiation rate), delta.t2 (type 2 cell differentiation rate)
#' @param parms.steady as parms.exp, but for the homeostatic phase
#' @param time.max maximal simulation time
#' @param time.samples time points at which simulation results are stored
#' @param N number of stem cells during homeostasis
#' @param NT1 number of type 1 cells during homeostasis
#' @param NT2 number of type 2 cells during homeostasis
#' @param mut.rate average number of mutation per division and daughter cell
#' @param mutation.mode should mutations be "constant" or "Binomial"ly distributed?
#' @param driver.mode string; "fixed_time" if driver is acquired at a fixed time point, "random", if driver is randomly acquired governed by the parameters mu_D and s
#' @param t.driver integer vector, the time point(s) at which the driver is acquired if `mode` is "fixed_time"
#' @param mother integer vector containing the mother clone of each daughter clone; will be ignored if `driver.mode` is "random".
#' @param mut.rate.D integer; driver mutation rate per cell division
#' @param s.shape shape parameter of the gamma distribution from which the selective advantage of a new driver is drawn; if `driver.mode`== "fixed_time" and `t.driver` is a vector of multiple events, `s.shape` must be a vector of equal length as `t.driver`
#' @param s.rate rate parameter of the gamma distribution from which the selective advantage of a new driver is drawn; if `driver.mode`== "fixed_time" and `t.driver` is a vector of multiple events, `s.shape` must be a vector of equal length as `t.driver`
#' @param tau tau leaping parameter: how many steps should be merged during homeostasis? Defaults to 1, no leaping.
#' @param report.at.f frequency of a selected clone at which the simulation should be stopped. Defaults to NA - don't stop at a certain frequency. 
#' @return a list of state.lists at the desired time samples
#' @export

gillespie.sim.s.t1.t2 <- function(parms.exp, parms.steady, time.max=50, time.samples=c(0, 5, 25, 50), N=1000, NT1=8000, NT2=2000, report.at.f = NA,
                              mut.rate = 3, mutation.mode="Binomial", driver.mode="random", t.driver = NA, mother = NULL, mut.rate.D = 0, s.shape = 1.5, s.rate = 35, tau = 1){
  
  if(max(time.samples) < time.max){
    time.samples <- c(time.samples, time.max, time.max + 1)
  }
  
  if(driver.mode=="fixed_time"){
    mut.rate.D=0 ## no random drivers
    if(length(s.shape) != length(s.rate) | length(t.driver) != length(s.rate) | length(t.driver) != length(s.shape) | length(t.driver) != length(mother)){
      stop("Error: please provide equal number of `t.driver`, `s.rate` and `s.shape`.")
    }
  }
  
  if(! any(c("lambda.s", "lambda.t1", "lambda.t2", "alpha.s1", "delta.s", "delta.t1", "delta.t2")   %in% names(parms.exp))){
    stop("Error: parms.exp must contain values for `lambda.s`, `lambda.t1`, `lambda.t2`, `alpha.s1`, `alpha.s2`, `delta.s`, `delta.t1` and `delta.t2`.")
  }
  if(! any(c("lambda.s", "lambda.t1", "lambda.t2", "alpha.s1", "alpha.s2", "delta.s", "delta.t1", "delta.t2")   %in% names(parms.steady))){
    stop("Error: parms.steady must contain values for `lambda.s`, `lambda.t1`, `lambda.t2`, `alpha.s1`, `alpha.s2`, `delta.s`, `delta.t1` and `delta.t2`.")
  }
  
  if(driver.mode!="fixed_time"){
    s.shape <- s.shape[1]
    s.rate <- s.rate[1]
  }
  
  #### output
  
  res.list <- list()
  
  #### initialize the system
  
  tree <- .initialize.sim.s.p()
  
  #### 1. expansion phase, until reaching 2*N stem cells, only simulate stem cells
  
  parms.exp.I <- parms.exp
  parms.exp.I[c("lambda.p")] <- 0
  #parms.exp.I[c("delta.s")] <- 0
  parms.exp.I[c("alpha.s1")] <- 0
  parms.exp.I[c("alpha.s2")] <- 0
  
  n.cells <-  c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2), "3"=sum(tree$tip.class==3))
  
  ## if a hierarchical system is simulated, produce 2*N cells of which 50% differentiate. Else, just N stem cells
  target.size <- ifelse(NT1>0 | NT2>0, 2*N, N)
  expansion.time <- 0
  while(n.cells["1"] < target.size ){
    
    if(sum(n.cells)%%1000==0){
      print(paste("Expansion phase: population has grown to", sum(n.cells), "cells"))
    }
    
    ## reinitialize if population goes extinct
    if(n.cells["1"]==0){
      tree <- .initialize.sim.s.p()
      n.cells <-  c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2),
                    "3"=sum(tree$tip.class==3))
    }
    
    ## sample the propensities for the reaction to happen
    props <- .props.het(n.cells, parms = parms.exp.I)
    
    # update expansion time
    dt <- rexp(1, rate = sum(props))
    expansion.time <- expansion.time + dt
    
    reaction <- which(cumsum(props)/sum(props) >= runif(1,0,1))
    reaction <- reaction[1]
    
    if(reaction %in% c(1,2,3,8)){
      cell.type <- 1
    }else if(reaction %in% c(4,5)){
      cell.type <- 2
    }else if(reaction %in% c(6,7)){
      cell.type <- 3
    }

    ## symmetric division
    if(reaction %in% c(1,4,6)){
      
      tree <- .division.s.p(cell.type=cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate)
      
    }else if(reaction==2){ ## differentiation of stem cell into type 1 cell
      
      tree <- .differentiation.s.p(tree = tree, type = 2)
      
    }else if(reaction==3){ ## differentiation of stem cell into type 2 cell
      
      tree <- .differentiation.s.p(tree = tree, type = 3)
      
    }else if(reaction %in% c(5, 7, 8)){## cell loss
      ## which mutations are lost?
      
      tree <- .loss(cell.type=cell.type, tree = tree)
      
    }
    n.cells <- c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2),
                 "3"=sum(tree$tip.class==3))
  }
  
  tree$expansion.time = expansion.time
  
  #### 2. simulate type 1 and type 2 cells by letting 50% of the stem cells differentiate
  
  parms.exp.II <- parms.exp
  parms.exp.II[c("lambda.s")] <- 0
  parms.exp.II[c("delta.t1")] <- 0
  parms.exp.II[c("lambda.t1")] <- 0
  parms.exp.II[c("delta.t2")] <- 0
  parms.exp.II[c("lambda.t2")] <- 0
  parms.exp.II[c("delta.s")] <- 0
  
  while(n.cells["1"] > N ){
    
    if(sum(n.cells[c(2,3)])%%1000==0){
      print(paste("Expansion phase: Progenitor population has grown to", sum(n.cells[c(2,3)]), "cells"))
    }
    
    ## reinitialize if population goes extinct
    if(n.cells[1]==0){
      tree <- .initialize.sim.s.p()
    }
    
    ## sample the propensities for the reaction to happen
    props <- .props.het(n.cells, parms = parms.exp.II)
    
    reaction <- which(cumsum(props)/sum(props) >= runif(1,0,1))
    reaction <- reaction[1]
    
    if(reaction %in% c(1,2,3,8)){
      cell.type <- 1
    }else if(reaction %in% c(4,5)){
      cell.type <- 2
    }else if(reaction %in% c(6,7)){
      cell.type <- 3
    }    
    
    ## symmetric division
    if(reaction %in% c(1,4,6)){
      
      tree <- .division.s.p(cell.type=cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate)
      
    }else if(reaction==2){ ## differentiation of stem cell into type 1 cell
      
      tree <- .differentiation.s.p(tree = tree, type = 2)
      
    }else if(reaction==3){ ## differentiation of stem cell into type 2 cell
      
      tree <- .differentiation.s.p(tree = tree, type = 3)
      
    }else if(reaction %in% c(5, 7, 8)){## cell loss
      ## which mutations are lost?
      
      tree <- .loss(cell.type=cell.type, tree = tree)
      
    }    
    n.cells <- c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2),
                 "3"=sum(tree$tip.class==3))
    
  }
  
  #### 3. Let type 1 divide until reaching NT1
  
  parms.exp.III <- parms.exp
  parms.exp.III[c("lambda.s")] <- 0
  parms.exp.III[c("alpha.s1")] <- 0
  parms.exp.III[c("alpha.s2")] <- 0
  parms.exp.III[c("delta.s")] <- 0
  parms.exp.III[c("delta.t1")] <- 0
  parms.exp.III[c("delta.t2")] <- 0
  parms.exp.III[c("lambda.t2")] <- 0
  
  while( n.cells[2] < NT1 ){
    
    
    ## reinitialize if population goes extinct
    if(n.cells[1]==0){
      tree <- .initialize.sim.s.p()
    }
    
    ## sample the propensities for the reaction to happen
    props <- .props.het(n.cells, parms = parms.exp.III)
    
    reaction <- which(cumsum(props)/sum(props) >= runif(1,0,1))
    reaction <- reaction[1]
    
    if(reaction %in% c(1,2,3,8)){
      cell.type <- 1
    }else if(reaction %in% c(4,5)){
      cell.type <- 2
    }else if(reaction %in% c(6,7)){
      cell.type <- 3
    }    
    
    ## symmetric division
    if(reaction %in% c(1,4,6)){
      
      tree <- .division.s.p(cell.type=cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate)
      
    }else if(reaction==2){ ## differentiation of stem cell into type 1 cell
      
      tree <- .differentiation.s.p(tree = tree, type = 2)
      
    }else if(reaction==3){ ## differentiation of stem cell into type 2 cell
      
      tree <- .differentiation.s.p(tree = tree, type = 3)
      
    }else if(reaction %in% c(5, 7, 8)){## cell loss
      ## which mutations are lost?
      
      tree <- .loss(cell.type=cell.type, tree = tree)
      
    }  
    n.cells <- c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2),
                 "3"=sum(tree$tip.class==3))
  }
  
  #### 4. Let type 2 divide until reaching NT2
  
  parms.exp.III <- parms.exp
  parms.exp.III[c("lambda.s")] <- 0
  parms.exp.III[c("alpha.s1")] <- 0
  parms.exp.III[c("alpha.s2")] <- 0
  parms.exp.III[c("delta.s")] <- 0
  parms.exp.III[c("delta.t1")] <- 0
  parms.exp.III[c("delta.t2")] <- 0
  parms.exp.III[c("lambda.t1")] <- 0
  
  while( n.cells[3] < NT2 ){
    
    
    ## reinitialize if population goes extinct
    if(n.cells[1]==0){
      tree <- .initialize.sim.s.p()
    }
    
    ## sample the propensities for the reaction to happen
    props <- .props.het(n.cells, parms = parms.exp.III)
    
    reaction <- which(cumsum(props)/sum(props) >= runif(1,0,1))
    reaction <- reaction[1]
    
    if(reaction %in% c(1,2,3,8)){
      cell.type <- 1
    }else if(reaction %in% c(4,5)){
      cell.type <- 2
    }else if(reaction %in% c(6,7)){
      cell.type <- 3
    }    
    
    ## symmetric division
    if(reaction %in% c(1,4,6)){
      
      tree <- .division.s.p(cell.type=cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate)
      
    }else if(reaction==2){ ## differentiation of stem cell into type 1 cell
      
      tree <- .differentiation.s.p(tree = tree, type = 2)
      
    }else if(reaction==3){ ## differentiation of stem cell into type 2 cell
      
      tree <- .differentiation.s.p(tree = tree, type = 3)
      
    }else if(reaction %in% c(5, 7, 8)){## cell loss
      ## which mutations are lost?
      
      tree <- .loss(cell.type=cell.type, tree = tree)
      
    }  
    n.cells <- c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2),
                 "3"=sum(tree$tip.class==3))
  }
  
  res.list[["0"]] <- tree
  
  #### 5. steady state phase: now, compute at discrete times
  time <- 0
  next.time.check <- 0
  next.f.check <- report.at.f[1]
  
  ## the largest selected clone
  if(any(tree$sel.adv.>1)){
    max.clone <- max(table(tree$sel.adv.[tree$sel.adv.!=1]))
  }else{
    max.clone <- 0
  }
  
  # initialize the driver event such that in the case of multiple drivers each event will be associated with the correct selective advantage.
  if(driver.mode=="fixed_time"){
    t.this.driver <- t.driver[1]
    s.shape.this.driver <- s.shape[1]
    s.rate.this.driver <- s.rate[1]
    mother.this.driver <- mother[1]
  }
  
  while(time < time.max){
    
    if(sum(n.cells)==0){break}
    ## sample the propensities for the reaction to happen
    props <- .props.het(n.cells, parms = parms.steady)
    
    ## adjust dt according to tau
    dt <- rexp(1, rate = sum(props)/tau)
    time <- time + dt
    
    if(round(time) %in% time.samples & time > next.time.check){
      res.list[[as.character(round(time))]] <- tree
      print(paste("Simulation time:", time))
      next.time.check <- time.samples[time.samples > time][1]
    }
    
    if(!is.na(next.f.check)){
      if(max.clone/sum(n.cells) > next.f.check){
        res.list[[as.character(round(time))]] <- tree
        next.f.check <- report.at.f[report.at.f > next.f.check][1]
        if(length(next.f.check)==0){
          next.f.check <- NA
        }
      }
    }
    
    reactions <- table(sapply(runif(tau, 0, 1), function(x){which(cumsum(props)/sum(props) >= x)[1]}))
    
    for(i in 1:length(reactions)){
      reaction <- names(reactions)[i]
      nr.of.reactions <- reactions[i]

      if(reaction %in% c(1,2,3,8)){
        cell.type <- 1
      }else if(reaction %in% c(4,5)){
        cell.type <- 2
      }else if(reaction %in% c(6,7)){
        cell.type <- 3
      }   
      
      ## symmetric division
      if(reaction %in% c(1,4,6)){
        
        tree <- .division.s.p(cell.type=cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate, nr = nr.of.reactions)
        
      }else if(reaction==2){ ## differentiation of stem cell into type 1 cell
        
        tree <- .differentiation.s.p(tree = tree, type = 2, nr = nr.of.reactions)
        
      }else if(reaction==3){ ## differentiation of stem cell into type 2 cell
        
        tree <- .differentiation.s.p(tree = tree, type = 3, nr = nr.of.reactions)
        
      }else if(reaction %in% c(5, 7, 8)){## cell loss
        ## which mutations are lost?
        
        tree <- .loss(cell.type=cell.type, tree = tree, nr = nr.of.reactions)
        
      }  
    
    }
    
    ## if drivers are introduced at a fixed time point rather than randomly: 
    if(driver.mode=="fixed_time"){
      
      for(birth.date in t.driver){
        if(birth.date < time & sum(tree$clone.id == (which(t.driver == birth.date) + 1)) == 0){
          mother.this.driver <- mother[which(t.driver == birth.date)]
          s.shape.this.driver <- s.shape[which(t.driver == birth.date)]
          s.rate.this.driver <- s.rate[which(t.driver == birth.date)]
          ## randomly select the cell to acquire the driver mutant
          if(length(tree$tip.class[tree$tip.class==1 & 
                                   tree$clone.id == mother.this.driver]) == 0){next}
          random.cell <- sample(seq(1,length(tree$tip.class[tree$tip.class==1 & 
                                                              tree$clone.id == mother.this.driver])), 1)
          tree[tree$tip.class==1 & 
                 tree$clone.id == mother.this.driver]$sel.adv.[random.cell] <- 1 + s.shape.this.driver/s.rate.this.driver
          tree[tree$tip.class==1 & 
                 tree$clone.id == mother.this.driver]$clone.id[random.cell]  <- 1 + which(t.driver == birth.date)
          
        }
      }
      
    }
    
    if(any(tree$sel.adv.>1)){
      max.clone <- max(table(tree$sel.adv.[tree$sel.adv.!=1]))
    }else{
      max.clone <- 0
    }
    
    n.cells <- c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2),
                 "3"=sum(tree$tip.class==3))
  }
  
  return(res.list)
}

