#' Propensity function for 2 compartments ("S"tem cells and "P"rogenitors)
#' 
#' This function computes the propensities for a linear, hierarchical system of cell division, death and differentiation
#' @param cell.count a vector with cell numbers (S, P)
#' @param parms a named parameter vector containing the rates of cell division ("lambda.s"/"lambda.p") and differentiation ("alpha.s") or loss ("delta.s"/"delta.p"), where "s" and "p" indicate stem and progenitor cells, respectively
#' @return the propensities of the possible reactions
#' @export

.props.2c <- function(cell.count, parms) {
  with(as.list(c(cell.count, parms)), {
    
    # division of S
    div.s <- lambda.s*cell.count[1]
    # differentiation of S
    loss.s <- alpha.s*cell.count[1]
    # division of P
    div.p <- lambda.p*cell.count[2]
    # differentiation of P
    loss.p <- delta.p*cell.count[2]
    # death of S
    death.s <- delta.s*cell.count[1]
    # asymmetric division (only if lambda.a is defined)
    if(exists("lambda.a")){
      as.div.s <- lambda.a*cell.count[1]
      res <- c(div.s,  loss.s, div.p, loss.p, death.s, as.div.s)
      return(res)
    }else{
      res <- c(div.s,  loss.s, div.p, loss.p, death.s)
      return(res)
    }
    

  })
}

#' Acquisition of neutral mutations
#' 
#' This function simulates the acquisition of new mutations
#' @param mut.rate integer; average number of mutations per division and daughter cell
#' @param mutation.mode character vector; should the number be based on a binomial distribution ("Binomial") or should a constant number be introduced ("constant"); defaults to "Binomial".
#' @return the number of new mutations per daughter cell
#' @export

.mutation.acquisition <- function(mut.rate, mutation.mode = "Binomial"){
  if(mutation.mode=="constant"){
    n.mut <- c(mut.rate,mut.rate)
  }else{
    n.mut <- rbinom(n = 2, size = 10^9, prob = mut.rate*10^-9)
  }
  return(n.mut)
}

#' System initialization
#' 
#' This function initializes the simulation
#' @param N integer; number of stem cells at homeostasis
#' @param NP integer; number of progenitor cells at homeostasis
#' @param mut.rate integer; number of mutations per cell division and daughter cell
#' @param time.max time of simulation
#' @param parms.steady named vector of steady state parameters; must contain "lambda.s"/"lambda.p" the stem and progenitor division rate; "delta.s"/"delta.p" the stem and progenitor loss rate and "alpha.s" the stem cell differentiation rate
#' @return returns the state list; initialized with a single stem cell and a single mutation
#' @export

.initialize.sim.s.p <- function(){
  
  tree <- list(edge=matrix(c(2, 1), nrow=1),
               edge.length=2,
               tip.label=1,
               Nnode=1,
               tip.class=1,
               sel.adv.=1)
  
  class(tree) <- "phylo"
  
  return(tree)
}

#' Cell division
#' 
#' Simulate division of a cell together with acquisition of neutral and driver mutations
#' @param cell.type the cell type that is to divide (1, stem cell; 2, progenitor cell)
#' @param tree object of class phylo; the current tree. Needs the following additional list elements: `tip.class`, specifying the cell type and `sel.adv.`, specifying the selective advantage of each tip
#' @param mut.rate.D integer; driver mutation rate per cell division
#' @param s.shape shape parameter of the gamma distribution from which the selective advantage of a new driver is drawn
#' @param s.rate rate parameter of the gamma distribution from which the selective advantage of a new driver is drawn
#' @param mutation.mode should the number of new mutations be "constant" or "Binomial"ly distributed?
#' @param mut.rate average number of neutral mutations per daughter cell and division
#' @param symmetric logical, is the division a symmetric division (2 new daughter cells) or an asymmetric division (1 daughter cell differentiates); defaults to T
#' @param nr number of reactions to simulate; defaults to 1
#' @return the updated tree

.division.s.p <- function(cell.type, tree, mut.rate.D = 0, s.shape = 1.5, s.rate = 35, mutation.mode, mut.rate, symmetric = T, nr = 1){
  
  if(symmetric==F & cell.type==2){ ## asymmetric division can only take place in stem cells
    return(tree)
  }
  
  ## the living cells are the tips
  living.cells <- tree$tip.label[tree$tip.class==cell.type]
  living.cells.sel.adv <- tree$sel.adv.[tree$tip.class==cell.type]

  dividing.cell <- sample(living.cells, nr, prob = living.cells.sel.adv/sum(living.cells.sel.adv))
  
  ## sample the new mutations
  n.mut <- as.vector(sapply(1:nr, function(x){.mutation.acquisition(mut.rate = mut.rate, mutation.mode = mutation.mode)}))
  
  ## selective advantage of daughters
  new.s <- as.vector(sapply(1:nr, function(x){.driver.mutation.acquisition(mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate)}))
  
  ## add the 2 daughters together with their mutations and selective advantage
  
  new.tree <- phangorn::add.tips(tree = tree, tips = max(tree$tip.label) + seq(1, 2*nr), where = dividing.cell, edge.length = n.mut)
  
  new.tree$sel.adv. <- c(tree$sel.adv., c(tree$sel.adv.[dividing.cell]+new.s[1:nr], tree$sel.adv.[dividing.cell]+new.s[(nr+1):(2*nr)]))
  
  if(symmetric){
    new.tree$tip.class <- c(tree$tip.class, rep(cell.type, 2*nr))
  }else if (cell.type==1){
    new.tree$tip.class <- c(tree$tip.class, rep(c(cell.type, 2), each=nr))
  }
  

  ## remove the dividing cell from the tips
  
  new.tree$sel.adv. <- new.tree$sel.adv.[-which(new.tree$tip.label %in% dividing.cell)]
  new.tree$tip.class <- new.tree$tip.class[-which(new.tree$tip.label %in% dividing.cell)]
  
  new.tree <- ape::drop.tip(new.tree, dividing.cell, trim.internal = T, collapse.singles = F)
  
  new.tree$tip.label <- sort(setdiff(new.tree$edge[,2], new.tree$edge[,1]))
  
  tree <- new.tree
  
  return(tree)
  
}

#' Cell differentiation
#' 
#' Simulate differentiation of a stem cell into a progenitor cell
#' @param tree object of class phylo; the current tree. Needs the following additional list elements: tip.class, specifying the cell type and sel.adv., specifying the selective advantage of each tip
#' @param nr number of reactions to simulate; defaults to 1
#' @return the updated tree

.differentiation.s.p <- function(tree, nr=1){
  
  ## the living stem cells are the tips of class 1
  living.cells <- tree$tip.label[tree$tip.class==1]

  ## randomly sample a stem cell that differentiates
  differentiating.cell <- sample(living.cells, nr)
  
  ## update the cell's type to progenitor cell
  tree$tip.class[tree$tip.label %in% differentiating.cell] <- 2
  
  return(tree)
  
}

#' Cell loss
#' 
#' Simulate loss of a stem cell or a progenitor cell
#' @param cell.type integer; which cell type is dividing? 1, stem cell, 2, progenitor cell
#' @param tree object of class phylo; the current tree. Needs the following additional list elements: `tip.class`, specifying the cell type and `sel.adv.`, specifying the selective advantage of each tip
#' @param nr number of reactions to simulate; defaults to 1
#' @return the updated tree

.loss <- function(cell.type, tree, nr=1){
  
  ## the living cells are the tips
  living.cells <- tree$tip.label[tree$tip.class==cell.type]
  dying.cell <- sample(living.cells, nr)

  ## remove the dying cell from the tips
  new.tree <- tree
  new.tree$sel.adv. <- new.tree$sel.adv.[-which(new.tree$tip.label %in% dying.cell)]
  new.tree$tip.class <- new.tree$tip.class[-which(new.tree$tip.label %in% dying.cell)]
  
  new.tree <- ape::drop.tip(new.tree, dying.cell, trim.internal = T, collapse.singles = F)
  
  new.tree$tip.label <- sort(setdiff(new.tree$edge[,2], new.tree$edge[,1]))
  tree <- new.tree
  
  return(tree)
}

#' Driver mutation acquisition
#' 
#' Simulate the acquisition of new driver mutations in both daughter cells
#' @param mut.rate.D integer; driver mutation rate per cell division
#' @param s.shape shape parameter of the gamma distribution from which the selective advantage of a new driver is drawn
#' @param s.rate rate parameter of the gamma distribution from which the selective advantage of a new driver is drawn
#' @return the additional selective advantage acquired during the cell division

.driver.mutation.acquisition <- function(mut.rate.D, s.shape, s.rate){
  
  ## initialize the selective advantage of both daughters with 0
  s <- c(0, 0)
  
  for(i in c(1,2)){
    ## sample for both daughter cells whether they obtain a new driver
    new.driver <- sample(x=c(0, 1), size=1, prob=c(1-mut.rate.D, mut.rate.D))
    
    if(new.driver){
      ## sample s such that it's at least 5% with parameters taken from Mitchell et al.
      while(s[i] < 0.05){
        s[i] <- rgamma(1, shape=s.shape, rate=s.rate)
      }
    }
  }
  
  return(s)   
}


#' Tree simulation
#' 
#' Simulate the phylogenetic tree of a physiological population that grows in 2 regimes: initial exponential expansion and subsequent homeostasis. The function can simuate stem cells only or stem and progenitor cells
#' @param parms.exp expansion parameters, vector that must contain `lambda.s` (stem cell division rate), `lambda.p` (progenitor cell division rate), `alpha.s` (stem cell differentiation rate), `delta.s` (stem cell loss rate), `delta.p` (progenitor cell differentiation rate)
#' @param parms.steady same as `parms.exp`, but for the homeostatic phase
#' @param time.max maximal simulation time
#' @param time.samples time points at which simulation results are stored
#' @param N number of stem cells during homeostasis
#' @param NP number of progenitor cells during homeostasis
#' @param mut.rate average number of mutation per division and daughter cell
#' @param mutation.mode should mutations be "constant" or "Binomial"ly distributed?
#' @param driver.mode string; "fixed_time" if driver is acquired at a fixed time point, "random", if driver is randomly acquired governed by the parameters `mu_D` and `s`
#' @param t.driver integer, the time point at which the driver is acquired if mode is "fixed_time"
#' @param mut.rate.D integer; driver mutation rate per cell division
#' @param s.shape shape parameter of the gamma distribution from which the selective advantage of a new driver is drawn. If the function is run with `driver.mode="fixed"`, the selective advantage is not randomyl drawn but computed as s.shape/d.rate.
#' @param s.rate rate parameter of the gamma distribution from which the selective advantage of a new driver is drawn
#' @param tau tau leaping parameter: how many steps should be merged during homeostasis? Defaults to 1, no leaping.
#' @param report.at.f frequency of a selected clone at which the simulation should be stopped. Defaults to NA - don't stop at a certain frequency. 
#' @return a list of state.lists at the desired time samples
#' @export

gillespie.sim.s.p <- function(parms.exp, parms.steady, time.max=50, time.samples=c(0, 5, 25, 50), N=1000, NP=10000, report.at.f = NA,
                              mut.rate = 3, mutation.mode="Binomial", driver.mode="random", t.driver = NA, mut.rate.D = 0, s.shape = 1.5, s.rate = 35, tau = 1){
  
  if(max(time.samples) < time.max){
    time.samples <- c(time.samples, time.max, time.max + 1)
  }
  
  if(driver.mode=="fixed_time"){
    mut.rate.D=0 ## no random drivers
  }
  
  #### output
  
  res.list <- list()
  
  #### initialize the system
  
  tree <- .initialize.sim.s.p()
  
  #### 1. expansion phase, until reaching 2*N stem cells, only simulate stem cells
  
  parms.exp.I <- parms.exp
  parms.exp.I[c("lambda.p")] <- 0
  #parms.exp.I[c("delta.s")] <- 0
  parms.exp.I[c("alpha.s")] <- 0
  
  n.cells <-  c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2))
  
  ## if a hierarchical system is simulated, produce 2*N cells of which 50% differentiate. Else, just N stem cells
  target.size <- ifelse(NP>0, 2*N, N)
  expansion.time <- 0
  while(n.cells["1"] < target.size ){
    
    if(sum(n.cells)%%1000==0){
      print(paste("Expansion phase: population has grown to", sum(n.cells), "cells"))
    }
    
    ## reinitialize if population goes extinct
    if(n.cells["1"]==0){
      tree <- .initialize.sim.s.p()
      n.cells <-  c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2))
    }
    
    ## sample the propensities for the reaction to happen
    props <- .props.2c(n.cells, parms = parms.exp.I)
    
    # update expansion time
    dt <- rexp(1, rate = sum(props))
    expansion.time <- expansion.time + dt
    
    reaction <- which(cumsum(props)/sum(props) >= runif(1,0,1))
    reaction <- reaction[1]
    
    cell.type <- ifelse(reaction %in% c(1,2,5,6), 1, 2)
    
    ## symmetric division
    if(reaction %in% c(1,3)){
      
      tree <- .division.s.p(cell.type=cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate)
      
    }else if(reaction==6){ ## asymmetric division of a stem cell
      
      tree <- .division.s.p(cell.type = cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate, symmetric = F)
      
    }else if(reaction==2){ ## differentiation of stem cell into progenitor cell
      
      tree <- .differentiation.s.p(tree = tree)
      
    }else if(reaction==4){## loss of progenitor cell
      ## which mutations are lost?
      
      tree <- .loss(cell.type=cell.type, tree = tree)
      
    }else if(reaction==5){## loss of stem cell
      ## which mutations are lost?
      
      tree <- .loss(cell.type=cell.type, tree = tree)
      
    }
    n.cells <- c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2))
  }
  
  tree$expansion.time = expansion.time
  
  #### 2. simulate progenitor cells by letting 50% of the stem cells differentiate
  
  parms.exp.II <- parms.exp
  parms.exp.II[c("lambda.s")] <- 0
  parms.exp.II[c("delta.p")] <- 0
  parms.exp.II[c("lambda.p")] <- 0
  parms.exp.II[c("delta.s")] <- 0
  
  while(n.cells["1"] > N ){
    
    if(sum(n.cells[2])%%1000==0){
      print(paste("Expansion phase: Progenitor population has grown to", n.cells[2], "cells"))
    }
    
    ## reinitialize if population goes extinct
    if(n.cells[1]==0){
      tree <- .initialize.sim.s.p()
    }
    
    ## sample the propensities for the reaction to happen
    props <- .props.2c(n.cells, parms = parms.exp.II)
    
    reaction <- which(cumsum(props)/sum(props) >= runif(1,0,1))
    reaction <- reaction[1]

    cell.type <- ifelse(reaction %in% c(1,2,5,6), 1, 2)
    
    ## symmetric division
    if(reaction %in% c(1,3)){
      
      tree <- .division.s.p(cell.type=cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate)
      
    }else if(reaction==6){ ## asymmetric division of a stem cell
      
      tree <- .division.s.p(cell.type = cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate, symmetric = F)
      
    }else if(reaction==2){ ## differentiation of stem cell into progenitor cell
      
      tree <- .differentiation.s.p(tree = tree)
      
    }else if(reaction==4){## loss of progenitor cell
      ## which mutations are lost?
      
      tree <- .loss(cell.type=cell.type, tree = tree)
      
    }else if(reaction==5){## loss of stem cell
      ## which mutations are lost?
      
      tree <- .loss(cell.type=cell.type, tree = tree)
      
    }
    
    n.cells <- c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2))
    
  }
  
  #### 3. Let progenitors divide until reaching NP
  
  parms.exp.III <- parms.exp
  parms.exp.III[c("lambda.s")] <- 0
  parms.exp.III[c("alpha.s")] <- 0
  parms.exp.III[c("delta.s")] <- 0
  parms.exp.III[c("delta.p")] <- 0
  
  while( n.cells[2] < NP ){
    
    
    ## reinitialize if population goes extinct
    if(n.cells[1]==0){
      tree <- .initialize.sim.s.p()
    }
    
    ## sample the propensities for the reaction to happen
    props <- .props.2c(n.cells, parms = parms.exp.III)
    
    reaction <- which(cumsum(props)/sum(props) >= runif(1,0,1))
    reaction <- reaction[1]
    
    cell.type <- ifelse(reaction %in% c(1,2,5,6), 1, 2)
    
    ## symmetric division
    if(reaction %in% c(1,3)){
      
      tree <- .division.s.p(cell.type=cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate)
      
    }else if(reaction==6){ ## asymmetric division of a stem cell
      
      tree <- .division.s.p(cell.type = cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate, symmetric = F)
      
    }else if(reaction==2){ ## differentiation of stem cell into progenitor cell
      
      tree <- .differentiation.s.p(tree = tree)
      
    }else if(reaction==4){## loss of progenitor cell
      ## which mutations are lost?
      
      tree <- .loss(tree = tree)
      
    }else if(reaction==5){## loss of stem cell
      ## which mutations are lost?
      
      tree <- .loss(cell.type = cell.type, tree = tree)
      
    }
    n.cells <- c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2))
  }
  
  res.list[["0"]] <- tree
  
  #### 4. steady state phase: now, compute at discrete times
  time <- 0
  next.time.check <- 0
  next.f.check <- report.at.f[1]
  
  ## the largest selected clone
  if(any(tree$sel.adv.>1)){
    max.clone <- max(table(tree$sel.adv.[tree$sel.adv.!=1]))
  }else{
    max.clone <- 0
  }
  
  while(time < time.max){

    if(sum(n.cells)==0){break}
    
    ## sample the propensities for the reaction to happen
    props <- .props.2c(n.cells, parms = parms.steady)
    
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

    ## asymmetric division must take place prior to death, otherwise run into problems
    reactions <- reactions[intersect(c("1", "6", "2", "3", "4", "5"), names(reactions))]

    for(i in 1:length(reactions)){
      reaction <- names(reactions)[i]
      nr.of.reactions <- reactions[i]
      cell.type <- ifelse(reaction %in% c(1,2,5,6), 1, 2)
      
      
      ## symmetric division
      if(reaction %in% c(1,3)){
        
        tree <- .division.s.p(cell.type=cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate, nr = nr.of.reactions)
        
      }else if(reaction==6){ ## asymmetric division of a stem cell

        tree <- .division.s.p(cell.type = cell.type, tree = tree, mut.rate.D = mut.rate.D, s.shape = s.shape, s.rate = s.rate, mutation.mode = mutation.mode, mut.rate = mut.rate, nr = nr.of.reactions, symmetric = F)
        
      }else if(reaction==2){ ## differentiation of stem cell into progenitor cell
        
        tree <- .differentiation.s.p(tree = tree, nr = nr.of.reactions)
        
      }else if(reaction==4){## loss of progenitor cell
        
        tree <- .loss(cell.type=cell.type, tree = tree, nr = nr.of.reactions)
        
      }else if(reaction==5){## loss of progenitor cell
        
        tree <- .loss(cell.type=cell.type, tree = tree, nr = nr.of.reactions)
        
      }
      
    }
    
    ## if drivers are introduced at a fixed time point rather than randomly: 
    if(driver.mode=="fixed_time" & time > t.driver & all(tree$sel.adv.==1)){
      ## randomly select the cell to acquire the driver mutant
      random.cell <- sample(seq(1,length(tree$tip.class[tree$tip.class==1])), 1)
      tree[tree$tip.class==1]$sel.adv.[random.cell] <- 1 + s.shape/s.rate
      
    }
    
    if(any(tree$sel.adv.>1)){
      max.clone <- max(table(tree$sel.adv.[tree$sel.adv.!=1]))
    }else{
      max.clone <- 0
    }
    
    n.cells <- c("1"=sum(tree$tip.class==1), "2"=sum(tree$tip.class==2))
  }
  
  return(res.list)
}

