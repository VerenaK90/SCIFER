#' Subset a phylogenetic tree on a particular cell type
#' 
#' @param tree object of class `phylo`; the current tree. Needs the following additional list elements: `tip.class`, specifying the cell type and `sel.adv.`, specifying the selective advantage of each tip.
#' @param cell.type cell type of interest
#' @return the updated stat list
#' @export

.subset.cell.type <- function(tree, cell.type){
  
  if(!any(tree$tip.class==cell.type)){
    stop(paste("ERROR: no cell of type", cell.type, "in tree"))
  }
  if(all(tree$tip.class==cell.type)){
    stop(paste("ERROR: tree consists only of type", cell.type, "cells"))
  }
  ## the living cells are the tips
  dying.cells <- tree$tip.label[tree$tip.class!=cell.type]
  
  ## remove the dying cell from the tips
  new.tree <- tree
  new.tree$sel.adv. <- new.tree$sel.adv.[-which(new.tree$tip.label %in% dying.cells)]
  new.tree$tip.class <- new.tree$tip.class[-which(new.tree$tip.label %in% dying.cells)]
  
  new.tree <- ape::drop.tip(new.tree, dying.cells, trim.internal = T, collapse.singles = F)
  
  new.tree$tip.label <- setdiff(new.tree$edge[,2], new.tree$edge[,1])
  tree <- new.tree
  
  return(tree)
  
}

#' Sample a random subset of a phylogenetic tree
#' 
#' @param tree object of class `phylo`; the current tree. Needs the following additional list elements: `tip.class`, specifying the cell type and `sel.adv.`, specifying the selective advantage of each tip.
#' @param sample.size the number of cells to be sampled
#' @return the updated state list
#' @export

.simulate.sampling <- function(tree, sample.size){
  
  if(length(tree$tip.label) < sample.size){
    stop("ERROR: tree is smaller than given sample size")
  }
  
  ## the living cells are the tips
  living.cells <- tree$tip.label
  surviving.cells <- sample(living.cells, sample.size)
  dying.cells <- setdiff(living.cells, surviving.cells)
  
  ## remove the dying cell from the tips
  new.tree <- tree
  new.tree$sel.adv. <- new.tree$sel.adv.[-which(new.tree$tip.label %in% dying.cells)]
  new.tree$tip.class <- new.tree$tip.class[-which(new.tree$tip.label %in% dying.cells)]
  
  new.tree <- ape::drop.tip(new.tree, dying.cells, trim.internal = T, collapse.singles = F)
  
  new.tree$tip.label <- setdiff(new.tree$edge[,2], new.tree$edge[,1])
  tree <- new.tree
  
  return(tree)
}


#' Compute the VAF of a mutation in the population from a phylogenetic tree
#' 
#' @param tree object of class phylo; the current tree. Needs the following additional list elements: `tip.class`, specifying the cell type and `sel.adv.`, specifying the selective advantage of each tip.
#' @return a vector of VAFs
#' @export

get_vaf_from_tree <- function(tree){

  #tree <- Preorder(tree)
  tree <- ape::collapse.singles(tree)
  tree$tip.label <- sort(setdiff(tree$edge[,2], tree$edge[,1]))
  ## need to count for each mutation the number of cells it's present in
  vaf <- c()
  all.nodes <- unique(tree$edge[,2])

  for(i in all.nodes){
    parent <- tree$edge[tree$edge[,2]==i,1]
    mutations.this.node <- tree$edge.length[which(tree$edge[,2]==i)]
    children <- phytools::getDescendants(tree, i)
    tip.children <- intersect(tree$tip.label, children)
    n.children <- length(tip.children)
    vaf <- c(vaf, rep(n.children, mutations.this.node))
  }
  vaf <- vaf/length(tree$tip.label)/2
  vaf
}

#' Simulate the measured VAF after sequencing
#' 
#' @param vaf a vector of true VAFs in the population
#' @param depth average coverage in sequencing
#' @return a vector of simulated VAFs after sequencing
#' @export


simulate_vaf_upon_sequencing <- function(vaf, depth){
  sapply(vaf, function(x){
    cov <- rpois(1, depth)
    reads <- rbinom(1, size=cov, x)
    if(reads < 3){reads <- 0}
    reads/cov
  })
}

#' Compute the number of mutations per tip cell in a phylogenetic tree
#' 
#' @param tree object of class `phylo`; the current tree. Needs the following additional list elements: `tip.class`, specifying the cell type and `sel.adv.`, specifying the selective advantage of each tip.
#' @return a vector of mutation counts
#' @export

get_mutations_per_tip <- function(tree){
  
  ## need to count for each mutation the number of cells it's present in
  
  all.tips <- as.numeric(tree$tip.label)
  mutations.at.tips <- c()
  
  for(i in all.tips){
    ancestors.this.node <- phangorn::Ancestors(x = tree, node = i, "all")
    mutations.at.parents.and.node <- sum(tree$edge.length[which(tree$edge[,2] %in% c(i, ancestors.this.node))], na.rm=T)
    mutations.at.tips <- c(mutations.at.tips, mutations.at.parents.and.node)
  }
  
  mutations.at.tips
}



#' Get a binary mutation matrix from a phylogenetic tree
#' 
#' @param tree object of class `phylo`; the current tree. Needs the following additional list elements: `tip.class`, specifying the cell type and `sel.adv.`, specifying the selective advantage of each tip.
#' @return a matrix with rows corresponding to mutations and columns to cells; entries are binaries of 0 or 1, indicating, respectively, absence or presence of the mutation.
#' @export


get_matrix_from_tree <- function(tree){
  
  tree <- ape::collapse.singles(tree)
  tree$tip.label <- sort(setdiff(tree$edge[,2], tree$edge[,1]))
  
  ## need to assign to each mutation the cells it's present in
  
  mut.mat <- matrix(0, sum(tree$edge.length), length(tree$tip.label))
  
  all.nodes <- unique(tree$edge[,2])
  empty.muts <- which(rowSums(mut.mat)==0)
  
  for(i in all.nodes){
    mutations.this.node <- tree$edge.length[which(tree$edge[,2]==i)]
    children <- phytools::getDescendants(tree, i)
    tip.children <- intersect(tree$tip.label, children)

    mut.mat[empty.muts[1:mutations.this.node],tip.children] <- 1
    empty.muts <- empty.muts[-(1:mutations.this.node)]
  }
  
  mut.mat
}