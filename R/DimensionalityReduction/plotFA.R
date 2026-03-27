################################################################################
# Plot function for factor analysis models                                       
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Multivariate Analysis
################################################################################

################################################################################
# plotFA uses qgraph library to plot the directed graph with the factors
# and their relationships to the original variables. The function receives
# an fa() model, the threshold to decide which edges to plot (by default it
# plots all of them), a boolean telling the function if an oblique rotation
# has been used to plot the correlation among factors and the names to appear
# in the variable nodes. Aesthetic variables have been hard-coded.
plotFA=function(fa_mod,threshold=0,oblique=F,nodesNames=NULL){
  
  # Obtain the dimensions of the model and the loadings  
  L = as.matrix(fa_mod$loadings)
  p = nrow(L)
  k = ncol(L)
  
  
  # Assign names to the variable nodes
  if(!is.null(nodesNames)){
    if(length(nodesNames)!=p){
      stop("length of nodesNames differs from number of variables")
    }else{
      rownames(L)=nodesNames
    }
  }
  
  # Creates the labels for the factor nodes
  factorNames = paste0("F", 1:k)
  colnames(L) = factorNames
  
  # Create adjacency matrix
  adj = matrix(0, nrow = p + k, ncol = p + k)
  adj[(p+1):(p+k), 1:p] = t(L)
  
  if(oblique){
    Phi = fa_mod$Phi
    diag(Phi)=0
    Phi[lower.tri(Phi)] = 0
    adj[(p+1):(p+k), (p+1):(p+k)] = Phi
    
    # Matrix indicating which edges are directed
    directed_mat = matrix(FALSE, nrow = p + k, ncol = p + k)
    
    # Only loadings are directed
    directed_mat[(p+1):(p+k), 1:p] = T
  }else{
    directed_mat =T
  }
  
  # Keep the edges with weights > threshold
  adj[abs(adj) <= threshold] = 0
  
  # Labels
  labels = c(rownames(L), colnames(L))
  
  # Node colors
  node_colors = c(rep("lightblue", p),rep("tomato", k))
  
  # Edge colors
  edge_colors = ifelse(adj > 0, "darkgreen", "darkred")
  if(oblique){
    edge_colors[(p+1):(p+k), (p+1):(p+k)] <- "darkblue"
  }
  
  # Plot
  qgraph(adj,
         layout = "spring",
         labels = labels,
         color = node_colors,
         esize = 1,
         edge.labels = round(adj, 2),
         edge.label.cex = 1.25,
         edge.color = edge_colors,
         edge.label.color='black',
         directed = directed_mat,
         arrows = TRUE)
}
################################################################################