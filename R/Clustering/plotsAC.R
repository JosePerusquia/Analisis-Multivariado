##################################################################
# Cluster analysis auxiliary functions                           
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Multivariate Analysis
##################################################################

##################################################################
# Plot dendrogram function
# Args:
#   - hc_obj : 
#   - title  : title of the plot
plot_dendrogram = function(hc_obj, title = "") {
  dend = as.dendrogram(hc_obj)
  dend_data = dendro_data(dend, type = "rectangle")
  
  p= ggplot(dend_data$segments) + 
     geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
     geom_text(data = dend_data$labels, 
              aes(x, y, label = label),
              hjust = 1, angle = 90, size = 2) +
     ylim(-0.15, max(hc_obj$height)) +
     theme_minimal() +
     labs(x = "", y = "Height", title = title)
  
  print(p)
  return(p)
}
##################################################################

##################################################################
# Plots the silhouette and returns potential misclassifications
# Args: 
#   - clusters 
#   - dist_mat : distance matrix used
#   - title
#   - col      : vector of colors to be used (one for each cluster)
plot_silhouette = function(clusters, dist_mat,
                           title = "Silhouette", col){
  
  si = silhouette(clusters, dist_mat)
  plot(si, col = col , main = title)
  
  indices = which(si[,3] < 0)
  
  if(length(indices) == 0){
    
    misclass = NULL
    
  } else {
    
    misclass = data.frame(
      index     = indices,
      cluster   = si[indices,1],
      neighbor  = si[indices,2],
      sil_width = si[indices,3]
    )
  }
  
  return(list(
    silhouette = si,
    MC = misclass
  ))
}
##################################################################

##################################################################
# Projects the centroids or medoids into the pca plot
# Args:
#   - points  : points to be projected
#   - pca_obj : a prcomp object
project_to_pca = function(points, pca_obj) {
  (points - pca_obj$center) %*% pca_obj$rotation[, 1:2]
}
##################################################################
