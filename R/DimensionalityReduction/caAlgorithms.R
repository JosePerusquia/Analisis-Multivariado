############################################################
# Correspondence analysis                                        
# Author: Jose Antonio Perusquia Cortes
# Afil : Facultad de Ciencias - UNAM
# Module : Multivariate analysis
############################################################

############################################################
# Correspondence analysis using the GSVD decomposition on 
# the rows of a contingency table, for the columns the same
# algorithm should be used in the transposed table.
# The arguments are:
# x : contigency table in a matrix format
ca_gsvd = function(x){
  
  # Dimensions 
  dims = dim(x)
  r = dims[1]
  c = dims[2]
  
  # Total number of observations
  n = sum(x)
  
  # Sums by row and columns
  row_sum=rowSums(x)
  col_sum=colSums(x)
  
  # Masses and centroid
  row_masses=row_sum/n
  centroid=col_sum/n
  
  # Diagonal matrices with masses and centroid
  Dr=diag(row_masses)
  Dc=diag(centroid)
  
  # Square root of these diagonal matrices
  Dr_sq=sqrtm(Dr)
  Dc_sq=solve(sqrtm(Dc))
  
  # Profiles by row
  profiles = x / row_sum
  
  # Take out the centroid
  R = sweep(profiles, 2, centroid)
  
  # GSVD decomposition
  res=svd(Dr_sq%*%R%*%Dc_sq)
  
  # Eigenvalues
  values=res$d;values
  
  # Resulting matrices
  U=res$u;U
  V=res$v;V
  N=solve(Dr_sq)%*%U;N
  M=solve(Dc_sq)%*%V;M
  
  # Coordinates of the rows
  f=N[,c(1:2)]%*%diag(values[c(1:2)])
  
  return(f)
}
############################################################

############################################################
# Correspondane analysis using the SVD decomposition on the 
# correspondence matrix.
# The arguments are:
# x : contigency table in a matrix format
# The function returns:
# f = row coordinates
# g = col coordinates
# sv = singular values
ca_svd = function(x){
  
  # Dimensions 
  dims = dim(x)
  r = dims[1]
  c = dims[2]
  
  # Total number of observations
  n = sum(x)
  
  # Sums by row and columns
  row_sum=rowSums(x)
  col_sum=colSums(x)
  
  # Masses and centroid
  row_masses=row_sum/n
  centroid=col_sum/n
  
  # Diagonal matrices with masses and centroid
  Dr=diag(row_masses)
  Dc=diag(centroid)
  
  # Inverse of squared root
  Dr_sq=solve(sqrtm(Dr))
  Dc_sq=solve(sqrtm(Dc))
  
  # Transform the masses and centroid to matrix
  row_m=as.matrix(row_masses)
  centroid_m=as.matrix(centroid)
  
  # Correspondence matrix
  P=x/n
  
  # Standardised matrix
  A=Dr_sq%*%(P-(row_m%*%t(centroid_m)))%*%Dc_sq
  
  # SVD decomposition
  res=svd(A)
  
  # Eigenvalues and eigenvectors
  values=res$d
  U=res$u
  V=res$v
  
  # Standard coordinates
  X=Dr_sq%*%U
  Y=Dc_sq%*%V
  
  # Principal coordinates
  f=X%*%diag(res$d)
  g=Y%*%diag(res$d)
  
  l = list(f=f,g=g,sv=values)
  return(l)
}
############################################################

