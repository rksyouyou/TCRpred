
psd <- function(S,M){
  if (isSymmetric(S)==FALSE){
    stop("S should be symmetric.")
  }
  es = eigen(S)
  S <- as.matrix(S)
  Q <- es$values
  m <- sum(Q>=0)
  P <- es$vectors
  n <- dim(S)[1]
  my.psd <- matrix(0, n, n)
  rownames(my.psd) <- rownames(S)
  colnames(my.psd) <- colnames(S)
  for (j in 1: min(m,M)){
    my.psd <- my.psd + Q[j]*P[,j]%*%t(P[,j])
  }
  return(my.psd)
}
