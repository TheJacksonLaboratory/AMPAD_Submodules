# Custom functions that are helpful over base R

lapply.n <- function(X, FUN) {
  
  tmp <- lapply(X, FUN)
  names(tmp) <- X
  return(tmp)
}
