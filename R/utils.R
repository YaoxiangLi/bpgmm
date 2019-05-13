#' Tool for vector to matrix
#'
#' @param ZOneDim a vector.
#' @param m the number of cluster.
#' @param n sample size.
#' @return adjacency matrix
#'
#' @export
getZmat <- function(ZOneDim, m, n) {
  Zmat <- matrix(NA, m, n)
  for (i in 1:m) {
    Zmat[i,] = as.numeric(ZOneDim == i)
  }
  Zmat
}

#' Log scale ratio calculation
#'
#' @param deno denominator.
#' @param nume numerator.
#' @return result of ratio
#'
#' @export
calculateRatio <- function(deno, nume) {
  ## deno nume both in log scale
  maxNume = max(nume)
  transDeno = deno - maxNume
  transNume = nume - maxNume
  res = exp(transDeno) / (sum(exp(transNume)))
  res
}


#' Convert list of string to vector of string
#'
#' @param stringList list of string
#' @return vector of string
#'
#' @export
listToStrVec <- function(stringList) {
  for (i in 1:length(stringList)) {
    stringList[[i]] <-
      paste0("(", paste0(stringList[[i]], collapse = ""), ")")
  }
  unlist(stringList)
}

