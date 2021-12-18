library(Rcpp)

sourceCpp("Sling.cpp")

#' An implementation of Lasso regression via Coordinate Descent
#' @param X A size n*p matrix
#' @param y A size n vector
#' @param l A number
#' @param max_iter A integer number
#' @return A list containing w, itr
standardLassoRcpp = function(X, y, l, max_iter=1000){
  return(standardLasso(X, y, l, max_iter))
}

#' An implementation of Sling algorithm
#' @param X A size n*p matrix
#' @param y A size n vector
#' @param l A number
#' @param max_iter A integer number
#' @return A list containing w, itr
slingRcpp = function(X, y, l, max_iter=1000){
  return(sling(X, y, l, max_iter))
}