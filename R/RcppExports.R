# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

vecSum <- function(vec) {
    .Call(`_Sling_vecSum`, vec)
}

Sum <- function(vec) {
    .Call(`_Sling_Sum`, vec)
}

Sum2 <- function(vec) {
    .Call(`_Sling_Sum2`, vec)
}

Sum3 <- function(vec) {
    .Call(`_Sling_Sum3`, vec)
}

Sum4 <- function(vec) {
    .Call(`_Sling_Sum4`, vec)
}

Sum5 <- function(vec) {
    .Call(`_Sling_Sum5`, vec)
}

Inner_Product1 <- function(vec1, vec2) {
    .Call(`_Sling_Inner_Product1`, vec1, vec2)
}

Inner_Product2 <- function(vec1, vec2) {
    .Call(`_Sling_Inner_Product2`, vec1, vec2)
}

mmult <- function(m, v) {
    .Call(`_Sling_mmult`, m, v)
}

mvmult <- function(m, v) {
    .Call(`_Sling_mvmult`, m, v)
}

standardLasso <- function(X, y, l, max_iter = 100L) {
    .Call(`_Sling_standardLasso`, X, y, l, max_iter)
}

sling <- function(X, y, l, max_iter = 100L) {
    .Call(`_Sling_sling`, X, y, l, max_iter)
}

