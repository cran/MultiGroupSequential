##' Hommel procedure
##'
##' `hommeld()` implement the Hommel procedure with different alphas for
##' different endpoints.
##'
##' @details The package [hommel](https://cran.r-project.org/package=hommel)
##'   can handle Hommel procedure with different alpha's for different endpoints,
##'   the function `hommeld()` is just a wrapper of [hommel::hommel()].
##'
##' @inherit hochbergd
##' @author Xiaodong Luo
##' @concept Hommel procedure
##' @examples
##' hommeld(
##'   pvalues = runif(5),
##'   alpha = seq(0.01, 0.025, len = 5),
##'   epsilon = 1.0e-10,
##'   precision = 10
##' )
##' @export
hommeld=function(pvalues,alpha,epsilon=1.0e-10,precision=10){
  ## Hommel procedure that can handle different alphas for different endpoints
  alpha1=pmax(alpha,epsilon)
  qvalues=pvalues/alpha1
  #qvalues=round(qvalues,digits=precision)
  hom2 <- hommel::hommel(qvalues, simes = TRUE)
  decisions=(hommel::p.adjust(hom2) <= 1)
  list(decisions=decisions)
}
