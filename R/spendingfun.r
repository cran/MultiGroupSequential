#' Calculate alpha spending function
#'
#' `spendingfun()` calculates the alpha spending function.
#'
#' @details
#'   * `"OBF"`: O'Brien-Fleming family; \eqn{2\{1-\Phi(\Phi^{-1}(1-\alpha/2)/t^{\rho/2})\}};
#'   * `"pocock"`: Pocock family; \eqn{\alpha \log\{1+(e-1)*t\}};
#'   * `"power"`: Power family; \eqn{\alpha*t^{\rho}}
#'
#'   Note that the OBF and Pocock spending functions are not the originally
#'   proposed ones, they are the modified ones that closely resemble the original
#'   versions. That being said, you might still see some differences.
#' @param alpha Numeric scalar of the overall alpha to be spent.
#' @param fractions Numeric vector of the sequence of  information fractions.
#'   All elements should be between 0 and 1 with the last one being exactly 1.
#' @param family Character scalar for the family of spending functions, one of
#'   `"OBF"`, `"pocock"`, `"power"`.
#' @param rho Numeric scalar of auxiliary parameter for O'Brien-Fleming and
#'   power family.
#' @return List with an element named `aseq` for the alpha spending sequence.
#' @author Xiaodong Luo
#' @concept group-sequential
#' @concept alpha-spending
#' @examples
#' spendingfun(
#'   alpha = 0.025,
#'   fractions = seq(0.2, 1, by = 0.2),
#'   family = "OBF",
#'   rho = 1
#' )
#' @export
spendingfun=function(alpha,fractions=seq(0.2,1,by=0.2),family="OBF",rho=1){
  if (family=='OBF'){
    qa=qnorm(1-alpha/2)/fractions^(rho/2)
    aseq=2*(1-pnorm(qa))
  }
  else if (family=='pocock'){
    qa=1+(exp(1)-1)*fractions
    aseq=alpha*log(qa)
  }
  else if (family=='power'){
    aseq=alpha*fractions^rho
  }
  list(aseq=aseq)
}
