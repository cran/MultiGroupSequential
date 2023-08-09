#' Update graph
#'
#' `updategraph()` updates the graph when only a subset of original hypotheses
#' is concerned.
#'
#' @param S1 Integer indices of the subset of hypotheses, S1 must be a non-empty
#'   subset of `S0` and must be sorted increasingly.
#' @param W0 Numeric vector for the initial weights of the graph.
#' @param G0 Numeric matrix of dimesion `length(W0)` by `length(W0)` for the
#'   initial transition matrix of the graph.
#' @param S0 Integer indices for the set of hypotheses from 1 to length of `W0`.
#' @return List with the following elements
#'   * `S1`: Integer indices the same as the input `S1`.
#'   * `W1`: Numeric vector for weights of the updated graph.
#'   * `G1`: Numeric transition of the updated graph.
#' @author Xiaodong Luo
#' @concept graphical procedure
#' @concept multiple comparison
#' @concept closed testing
#' @examples
#' ## We can use the function to produce a closed testing tree
#' ## A function to create power set
#' powerset <- function(x) {
#'   sets <- lapply(1:(length(x)), function(i) combn(x, i, simplify = FALSE))
#'   unlist(sets, recursive = FALSE)
#' }
#'
#' n <- 3    # number of hypotheses
#' pn <- 2^n-1
#' pset <- powerset(seq(1, n, by = 1))    # create the power set
#' df <- data.frame(matrix(ncol = 1+n, nrow = 0))    # create the dataset
#' colnames(df) <- c("Test", paste0("H", seq(1, n, by = 1), sep = ""))
#'
#' W0 <- c(1/3, 1/3, 1/3)    # the weights of the graph
#' m <- rbind(H1 = c(0, 1/2, 1/2),
#'            H2 = c(1/2, 0, 1/2),
#'            H3 = c(1/2, 1/2, 0))
#' G0 <- matrix(m, nrow = 3, ncol = 3)    # the transition matrix of the graph
#'
#' for (j in 1:pn){
#'     abc <- updategraph(S1 = pset[[j]], W0 = W0, G0 = G0)
#'     temp <- rep("-", n)
#'     temp[pset[[j]]] <- abc$W1
#'     temp <- c(paste(pset[[j]], collapse = ""), temp)
#'     df[j, ] <- temp
#' }
#' df    # the dataframe lists the closed testing tree
#' @export
updategraph=function (S1 = c(2,3),W0 = c(0.5,0.5,0,0),G0 = rbind(c(0,0,1,0),c(0,0,0,1),c(0,1,0,0),c(1,0,0,0)),S0 = seq(1,length(W0), by = 1)){
  SS = setdiff(S0, S1)
  nss = length(SS)
  WT = W0
  GT = G0
  if (nss > 0) {
    for (j in 1:nss) {
      sj = SS[j] + 1 - j  # this ensures we select the correct index after the graph is updated
      WT = WT + WT[sj] * GT[sj, ] # update W
      GT.new = GT
      for (l in 1:nrow(GT.new)) {  # update G
        temp=(GT[l,sj] * GT[sj,l])
        for (k in 1:ncol(GT.new)) {
          if (l==k | temp>=1){GT.new[l,k]=0}
          else {GT.new[l, k]=(GT[l, k] + GT[l, sj]*GT[sj, k])/(1 - temp)}
        }
      }
      WT <- WT[-sj]
      GT <- as.matrix(GT.new[-sj, -sj])
    }
  }
  list(S1 = S1, W1 = WT, G1 = GT)
}
