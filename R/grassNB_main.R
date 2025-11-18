#'  @title Main parent function for GRASS-NB
#' @description Bayesian variable selection in spatially indexed count data with flexible choices of priors
#'
#' @param X is the matrix of features
#' @param y is the vector of the count outcome
#' @param group_ind is the vector of feature grouping indicators
#' @param NeighborhoodList is an object of class \code{nb}, typically created with \code{\link[spdep]{poly2nb}}. Each element of the list corresponds to a spatial unit and contains the indices of its neighboring units, as defined by the chosen contiguity rile (e.g queen or rook). Used to encode the spatial adjacency structure for analysis.
#' @param which.prior denotes the model to be used, "HS" or "SS", horseshoe and spike-and-slab respectively
#' @param niter is the number of MCMC iterations
#' @param verbose is TRUE or FALSE, to show or suppress progression updates
#' @param pop_col is a vector of the population of each spatial unit for the model offset, if modeling with space
#' @param scale is an integer of the desired scale for outcome rates (e.g 100000 for per 100,000 rates)
#' @param pi_star is the prior probability of there being any active feature in a group (when grouping features and using spike and slab)
#' @param alpha0 is the concentration parameter for the Beta hyper parameter for variable inclusion probability (when grouping features and using spike and slab). Lower values allow for more uncertainty and flexibility.

#' @importFrom stats var rbinom rgamma rbeta runif

#' @return It returns the estimated beta's, r, shrinkage parameters, phi, Delta and pDelta (if applicable )
#' @export


grassNB <- function(X, y, group_ind = NULL, NeighborhoodList = NULL, which.prior = c("HS", "SS"), niter=10000,
                    verbose = "TRUE", pop_col = NULL, scale = NULL, pi_star = 0.2, alpha0 = 10){

  # --- normalize ---
  y <- if (is.matrix(y)) as.numeric(y[,1]) else as.numeric(y)
  X <- as.matrix(X)
  if (!is.null(pop_col)) pop_col <- as.numeric(pop_col)

  N <- nrow(X); P <- ncol(X)
  which.prior <- match.arg(which.prior)
  if (!is.numeric(niter) || length(niter) != 1L || !is.finite(niter) || niter < 1 || niter %% 1 != 0) {
    stop("'niter' must be a single positive integer.")
  }

  # rate arguments must travel together
  if (xor(is.null(pop_col), is.null(scale))) {
    stop("Both 'pop_col' and 'scale' must be provided together when modeling rates, or both should be NULL when modeling counts.")
  }
  if (!is.null(scale) && (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) ||
                          scale <= 0 || scale %% 1 != 0)) {
    stop("'scale' must be a single positive integer (e.g., 100000).")
  }
  if (!is.null(pop_col)) {
    if (length(pop_col) != N) stop("Length of 'pop_col' must equal the number of rows of X (and length of y).")
    if (any(!is.finite(pop_col)) || any(pop_col <= 0)) stop("'pop_col' must be positive and finite for all units.")
  }

  # dimension checks (before any cbind/intercept)
  if (length(y) != N) stop("Length of y and number of rows of X do not match.")
  if (!is.null(NeighborhoodList) && length(NeighborhoodList) != N) stop("Length of y and length of NeighborhoodList do not match.")
  if (!is.null(group_ind) && length(group_ind) != ncol(X)) stop("Length of group_ind and number of columns of X do not match.")

  # class check for NeighborhoodList
  if (!is.null(NeighborhoodList) & !inherits(NeighborhoodList, "nb")) {
    stop("'NeighborhoodList' must be an object of class 'nb' (e.g., from spdep::poly2nb).")
  }

  # optional count-data warning
  if (isTRUE(verbose)) {
    if (any(!is.finite(y)) || any(y < 0) || any(y != floor(y))) {
      warning("Outcome 'y' should be nonnegative integers for counts; found non-integer or negative values.")
    } else if (length(unique(y)) > length(y) * 0.75) {
      warning("The data is possibly continuous; please double-check that the outcome is a discrete count.")
    }
  }


  #creating the design matrix for child functions
  K <- cbind(rep(1,N), X) #design matrix


  if (is.null(group_ind) & is.null(pop_col)) { # standard update w/o offset
    fitted_model <- VS_Standard(
      K = K, y = y, NeighborhoodList = NeighborhoodList,
      which.prior = which.prior,
      niter = niter, verbose = verbose,
      pop_col = pop_col, scale = scale
    )
  } else if (is.null(group_ind) & !is.null(pop_col)) { # standard update w/ offset
    fitted_model <- VS_Standard_offset(
      K = K, y = y, NeighborhoodList = NeighborhoodList,
      which.prior = which.prior,
      niter = niter, verbose = verbose,
      pop_col = pop_col, scale = scale
    )
  } else if (!is.null(group_ind) & is.null(pop_col)) { # group update w/o offset
    fitted_model <- VS_Group(
      K = K, y = y, group_ind = group_ind,
      NeighborhoodList = NeighborhoodList,
      which.prior = which.prior,
      niter = niter, verbose = verbose,
      pop_col = pop_col, scale = scale,
      pi_star = pi_star, alpha0 = alpha0
    )
  }else { # group update w/ offset
    fitted_model <- VS_Group_offset(
      K = K, y = y, group_ind = group_ind,
      NeighborhoodList = NeighborhoodList,
      which.prior = which.prior,
      niter = niter, verbose = verbose,
      pop_col = pop_col, scale = scale,
      pi_star = pi_star, alpha0 = alpha0
    )
  }

  return(fitted_model)

}




