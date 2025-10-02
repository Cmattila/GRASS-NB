#' @title AR(1) correlation matrix
#' @description
#' Construct an \eqn{n \times n} AR(1) correlation matrix with parameter \code{rho}.
#'
#' @param n Integer \eqn{\ge 1}. Matrix dimension.
#' @param rho Numeric in \eqn{(-1, 1)}. AR(1) correlation parameter.
#'
#' @return A numeric \eqn{n \times n} correlation matrix with entries \eqn{\rho^{|i-j|}}.
#' @export
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho^exponent
}


#' @title Generate indices of non-null predictors under grouped structures
#'
#' @description This function selects indices of non-null predictors based on group membership.
#' Two group structures are supported:
#' \describe{
#'   \item{\code{"Case1"}}{Selects \strong{two entire groups at random} whose sizes sum exactly to \code{num_nonnull}.
#'   All variables from those two groups are marked as non-null.}
#'   \item{\code{"Case2"}}{Distributes \code{num_nonnull} non-nulls across randomly chosen groups,
#'   assigning approximately \code{per_group_nonnull} non-nulls per group (with remainder distributed randomly).}
#' }
#'
#' @param p Integer. Total number of predictors including the intercept.
#'   Since the intercept is excluded, there are \code{p - 1} variables available for selection.
#' @param group_ind Integer or factor vector of length \code{p - 1}.
#'   Group membership indicator for each predictor (excluding the intercept).
#' @param num_nonnull Integer. Total number of predictors to set as non-null.
#' @param nonnull_group_structure Character string, either \code{"Case1"} or \code{"Case2"}.
#'   Controls how non-null variables are selected across groups.
#' @param per_group_nonnull Integer. (Case2 only) Number of non-nulls to allocate per group on average.
#'
#' @return An integer vector of indices (1-based, matching \code{beta_true} indexing) of selected non-null predictors.
#'
#' @examples
#' set.seed(123)
#' p <- 11
#' group_ind <- rep(1:5, each = 2)  # 10 predictors, 5 groups of size 2
#'
#' # Case1: choose two entire groups summing to num_nonnull
#' generate_nonnull_indices(p, group_ind, num_nonnull = 4, nonnull_group_structure = "Case1")
#'
#' # Case2: assign ~2 per group
#' generate_nonnull_indices(p, group_ind, num_nonnull = 6, nonnull_group_structure = "Case2", per_group_nonnull = 2)
#'
#' @export


generate_nonnull_indices <- function(
    p, group_ind, num_nonnull,
    nonnull_group_structure = c("Case1", "Case2"),
    per_group_nonnull = 3
){
  nonnull_group_structure <- match.arg(nonnull_group_structure)

  if ((p - 1) < num_nonnull) stop("num_nonnull must be <= p - 1")
  if (length(group_ind) != (p - 1)) stop("group_ind must be of length p - 1")

  #  group ids as character keys
  group_ids_chr <- as.character(sort(unique(group_ind)))
  G  <- length(group_ids_chr)
  Mg <- table(as.character(group_ind))              # named by group id
  group_indices <- split(1:(p - 1), as.character(group_ind))

  sample_nonNULL <- integer(0)

  if (nonnull_group_structure == "Case1") {
    if (G < 2) stop("Need at least two groups for Case1.")

    # all unordered pairs of groups
    all_pairs <- utils::combn(group_ids_chr, 2, simplify = FALSE)

    # keep pairs whose total size equals num_nonnull
    matching_pairs <- Filter(function(pr) sum(Mg[pr]) == num_nonnull, all_pairs)

    if (!length(matching_pairs)) {
      # diagnostics
      pair_totals <- vapply(all_pairs, function(pr) sum(Mg[pr]), numeric(1))
      stop(
        paste0(
          "No pair of full groups sums to num_nonnull = ", num_nonnull, ". ",
          "Available pair totals: ", paste(sort(unique(pair_totals)), collapse = ", "), "."
        )
      )
    }

    # randomly choose one pair and take ALL indices from those two groups
    chosen_pair <- matching_pairs[[sample.int(length(matching_pairs), 1)]]
    sample_nonNULL <- unlist(group_indices[chosen_pair], use.names = FALSE)

    # sanity check
    if (length(sample_nonNULL) != num_nonnull) {
      stop("Internal error: selected pair does not match num_nonnull after selection.")
    }

  } else if (nonnull_group_structure == "Case2") {
    num_groups_needed <- floor(num_nonnull / per_group_nonnull)
    remainder <- num_nonnull %% per_group_nonnull

    if (num_groups_needed > G) stop("Not enough groups to allocate nonzero betas.")

    # sample by group ids
    selected_groups <- sample(group_ids_chr, num_groups_needed)
    extra_groups <- if (remainder > 0) sample(selected_groups, remainder) else character(0)

    for (gid in selected_groups) {
      group_vars <- group_indices[[gid]]
      k <- per_group_nonnull + as.integer(gid %in% extra_groups)

      if (length(group_vars) < k) {
        stop(paste("Group", gid, "does not have enough variables for", k, "nonzeros."))
      }
      sample_nonNULL <- c(sample_nonNULL, sample(group_vars, k))
    }
  }

  # shift by +1 to match beta_true indexing convention
  nonnull_locs <- sort(sample_nonNULL + 1)
  return(nonnull_locs)
}

#' @title Simulate grouped features with AR(1) within-group correlation
#' @description
#' Generates block-diagonal Gaussian features by group using an AR(1) correlation
#' within each group. Returns design matrices \code{X} and \code{K}, plus a
#' \code{beta_true} vector with signals at \code{nonNull_indices}.
#'
#' @param N Sample size.
#' @param p Total number of coefficients including intercept.
#' @param group_ind Grouping vector for predictors (length \code{p-1}).
#' @param AR AR(1) correlation parameter within groups.
#' @param nonNull_indices Indices in \code{2:p} set to non-null.
#' @param signal_vec Signal values to assign at \code{nonNull_indices}.
#'
#' @importFrom mvnfast rmvn
#' @return A list containing \code{beta_true}, \code{X}, \code{K}, and related metadata.
#' @export
feature_sim_grouped <- function(N, p, group_ind, AR = 0.25,
                              nonNull_indices, signal_vec) {
  if (length(group_ind) != (p - 1)) stop("group_ind must be of length p - 1")
  if (any(nonNull_indices < 2 | nonNull_indices > p)) stop("nonNull_indices must be in 2:p")

  G <- length(unique(group_ind))
  Mg <- table(group_ind)

  group_indices <- split(1:(p-1), group_ind)

  # Simulate block-diagonal X matrix
  Xall <- vector("list", G)
  for (g in 1:G) {
    Sigma <- ar1_cor(Mg[g], AR)
    Xall[[g]] <- mvnfast::rmvn(N, rep(0, Mg[g]), 0.1 * Sigma)
  }
  X <- do.call(cbind, Xall)
  K <- cbind(1, X)  # Add intercept

  beta_true <- rep(0, p)
  beta_true[nonNull_indices] <- signal_vec  # assign signal values
  beta_true[1] <- 0                    # intercept remains 0

  sample_null <- setdiff(2:p, nonNull_indices)

  return(list(
    beta_true = beta_true,
    X = X,
    K = K,
    signal_values = signal_vec,
    N = N,
    p = p,
    sample_nonNULL = nonNull_indices,
    sample_null = sample_null,
    group_ind = group_ind
  ))
}

#' @title Generate negative-binomial outcomes
#' @description
#' Simulate negative-binomial outcomes with logit link. Can include an offset
#' with \code{pop_col} and \code{Scale}, plus spatial/random effects \code{phi_true}.
#'
#' @param N Sample size.
#' @param K Design matrix (with intercept in column 1).
#' @param r Negative-binomial dispersion parameter.
#' @param beta_true Coefficient vector.
#' @param Scale Optional scale for offset.
#' @param phi_true Spatial/random effect vector.
#' @param offset Logical or NULL; if non-NULL, include offset term.
#' @param pop_col Optional population offset vector.
#'
#' @return A list with elements \code{y} and \code{Scale}.
#' @export
y_gen_fun <- function(N, K, r=1, beta_true, Scale = NULL, phi_true, offset = NULL, pop_col = NULL){
  nis <- rep(1,N)
  Scale <- Scale
  if(is.null(offset)==FALSE){
    eta <- K%*%c(beta_true)+ log(pop_col) - log(Scale) +
      rep(phi_true,nis)} #see if wrap works
  else{eta <- K%*%c(beta_true)+ rep(phi_true,nis)}
  psi <- exp(eta)/(1+exp(eta)) # Prob of success
  q <- 1-psi                   # Note: rnbinom models q=1-p
  q[q < 1e-5] <- 1e-5
  y <- stats::rnbinom(N,r,q)
  return(list(y=y, Scale = Scale))
}

#' @title Create a grouping index
#' @description
#' Split \code{p - 1} non-intercept predictors into \code{G} groups, as evenly as possible.
#' If division is uneven, earlier groups receive one extra predictor.
#'
#' @param p Total number of coefficients including intercept.
#' @param G Number of groups.
#'
#' @return Integer vector of group indices of length \code{p - 1}.
#' @export
generate_group_ind <- function(p, G) {
  if ((p - 1) %% G != 0) {
    warning("p - 1 is not divisible by G. Groups will be uneven.")
  }
  split_sizes <- rep(floor((p - 1) / G), G)
  remainder <- (p - 1) %% G
  if (remainder > 0) split_sizes[1:remainder] <- split_sizes[1:remainder] + 1
  group_ind <- rep(1:G, times = split_sizes)
  return(group_ind)
}

#' @title Sample spatial random effects
#' @description
#' Draw spatial random effects \eqn{\phi} from a canonical MVN with precision
#' \eqn{\Lambda = \nu (Q + \epsilon I)} for each value of \code{nu}. Enforces a
#' sum-to-zero constraint.
#'
#' @param nu_vals Vector of precision parameters.
#' @param Q Precision matrix (e.g., ICAR).
#' @param n_space Dimension of the spatial domain.
#'
#' @importFrom spam rmvnorm.canonical
#' @return A list of centered \code{phi} vectors, one per \code{nu}.
#' @export
phi_True_func <- function(nu_vals, Q, n_space) {
  result <- list()

  for (nu in nu_vals) {
    # Construct the precision matrix
    Lambda <- nu * (Q + diag(10^{-8},n_space))

    # Sample one phi ~ MVN_canonical(0, Lambda)
    phi <- rmvnorm.canonical(1, b = rep(0, n_space), Q = as.matrix(Lambda))
    phi <- spam::rmvnorm.canonical(1, b = rep(0, n_space), Q = as.matrix(Lambda))
    phi <- as.numeric(phi - mean(phi))  # Sum-to-zero constraint

    # Store with informative name
    result[[paste0("phi_nu_", nu)]] <- phi
    cat("nu =", nu, "phi range =", range(phi), "actual SD =", stats::sd(phi), "\n")
  }

  return(result)
}


