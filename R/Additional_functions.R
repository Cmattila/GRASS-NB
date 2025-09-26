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


#' @title Generate indices of non-null coefficients under grouped structures
#' @description
#' Selects positions (in \code{2:p}) to be non-zero according to a grouped scheme.
#' In \code{Case1}, picks a combination of whole groups whose sizes sum exactly to
#' \code{num_nonnull}. In \code{Case2}, spreads \code{num_nonnull} across randomly
#' chosen groups with \code{per_group_nonnull} per group (plus any remainder).
#'
#' @param p Total number of coefficients including intercept at index 1.
#' @param group_ind Integer vector of length \code{p - 1} giving group id (1..G) for each predictor.
#' @param num_nonnull Total number of non-null coefficients to allocate.
#' @param nonnull_group_structure One of \code{c("Case1","Case2")}. See details.
#' @param per_group_nonnull Target number of non-nulls per group in \code{Case2}.
#'
#' @return Integer vector of sorted non-null indices in \code{2:p}.
#' @export
generate_nonnull_indices <- function(p, group_ind, num_nonnull,
                                     nonnull_group_structure = c("Case1", "Case2"),
                                     per_group_nonnull = 3) {
  nonnull_group_structure <- match.arg(nonnull_group_structure)

  if ((p - 1) < num_nonnull) stop("num_nonnull must be <= p - 1")
  if (length(group_ind) != (p - 1)) stop("group_ind must be of length p - 1")

  G <- length(unique(group_ind))
  Mg <- table(group_ind)
  group_indices <- split(1:(p - 1), group_ind)
  sample_nonNULL <- c()

  if (nonnull_group_structure == "Case1") {
    found <- FALSE
    group_ids <- as.character(sort(unique(group_ind)))

    for (k in 1:length(group_ids)) {
      group_combos <- combn(group_ids, k, simplify = FALSE)
      for (grp_set in group_combos) {
        total_size <- sum(Mg[grp_set])
        if (total_size == num_nonnull) {
          sample_nonNULL <- unlist(group_indices[grp_set])
          found <- TRUE
          break
        }
      }
      if (found) break
    }

    if (!found) {
      stop("No combination of groups adds up exactly to num_nonnull.")
    }

  } else if (nonnull_group_structure == "Case2") {
    num_groups_needed <- floor(num_nonnull / per_group_nonnull)
    remainder <- num_nonnull %% per_group_nonnull

    if (num_groups_needed > G) stop("Not enough groups to allocate nonzero betas.")

    selected_groups <- sample(1:G, num_groups_needed)
    extra_groups <- if (remainder > 0) sample(selected_groups, remainder) else integer(0)

    for (g in selected_groups) {
      group_vars <- group_indices[[as.character(g)]]

      k <- per_group_nonnull + as.integer(g %in% extra_groups)  # add extra if this group is selected for remainder

      if (length(group_vars) < k) {
        stop(paste("Group", g, "does not have enough variables for", k, "nonzeros."))
      }

      sample_nonNULL <- c(sample_nonNULL, sample(group_vars, k))
    }
  }

  nonnull_locs <- sort(sample_nonNULL + 1)  # shift to match beta_true
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
  y <- rnbinom(N,r,q)
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
#' @return A list of centered \code{phi} vectors, one per \code{nu}.
#' @export
phi_True_func <- function(nu_vals, Q, n_space) {
  result <- list()

  for (nu in nu_vals) {
    # Construct the precision matrix
    Lambda <- nu * (Q + diag(10^{-8},n_space))

    # Sample one phi ~ MVN_canonical(0, Lambda)
    phi <- rmvnorm.canonical(1, b = rep(0, n_space), Q = as.matrix(Lambda))
    phi <- as.numeric(phi - mean(phi))  # Sum-to-zero constraint

    # Store with informative name
    result[[paste0("phi_nu_", nu)]] <- phi
    cat("nu =", nu, "phi range =", range(phi), "actual SD =", sd(phi), "\n")
  }

  return(result)
}


