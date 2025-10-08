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


#' @title Generate indices of non-null features for simulation with singleton consideration.
#'
#' @description This function selects indices of non-null features based on group structure.
#' If `singletons = TRUE`, half of the non-nulls are selected from grouped features
#' and the other half from singleton groups (groups of size 1).
#'
#' @param p Integer. Total number of features (excluding intercept).
#' @param group_ind Integer vector of length `p - 1` mapping each feature to its group.
#' @param num_nonnull Integer. Total number of non-null features to select.
#' @param nonnull_group_structure Character. Either `"Case1"` or `"Case2"`; determines how non-nulls are distributed across groups.
#' @param per_group_nonnull Integer. Number of non-nulls per group (used in `"Case2"`).
#' @param singletons Logical. If TRUE, half of the non-nulls are selected from singleton groups.
#'
#' @return Integer vector of indices (1-based) of non-null features.
#'
#' @examples
#' generate_nonnull_indices_singletons(p = 200, group_ind = rep(1:50, each = 2), num_nonnull = 20, singletons = TRUE)
#' @export
generate_nonnull_indices_singletons <- function(
    p, group_ind, num_nonnull,
    nonnull_group_structure = c("Case1", "Case2"),
    per_group_nonnull = 3,
    singletons = TRUE
) {
  nonnull_group_structure <- match.arg(nonnull_group_structure)

  if ((p - 1) < num_nonnull) stop("num_nonnull must be <= p - 1")
  if (length(group_ind) != (p - 1)) stop("group_ind must be of length p - 1")

  group_ids_chr <- as.character(sort(unique(group_ind)))
  Mg <- table(as.character(group_ind))
  group_indices <- split(1:(p - 1), as.character(group_ind))

  # Identify singleton vs grouped groups
  singleton_ids <- names(Mg)[Mg == 1]
  grouped_ids   <- setdiff(group_ids_chr, singleton_ids)

  sample_nonNULL <- integer(0)

  # ---- helper for Case1: find subset of full groups summing to target ----
  .find_subset_equal <- function(Mg_named, target, restrict_ids = names(Mg_named)) {
    if (target == 0L) return(list(character(0)))
    Mg_sub <- Mg_named[restrict_ids]
    gids <- names(Mg_sub); sizes <- as.integer(Mg_sub)
    sol <- list()
    recurse <- function(i, cur_ids, cur_sum) {
      if (cur_sum == target) {
        sol[[length(sol) + 1L]] <<- cur_ids; return()
      }
      if (cur_sum > target || i > length(gids)) return()
      recurse(i + 1L, c(cur_ids, gids[i]), cur_sum + sizes[i])
      recurse(i + 1L, cur_ids, cur_sum)
    }
    recurse(1L, character(0), 0L)
    sol
  }

  # ---- helper for Case2: distribute non-nulls with per_group_nonnull + remainders ----
  # ---- helper: Case2 allocation (base per_k per group; entire remainder into ONE distinct group) ----
  .alloc_case2 <- function(target, candidate_ids, per_k, Mg_named) {
    if (target == 0L) return(integer(0))
    if (length(candidate_ids) == 0L) stop("No candidate groups available for Case2.")

    # capacities
    cap <- as.integer(Mg_named[candidate_ids]); names(cap) <- candidate_ids

    # how many full groups (per_k each) and the remainder
    num_full <- target %/% per_k
    rem      <- target %%  per_k

    if (num_full > length(candidate_ids)) {
      stop("Not enough distinct groups to allocate base per_group_nonnull.")
    }

    # pick distinct groups for the base allocation
    selected_full <- if (num_full > 0L) sample(candidate_ids, num_full, replace = FALSE) else character(0)

    # base k per selected group
    k_by_group <- if (num_full > 0L) setNames(rep(per_k, num_full), selected_full) else integer(0)

    # check base capacity
    if (num_full > 0L && any(cap[selected_full] < per_k)) {
      bad <- names(cap[selected_full])[cap[selected_full] < per_k]
      stop(paste0("Group(s) ", paste(bad, collapse = ", "),
                  " lack capacity for per_group_nonnull = ", per_k, "."))
    }
    # consume capacity for base
    if (num_full > 0L) cap[selected_full] <- cap[selected_full] - per_k

    # put entire remainder into ONE group, distinct from base groups
    if (rem > 0L) {
      pool_new <- setdiff(candidate_ids, selected_full)
      if (length(pool_new) == 0L) {
        stop("No distinct group available to receive the remainder for Case2.")
      }
      # must have capacity>=rem
      pool_ok  <- pool_new[cap[pool_new] >= rem]
      if (length(pool_ok) == 0L) {
        # no single distinct group can accommodate the whole remainder -> be explicit
        stop(paste0(
          "No single distinct group has capacity ", rem, " for the remainder. ",
          "Remaining capacities among unused groups: ",
          paste(sprintf("%s:%d", pool_new, cap[pool_new]), collapse = ", "), "."
        ))
      }
      rem_group <- sample(pool_ok, 1L)
      # initialize entry if needed, then add remainder
      if (!(rem_group %in% names(k_by_group))) k_by_group[rem_group] <- 0L
      k_by_group[rem_group] <- k_by_group[rem_group] + rem
      cap[rem_group] <- cap[rem_group] - rem
    }

    # final sanity
    if (sum(k_by_group) != target) stop("Internal alloc error: sum(k_by_group) != target.")
    if (any(k_by_group > Mg_named[names(k_by_group)])) stop("Internal alloc error: over-capacity.")
    k_by_group
  }
  if (singletons) {
    num_grouped    <- floor(num_nonnull / 2)
    num_singletons <- num_nonnull - num_grouped

    if (nonnull_group_structure == "Case2") {
      k_by_group <- .alloc_case2(num_grouped, grouped_ids, per_group_nonnull, Mg)
      for (gid in names(k_by_group)) {
        k <- k_by_group[[gid]]
        group_vars <- group_indices[[gid]]
        if (length(group_vars) < k)
          stop(paste("Group", gid, "too small for", k, "nonzeros."))
        sample_nonNULL <- c(sample_nonNULL, sample(group_vars, k))
      }
    } else if (nonnull_group_structure == "Case1" && num_grouped > 0) {
      matches <- .find_subset_equal(Mg, num_grouped, grouped_ids)
      if (length(matches) == 0L) stop("No Case1 grouping possible for given target.")
      chosen <- matches[[sample.int(length(matches), 1L)]]
      sample_nonNULL <- unlist(group_indices[chosen], use.names = FALSE)
    }

    if (num_singletons > 0) {
      singleton_vars <- unlist(group_indices[singleton_ids], use.names = FALSE)
      if (length(singleton_vars) < num_singletons) stop("Not enough singleton vars.")
      sample_nonNULL <- c(sample_nonNULL, sample(singleton_vars, num_singletons))
    }

  } else {
    if (nonnull_group_structure == "Case2") {
      k_by_group <- .alloc_case2(num_nonnull, group_ids_chr, per_group_nonnull, Mg)
      for (gid in names(k_by_group)) {
        k <- k_by_group[[gid]]
        group_vars <- group_indices[[gid]]
        if (length(group_vars) < k)
          stop(paste("Group", gid, "too small for", k, "nonzeros."))
        sample_nonNULL <- c(sample_nonNULL, sample(group_vars, k))
      }
    } else if (nonnull_group_structure == "Case1") {
      matches <- .find_subset_equal(Mg, num_nonnull, group_ids_chr)
      if (length(matches) == 0L) stop("No Case1 grouping possible for given target.")
      chosen <- matches[[sample.int(length(matches), 1L)]]
      sample_nonNULL <- unlist(group_indices[chosen], use.names = FALSE)
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


#' @title Simulate grouped features with AR(1) within-group correlation with singleton consideration
#' @description
#' Generates block-diagonal Gaussian features by group using an AR(1) correlation
#' within each group. If `singletons = TRUE`, singleton features are simulated jointly
#' with a shared correlation structure.
#'
#' @param N Sample size.
#' @param p Total number of coefficients including intercept.
#' @param group_ind Grouping vector for predictors (length \code{p-1}).
#' @param AR AR(1) correlation parameter within groups and between singletons.
#' @param nonNull_indices Indices in \code{2:p} set to non-null.
#' @param signal_vec Signal values to assign at \code{nonNull_indices}.
#' @param singletons Logical. If TRUE, singleton features are simulated jointly with correlation.
#'
#' @importFrom mvnfast rmvn
#' @return A list containing \code{beta_true}, \code{X}, \code{K}, and related metadata.
#' @export

feature_sim_grouped_singletons <- function(
    N, p, group_ind, AR = 0.25,
    nonNull_indices, signal_vec,
    singletons = TRUE
) {
  # Input checks
  if (length(group_ind) != (p - 1)) {
    stop("group_ind must be of length p - 1")
  }
  if (any(nonNull_indices < 2 | nonNull_indices > p)) {
    stop("nonNull_indices must be in 2:p")
  }
  if (length(signal_vec) != length(nonNull_indices)) {
    stop("signal_vec must have the same length as nonNull_indices")
  }

  # Group bookkeeping
  group_ids     <- unique(group_ind)
  Mg            <- table(group_ind)                     # group sizes
  group_indices <- split(1:(p - 1), group_ind)          # original positions by group
  singleton_ids <- names(Mg)[Mg == 1]
  grouped_ids   <- setdiff(group_ids, singleton_ids)

  # We'll fill features back into their ORIGINAL positions
  Xall <- vector("list", p - 1)

  # Common variance scale: want marginal SD = 0.1 -> var = 0.01
  var_scale <- 0.01

  # Simulate grouped features (size >= 2)
  for (gid in grouped_ids) {
    idx <- group_indices[[gid]]        # original positions for this group
    m   <- length(idx)
    Sigma <- ar1_cor(m, AR) * var_scale
    Xblock <- mvnfast::rmvn(N, rep(0, m), Sigma)  # N x m

    # Place columns back into original positions (preserve order)
    # Convert to list of N x 1 matrices
    Xcols <- lapply(seq_len(m), function(j) matrix(Xblock[, j], ncol = 1))
    Xall[idx] <- Xcols
  }

  # Simulate singleton features
  if (length(singleton_ids) > 0) {
    # Figure out the original positions of singleton columns (numeric indices 1:(p-1))
    singleton_pos <- unlist(group_indices[singleton_ids], use.names = FALSE)

    if (singletons) {
      # Jointly correlated singletons: constant correlation "AR" off-diagonal
      m_singleton <- length(singleton_pos)
      # Build correlation matrix safely
      if (m_singleton == 1) {
        Sigma_singleton <- matrix(1 * var_scale, 1, 1)
        X_singletons <- matrix(rnorm(N, 0, sqrt(var_scale)), ncol = 1)
      } else {
        Sigma_singleton <- matrix(AR, m_singleton, m_singleton)
        diag(Sigma_singleton) <- 1
        Sigma_singleton <- Sigma_singleton * var_scale
        X_singletons <- mvnfast::rmvn(N, rep(0, m_singleton), Sigma_singleton)  # N x m_singleton
      }
      # Put each column back to its original place as N x 1 matrices
      for (j in seq_along(singleton_pos)) {
        Xall[[singleton_pos[j]]] <- matrix(X_singletons[, j], ncol = 1)
      }
    } else {
      # Independent singletons
      for (pos in singleton_pos) {
        Xall[[pos]] <- matrix(rnorm(N, mean = 0, sd = sqrt(var_scale)), ncol = 1)
      }
    }
  }

  # Assemble X and K
  X <- do.call(cbind, Xall)   # N x (p-1), in ORIGINAL order
  K <- cbind(1, X)            # add intercept as column 1  -> N x p

  # True coefficients
  beta_true <- rep(0, p)
  beta_true[nonNull_indices] <- signal_vec
  beta_true[1] <- 0  # intercept remains 0

  # Null bookkeeping
  sample_null <- setdiff(2:p, nonNull_indices)

  # Return
  return(list(
    beta_true      = beta_true,
    X              = X,
    K              = K,
    signal_values  = signal_vec,
    N              = N,
    p              = p,
    sample_nonNULL = nonNull_indices,
    sample_null    = sample_null,
    group_ind      = group_ind
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


#' @title Assign feature groups for structured variable selection  with singleton consideration.
#'
#' @description This function assigns features to groups based on a specified vector of group sizes (`Mg`).
#' If `singletons = TRUE`, any remaining features not covered by `Mg` are assigned to singleton groups.
#'
#' @param p Integer. Total number of features (including intercept).
#' @param G Number of groups (assumes equal group sizes).
#' @param singletons Logical. If TRUE, (p/2)-1 remaining features are assigned to singleton groups.
#'
#' @return A list containing:
#'   \item{group_id}{Integer vector of length `p - 1` mapping each feature to its group.}
#'
#' @examples
#' generate_group_ind_singletons(p = 200, G = 20, singletons = TRUE)
#' @export
generate_group_ind_singletons <- function(p, G, singletons = TRUE) {
  total_features <- p - 1

  if (singletons) {
    # Split first half into G groups, rest are singletons
    grouped_features <- floor(p / 2)
    singleton_features <- total_features - grouped_features

    split_sizes <- rep(floor(grouped_features / G), G)
    remainder <- grouped_features %% G
    if (remainder > 0) split_sizes[1:remainder] <- split_sizes[1:remainder] + 1

    group_ind <- c(rep(1:G, times = split_sizes),
                   (G + 1):(G + singleton_features))  # singleton groups get unique IDs
  } else {
    # divide all features into groups of equal size
    if (total_features %% G != 0) {
      warning("p - 1 is not divisible by G. Groups will be uneven.")
    }
    split_sizes <- rep(floor(total_features / G), G)
    remainder <- total_features %% G
    if (remainder > 0) split_sizes[1:remainder] <- split_sizes[1:remainder] + 1

    group_ind <- rep(1:G, times = split_sizes)
  }

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


