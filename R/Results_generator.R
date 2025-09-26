#' @title Local-level or cell-specific credible interval (CI)
#'
#' @param Beta is a dataframe of dimension (niter x N) with MCMC beta samples, for N spatial units
#' @param local.p.ths is the size of the CI for each location
#' @return a N x 1 vector of 0/1, depending on the CI for a location contains 0 or not.
#' @importFrom bayestestR ci
#'
#' @export
#'

local_credible <-function(Beta, local.p.ths = 0.9){
  CI_of_betas = t(apply(Beta, 2,
                        function(x) unlist(bayestestR::ci(x, ci = local.p.ths, method = "ETI"))[-1]))
  zero_finder = apply(CI_of_betas, 1,
                      function(x){ifelse(x[1] <= 0 & x[2] >= 0, 0, 1)})
  return(zero_finder)
}


#' @title MCMC convergence trace (vector input)
#' @description Plot a univariate MCMC trace and show the posterior mean,
#' Geweke two-sided p-value, and effective sample size (ESS) in the legend.
#'
#' @param x Numeric vector of MCMC samples (length = n_iter).
#' @param main Optional plot title.
#' @param legend_pos Legend position (e.g., "topright").
#' @param col_trace Line color for the trace.
#' @param col_mean Line color for the posterior mean.
#' @param lwd Line width for the trace.
#' @param cex_leg Legend text size.
#' @param ylim Optional y-limits; if NULL, set to +-1.2*max(|x|) with a fallback.
#'
#' @return Invisibly, a list with \code{mean}, \code{geweke_p}, \code{ess}, \code{n}, \code{ylim}.
#'
#'@examples
#' \dontrun{
#'  # simple vector
#'  # samples <- rnorm(2000, 0.5, 1)
#'  # MCMC_plot(samples, main = "Beta_4")
#'
#'  # with custom legend location and color
#'  # MCMC_plot(samples, main = "beta_4", legend_pos = "bottomright", col_trace = "blue")
#' }
#'
#' @importFrom graphics plot abline legend
#' @importFrom coda as.mcmc effectiveSize geweke.diag
#' @importFrom stats pnorm
#' @export
MCMC_plot <- function(
    x,
    main = NULL,
    legend_pos = "topright",
    col_trace = "green",
    col_mean  = "black",
    lwd = 1,
    cex_leg = 0.9,
    ylim = NULL
){
  # --- sanitize vector ---
  x_vec <- as.numeric(x)
  x_vec <- x_vec[is.finite(x_vec)]
  n <- length(x_vec)
  if (n == 0L) stop("x must contain at least one finite numeric value.")

  # --- diagnostics (convert vector -> coda::mcmc) ---
  x_mcmc <- coda::as.mcmc(matrix(x_vec, ncol = 1L))
  # Geweke z -> two-sided p
  gw <- try(coda::geweke.diag(x_mcmc), silent = TRUE)
  geweke_p <- if (inherits(gw, "try-error") || is.null(gw$z)) {
    NA_real_
  } else {
    2 * (1 - stats::pnorm(abs(as.numeric(gw$z))))
  }
  # ESS
  ess <- try(as.numeric(coda::effectiveSize(x_mcmc)), silent = TRUE)
  if (inherits(ess, "try-error") || length(ess) == 0L) ess <- NA_real_

  # --- y-limits ---
  if (is.null(ylim)) {
    m <- suppressWarnings(max(abs(x_vec), na.rm = TRUE))
    if (!is.finite(m) || m == 0) m <- 1
    ylim <- c(-1.2 * m, 1.2 * m)
  }

  # --- plot ---
  graphics::plot(x_vec, type = "l", col = col_trace, lwd = lwd,
       main = main, xlab = "Iteration", ylab = "Value",
       ylim = ylim)
  graphics::abline(h = mean(x_vec), col = col_mean, lwd = 2)

  # --- legend ---
  leg_txt <- c(
    sprintf("Geweke p = %s",
            ifelse(is.finite(geweke_p), format(round(geweke_p, 3), nsmall = 3), "NA")),
    sprintf("ESS = %s",
            ifelse(is.finite(ess), format(round(ess, 0), scientific = FALSE), "NA"))
  )
  graphics::legend(legend_pos, legend = leg_txt, bty = "n", cex = cex_leg)

  invisible(list(
    mean = mean(x_vec),
    geweke_p = geweke_p,
    ess = ess,
    n = n,
    ylim = ylim
  ))
}



#' @title Plot posterior credible intervals across outcomes
#'
#' @description
#' Generic plot of beta posterior ETI intervals across one or more outcomes
#' (e.g., cancers) for a shared predictor set. This function is mainly intended for the HS based methods (standard or goruped).
#'
#' @param outcomes Character vector of outcome keys (e.g., c("bladder","breast","prostate")).
#' @param outcome_labels Optional named character vector mapping key -> pretty name.
#'   If omitted, the keys themselves are used.
#' @param betas Optional named list of posterior draw matrices (one per outcome).
#'   Each matrix must be n_iter x p with column names = predictors (no intercept).
#'   If NULL, supply \code{load_beta()}.
#' @param load_beta Optional loader function \code{function(outcome_key, method) -> matrix or NULL}
#'   used when \code{betas} is NULL.
#' @param method Character label forwarded to \code{load_beta} (e.g., "HS" or "SS").
#' @param predictors Character vector of predictor names (defines order/inclusion; no intercept).
#' @param predictor_labels Optional named vector mapping predictor -> x-axis label.
#' @param group_map Optional named vector mapping predictor -> group name (for banding).
#' @param group_levels Optional character vector giving the order of groups under x-axis.
#' @param use_groups Logical; if FALSE, disables grouping entirely (ordering & brackets).
#' @param intercept_names Character vector of column names to drop if present.
#' @param cred_level Numeric ETI level in (0,1); default 0.95.
#' @param star "mean" or "median" for the point on each interval.
#' @param x_order "group" (grouped then predictor order) or "given" (exact \code{predictors} order).
#'   If \code{use_groups=FALSE}, this is coerced to "given".
#' @param colors Optional color spec. Either a single color (character scalar) or
#'   a named vector mapping pretty outcome names -> colors.
#' @param add_brackets Logical; add group brackets/labels beneath x-axis (only if grouping active).
#' @param legend_title Character; title for the outcome legend (default "Outcome").
#' @param show_legend NULL (auto: hide if one outcome), TRUE, or FALSE.
#' @param verbose Logical; print brief diagnostics.
#' @param base_size Base font size for theme (default 22).
#' @param bracket_text_size Font size for the group names (default 4.5, different scale then base_size)
#'
#' @return List with:
#' \item{plot}{ggplot object}
#' \item{data}{data.frame used to draw the plot}
#' \item{loaded}{diagnostics for loaded matrices}
#' \item{skipped}{diagnostics for skipped outcomes}
#'
#' @examples
#' \dontrun{
#' # Example with a single outcome (legend auto-hidden)
#' beta_mat <- matrix(rnorm(4000), ncol=4, dimnames=list(NULL, c("a","b","c","d")))
#' res1 <- plot_cri_across_outcomes(
#'   outcomes = "breast",
#'   outcome_labels = c(breast="Breast (female)"),
#'   betas = list(breast = beta_mat),
#'   predictors = c("a","b","c","d"),
#'   predictor_labels = c(a="A", b="B", c="C", d="D"),
#'   use_groups = FALSE,
#'   legend_title = "Cancer"
#' )
#' print(res1$plot)
#'
#' # Multiple outcomes with grouping and custom legend title
#' betas <- list(
#'   bladder  = beta_mat,
#'   breast   = beta_mat,
#'   prostate = beta_mat
#' )
#' res2 <- plot_cri_across_outcomes(
#'   outcomes = c("bladder","breast","prostate"),
#'   outcome_labels = c(bladder="Bladder", breast="Breast (female)", prostate="Prostate (male)"),
#'   betas = betas,
#'   predictors = c("a","b","c","d"),
#'   predictor_labels = c(a="A", b="B", c="C", d="D"),
#'   group_map = c(a="Group 1", b="Group 1", c="Group 2", d="Group 2"),
#'   group_levels = c("Group 1","Group 2"),
#'   use_groups = TRUE,
#'   colors = c("Bladder"="#fec601","Breast (female)"="#FF3399","Prostate (male)"="#69B3E7"),
#'   legend_title = "Cancer Type"
#' )
#' print(res2$plot)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_linerange geom_point geom_hline scale_color_manual scale_fill_manual
#' @importFrom ggplot2 labs theme_minimal theme element_text element_blank position_dodge
#' @importFrom ggplot2 element_rect annotate coord_cartesian expand_limits margin
#' @importFrom dplyr bind_rows distinct arrange mutate group_by summarise
#' @importFrom bayestestR ci
#' @importFrom stats setNames median
#' @importFrom rlang .data
#' @export
plot_cri_across_outcomes <- function(
    outcomes,
    outcome_labels = NULL,
    betas = NULL,
    load_beta = NULL,
    method = "HS",
    predictors,
    predictor_labels = NULL,
    group_map = NULL,
    group_levels = NULL,
    use_groups = TRUE,
    intercept_names = c("(Intercept)","Intercept","beta0","beta_0",""),
    cred_level = 0.95,
    star = c("mean","median"),
    x_order = c("group","given"),
    colors = NULL,
    add_brackets = TRUE,
    legend_title = "Outcome",
    show_legend = NULL,
    verbose = TRUE,
    base_size = 22,
    bracket_text_size =4.5
){
  if (!requireNamespace("bayestestR", quietly = TRUE))
    stop("Package 'bayestestR' is required. Please install it.")

  star    <- match.arg(star)
  x_order <- match.arg(x_order)
  use_groups <- isTRUE(use_groups)

  # Pretty outcome names
  pretty <- if (is.null(outcome_labels)) setNames(outcomes, outcomes) else outcome_labels
  if (!all(outcomes %in% names(pretty)))
    stop("All outcomes must have labels in 'outcome_labels' (or omit to use keys).")

  # Input sources
  if (is.null(betas)) {
    if (!is.function(load_beta))
      stop("Provide either 'betas' (named list) or a 'load_beta(outcome, method)' function.")
  } else {
    if (!all(outcomes %in% names(betas)))
      stop("Names(betas) must include all 'outcomes'.")
  }

  # ---- Helpers ----
  normalize_beta <- function(M, predictors, intercept_names) {
    M <- as.matrix(M)

    # If no names, infer from predictors (+ optional intercept)
    if (is.null(colnames(M))) {
      if (ncol(M) == length(predictors)) {
        colnames(M) <- predictors
      } else if (ncol(M) == length(predictors) + 1) {
        colnames(M) <- c("(Intercept)", predictors)
      } else {
        stop(sprintf("Cannot infer beta column names (ncol=%d; |predictors|=%d).",
                     ncol(M), length(predictors)))
      }
    }

    # orientation check still works
    rn <- rownames(M); cn <- colnames(M)
    rn_match <- if (!is.null(rn)) sum(rn %in% predictors) else 0
    cn_match <- if (!is.null(cn)) sum(cn %in% predictors) else 0
    if (rn_match > cn_match) {
      M <- t(M)
    }

    # drop intercept-like
    keep <- !(colnames(M) %in% intercept_names)
    M <- M[, keep, drop = FALSE]

    # align to predictors, fill gaps
    common  <- intersect(predictors, colnames(M))
    missing <- setdiff(predictors, colnames(M))
    M2 <- cbind(
      M[, common, drop = FALSE],
      if (length(missing)) matrix(NA_real_, nrow(M), length(missing),
                                  dimnames = list(NULL, missing)) else NULL
    )
    M2[, predictors, drop = FALSE]
  }


  add_group_brackets_scaled_ <- function(
    p, data_df,
    offset_bracket = 0.06,  # fraction of span below min
    offset_text    = 0.10,
    tick_frac      = 0.02,
    line_col = "grey35", line_size = 0.6, text_size = bracket_text_size
  ){
    yr <- range(c(data_df$q_low, data_df$q_high, data_df$mean, data_df$median), na.rm = TRUE)
    span <- diff(yr); if (!is.finite(span) || span <= 0) span <- max(1, abs(yr[1L]))
    y_b <- yr[1L] - offset_bracket * span
    y_t <- yr[1L] - offset_text    * span
    tick_h <- tick_frac * span

    df <- dplyr::distinct(data_df, .data$predictor, .data$group_lbl) |>
      dplyr::arrange(.data$predictor) |>
      dplyr::mutate(x = as.numeric(.data$predictor))
    ranges <- df |>
      dplyr::filter(!is.na(.data$group_lbl)) |>
      dplyr::group_by(.data$group_lbl) |>
      dplyr::summarise(xmin = min(.data$x), xmax = max(.data$x), .groups = "drop")

    p +
      ggplot2::geom_segment(
        data = ranges,
        ggplot2::aes(x = .data$xmin - 0.45, xend = .data$xmax + 0.45, y = y_b, yend = y_b),
        inherit.aes = FALSE, linewidth = line_size, color = line_col
      ) +
      ggplot2::geom_segment(
        data = ranges,
        ggplot2::aes(x = .data$xmin - 0.45, xend = .data$xmin - 0.45, y = y_b, yend = y_b + tick_h),
        inherit.aes = FALSE, linewidth = line_size, color = line_col
      ) +
      ggplot2::geom_segment(
        data = ranges,
        ggplot2::aes(x = .data$xmax + 0.45, xend = .data$xmax + 0.45, y = y_b, yend = y_b + tick_h),
        inherit.aes = FALSE, linewidth = line_size, color = line_col
      ) +
      ggplot2::annotate(
        "text",
        x = (ranges$xmin + ranges$xmax)/2, y = y_t,
        label = ranges$group_lbl, size = text_size, color = "grey20"
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::expand_limits(y = y_t) +
      ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 20))
  }

  # ---- Grouping gate ----
  if (!use_groups || is.null(group_map) || is.null(group_levels)) {
    group_map <- NULL
    group_levels <- NULL
    if (x_order == "group") {
      if (verbose) message("Grouping disabled: coercing x_order='group' to 'given'.")
      x_order <- "given"
    }
  }

  # X order
  if (x_order == "group") {
    ord_df <- data.frame(
      predictor = predictors,
      group_lbl = unname(group_map[predictors]),
      stringsAsFactors = FALSE
    )
    ord_df$group_lbl <- factor(ord_df$group_lbl, levels = group_levels)
    ord_df <- ord_df[order(ord_df$group_lbl, match(ord_df$predictor, predictors)), ]
    pred_order <- ord_df$predictor
  } else {
    pred_order <- predictors
  }

  # ---- Collect across outcomes ----
  rows <- list(); loaded <- list(); skipped <- list()
  for (k in outcomes) {
    M <- if (!is.null(betas)) betas[[k]] else load_beta(k, method)
    if (is.null(M)) { skipped[[length(skipped)+1]] <- data.frame(outcome=k, reason="not found"); next }
    if (nrow(M) < 20L) { skipped[[length(skipped)+1]] <- data.frame(outcome=k, reason="<20 iterations"); next }

    M <- try(normalize_beta(M, predictors, intercept_names), silent = TRUE)
    if (inherits(M, "try-error")) {
      skipped[[length(skipped)+1]] <- data.frame(outcome=k, reason="column name mismatch"); next
    }

    # Named ETIs per predictor
    cn <- colnames(M)
    get_ci_bounds <- function(v) {
      ci_df <- bayestestR::ci(v, ci = cred_level, method = "ETI")
      c(low  = min(ci_df$CI_low, ci_df$CI_high, na.rm = TRUE),
        high = max(ci_df$CI_low, ci_df$CI_high, na.rm = TRUE))
    }
    ci_mat <- sapply(cn, function(nm) get_ci_bounds(M[, nm]))  # 2 x p, colnames = cn
    q_low  <- stats::setNames(ci_mat["low",  ], cn)
    q_high <- stats::setNames(ci_mat["high", ], cn)

    lab <- pretty[[k]]
    rows[[length(rows)+1]] <- data.frame(
      Outcome   = lab,
      predictor = factor(pred_order, levels = pred_order),
      q_low     = as.numeric(q_low [pred_order]),
      q_high    = as.numeric(q_high[pred_order]),
      median    = apply(M[, pred_order, drop = FALSE], 2, median,   na.rm = TRUE),
      mean      = colMeans(M[, pred_order, drop = FALSE],           na.rm = TRUE),
      group_lbl = if (!is.null(group_map))
        factor(unname(group_map[pred_order]), levels = group_levels)
      else factor(NA_character_, levels = NULL),
      stringsAsFactors = FALSE
    )

    loaded[[length(loaded)+1]] <- data.frame(
      outcome_key = k, pretty = lab, n_iter = nrow(M), n_pred = ncol(M),
      stringsAsFactors = FALSE
    )
  }

  if (!length(rows)) stop("No outcomes loaded. Check inputs.")
  df_all <- dplyr::bind_rows(rows)
  diag_loaded  <- if (length(loaded))  dplyr::bind_rows(loaded)  else data.frame()
  diag_skipped <- if (length(skipped)) dplyr::bind_rows(skipped) else data.frame()

  if (verbose && nrow(diag_loaded))  print(diag_loaded)
  if (verbose && nrow(diag_skipped)) message("Skipped:\n", paste(utils::capture.output(print(diag_skipped)), collapse = "\n")
)

  # star statistic
  df_all$star_y <- if (star == "mean") df_all$mean else df_all$median

  # legend visibility (auto-hide if 1 outcome unless overridden)
  present_outcomes <- as.character(unique(df_all$Outcome))
  legend_on <- if (is.null(show_legend)) length(present_outcomes) > 1 else isTRUE(show_legend)

  # x-axis labels
  xlabs <- if (!is.null(predictor_labels)) {
    # fill any missing labels with the predictor name itself
    labvec <- predictor_labels[as.character(levels(df_all$predictor))]
    labvec[is.na(labvec)] <- as.character(levels(df_all$predictor))[is.na(labvec)]
    labvec
  } else {
    levels(df_all$predictor)
  }

  # ---- Plot ----
  pd <- ggplot2::position_dodge(width = 0.78)
  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = .data$predictor, color = .data$Outcome)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = .data$q_low, ymax = .data$q_high), position = pd, linewidth = 1.1) +
    ggplot2::geom_point(ggplot2::aes(y = .data$star_y, fill = .data$Outcome),
                        position = pd, shape = 8, size = 2.8, stroke = 0.7,
                        show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::labs(
      x = NULL,
      y = sprintf("Posterior Beta (%d%% ETI; star = %s)", round(cred_level*100), tolower(star)),
      color = legend_title
    ) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = if (legend_on) "bottom" else "none",
      legend.box = "horizontal",
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
      axis.ticks.x = ggplot2::element_line(color = "grey50"),
      axis.ticks.length.x = grid::unit(3, "pt"),
      panel.background   = ggplot2::element_rect(fill = "white", color = NA),
      plot.background    = ggplot2::element_rect(fill = "white", color = NA),
      legend.background  = ggplot2::element_rect(fill = "white", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::scale_x_discrete(labels = xlabs)

  # Colors (single or named)
  if (!is.null(colors)) {
    if (is.null(names(colors)) && length(colors) == 1L) {
      colors <- stats::setNames(colors, present_outcomes)
    } else {
      colors <- colors[names(colors) %in% present_outcomes]
    }
    if (length(colors)) {
      p <- p +
        ggplot2::scale_color_manual(values = colors, limits = names(colors), drop = FALSE, name = legend_title) +
        ggplot2::scale_fill_manual(values = colors, limits = names(colors), drop = FALSE, name = legend_title)
    }
  }

  # Group brackets (scaled to data range so they don't get clipped)
  if (add_brackets && use_groups && !is.null(group_map)) {
    p <- add_group_brackets_scaled_(p, df_all,text_size = bracket_text_size)
  }

  list(plot = p, data = df_all, loaded = diag_loaded, skipped = diag_skipped)
}



#' @title Plot PIPs (pDelta) from MCMC outputs
#' @description
#' Plot Posterior Inclusion Probabilities (PIPs) per predictor, across one or more
#' outcomes, computed as column means of MCMC indicator matrices. Supply either:
#' (i) a single matrix `pdelta_mat` (iterations × predictors) with an `outcome_name`,
#' (i) a single matrix `pdelta_mat` (iterations x predictors) with an `outcome_name`,
#' or (ii) a named list `pdelta_list` of such matrices (one per outcome).
#' The argument `predictors` fixes x-axis order and inclusion.
#'
#' @param predictors Character vector of predictor keys (defines x-axis order/inclusion).
#' @param pdelta_mat Optional numeric matrix (iterations × predictors) for a single outcome.
#' @param pdelta_mat Optional numeric matrix (iterations x predictors) for a single outcome.
#' @param outcome_name Optional single character label for `pdelta_mat` (required if `pdelta_mat` is used).
#' @param pdelta_list Optional **named list** of numeric matrices (iterations × predictors), one per outcome.
#' @param pdelta_list Optional **named list** of numeric matrices (iterations x predictors), one per outcome.
#'   List names are used as outcome labels. If provided, `pdelta_mat`/`outcome_name` are ignored.
#' @param outcomes Optional character vector to subset outcomes (must match names in `pdelta_list`
#'   or `outcome_name` for the single-matrix case).
#' @param outcome_labels Optional named vector mapping outcome key -> pretty label used in the legend.
#'   Defaults to identity mapping for present outcomes.
#' @param predictor_labels Optional named vector mapping predictor key -> pretty x-axis label.
#'   Defaults to the keys.
#' @param group_map Optional named vector mapping predictor key -> group label (for grouping bands).
#' @param group_levels Optional character vector giving the order of the groups (required if grouping used).
#' @param use_groups Logical; if FALSE, turns grouping off entirely (ordering + bands). Default TRUE.
#' @param x_order "group" (grouped then predictor order) or "given" (exact `predictors` order).
#'   If `use_groups=FALSE`, coerced to "given".
#' @param thresholds Numeric vector of horizontal reference lines (e.g., `c(0.5, 0.9)`).
#' @param threshold_colors Optional same-length vector of colors for `thresholds`.
#' @param colors Optional color spec for outcomes: either a single color, or a named vector
#'   mapping pretty outcome names -> colors.
#' @param legend_title Character; title for the outcome legend (default "Outcome").
#' @param show_legend NULL (auto: hide if only one outcome), TRUE, or FALSE.
#' @param bar_width Width of bars (default 0.75).
#' @param dodge_width Dodge width for grouped bars (default 0.8).
#' @param base_size Base font size for theme (default 22).
#' @param verbose Logical; print brief diagnostics (default TRUE).
#' @param bracket_text_size Font size for the group names (default 4.5, different scale then base_size)
#'
#' @return A list with \item{plot}{ggplot object} and \item{data}{data.frame used for plotting}.
#'
#' @examples
#' \dontrun{
#' # Single outcome
#' set.seed(1)
#' predictors <- paste0("x", 1:6)
#' pdelta_mat <- matrix(rbinom(2000, 1, prob = runif(6, 0.1, 0.9)), ncol = 6)
#' res1 <- plot_pdelta_bars_mcmc(
#'   predictors   = predictors,
#'   pdelta_mat   = pdelta_mat,
#'   outcome_name = "Breast (female)",
#'   thresholds   = c(0.5, 0.9)
#' )
#' print(res1$plot)
#'
#' # Multiple outcomes
#' pdelta_list <- list(
#'   "Breast (female)" = matrix(rbinom(3000, 1, prob = runif(6, 0.2, 0.8)), ncol = 6),
#'   "Prostate (male)" = matrix(rbinom(3000, 1, prob = runif(6, 0.1, 0.6)), ncol = 6)
#' )
#' res2 <- plot_pdelta_bars_mcmc(
#'   predictors       = predictors,
#'   pdelta_list      = pdelta_list,
#'   predictor_labels = setNames(toupper(predictors), predictors),
#'   group_map        = c(x1="Demog", x2="Demog", x3="Behavior", x4="Behavior", x5="Access", x6="Access"),
#'   group_map = c(x1="Demog", x2="Demog", x3="Behavior",
#'   x4="Behavior", x5="Access", x6="Access"),
#'   group_levels     = c("Demog","Behavior","Access"),
#'   thresholds       = c(0.5, 0.9),
#'   legend_title     = "Cancer"
#' )
#' print(res2$plot)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_col geom_hline scale_fill_manual
#' @importFrom ggplot2 scale_x_discrete labs theme_minimal theme element_text element_blank
#' @importFrom ggplot2 position_dodge coord_cartesian annotate expand_limits element_rect margin
#' @importFrom dplyr distinct arrange mutate group_by summarise
#' @importFrom stats setNames
#' @importFrom rlang .data
#' @export
plot_pdelta_bars_mcmc <- function(
    predictors,
    pdelta_mat      = NULL,
    outcome_name    = NULL,
    pdelta_list     = NULL,
    outcomes        = NULL,
    outcome_labels  = NULL,
    predictor_labels= NULL,
    group_map       = NULL,
    group_levels    = NULL,
    use_groups      = TRUE,
    x_order         = c("group","given"),
    thresholds      = 0.75,
    threshold_colors= NULL,
    colors          = NULL,
    legend_title    = "Outcome",
    show_legend     = NULL,
    bar_width       = 0.75,
    dodge_width     = 0.8,
    base_size       = 22,
    verbose         = TRUE,
    bracket_text_size = 4.5
){
  x_order   <- match.arg(x_order)
  use_groups <- isTRUE(use_groups)

  if (missing(predictors) || !length(predictors)) stop("`predictors` must be provided and non-empty.")

  # ---- Build tidy data from matrices (colMeans of pDelta)
  if (!is.null(pdelta_list)) {
    if (is.null(names(pdelta_list)) || any(names(pdelta_list) == "")) {
      stop("`pdelta_list` must be a *named* list; names are used as outcome labels.")
    }
    # Subset outcomes if requested
    list_use <- if (is.null(outcomes)) pdelta_list else {
      miss <- setdiff(outcomes, names(pdelta_list))
      if (length(miss)) stop("Requested outcomes not in `pdelta_list`: ", paste(miss, collapse = ", "))
      pdelta_list[outcomes]
    }
    built <- lapply(names(list_use), function(out) {
      mat <- list_use[[out]]
      if (!is.matrix(mat)) stop("Each element of `pdelta_list` must be a matrix (iterations × predictors).")
      if (!is.matrix(mat)) stop("Each element of `pdelta_list` must be a matrix (iterations x predictors).")
      if (ncol(mat) != length(predictors)) {
        stop(sprintf("Outcome '%s': ncol(matrix) = %d != length(predictors) = %d",
                     out, ncol(mat), length(predictors)))
      }
      data.frame(
        Outcome   = out,
        predictor = predictors,
        pDelta    = colMeans(mat, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
    d <- do.call(rbind, built)

  } else if (!is.null(pdelta_mat)) {
    if (is.null(outcome_name) || length(outcome_name) != 1L) {
      stop("Provide a single `outcome_name` when using `pdelta_mat`.")
    }
    if (!is.matrix(pdelta_mat)) stop("`pdelta_mat` must be a matrix (iterations × predictors).")
    if (!is.matrix(pdelta_mat)) stop("`pdelta_mat` must be a matrix (iterations x predictors).")
    if (ncol(pdelta_mat) != length(predictors)) {
      stop(sprintf("ncol(pdelta_mat) = %d != length(predictors) = %d",
                   ncol(pdelta_mat), length(predictors)))
    }
    if (!is.null(outcomes) && !(outcome_name %in% outcomes)) {
      stop("`outcome_name` not included in requested `outcomes`.")
    }
    d <- data.frame(
      Outcome   = outcome_name,
      predictor = predictors,
      pDelta    = colMeans(pdelta_mat, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

  } else {
    stop("Provide either `pdelta_mat` (with `outcome_name`) or `pdelta_list`.")
  }

  # ---- Pretty outcome names
  present_keys <- unique(d$Outcome)
  pretty <- if (is.null(outcome_labels)) stats::setNames(present_keys, present_keys) else outcome_labels
  if (!all(present_keys %in% names(pretty))) {
    missing <- setdiff(present_keys, names(pretty))
    stop("Missing labels for outcomes: ", paste(missing, collapse = ", "))
  }
  d$Outcome <- pretty[match(d$Outcome, names(pretty))]

  # NA/Inf pDelta -> 0 for plotting
  d$pDelta[!is.finite(d$pDelta)] <- NA_real_
  d$pDelta <- ifelse(is.na(d$pDelta), 0, d$pDelta)

  # ---- Grouping gate & x order
  if (!use_groups || is.null(group_map) || is.null(group_levels)) {
    group_map <- NULL
    group_levels <- NULL
    if (x_order == "group") {
      if (verbose) message("Grouping disabled: coercing x_order='group' to 'given'.")
      x_order <- "given"
    }
  }

  # predictor order
  if (x_order == "group") {
    ord_df <- data.frame(
      predictor = predictors,
      group_lbl = unname(group_map[predictors]),
      stringsAsFactors = FALSE
    )
    ord_df$group_lbl <- factor(ord_df$group_lbl, levels = group_levels)
    ord_df <- ord_df[order(ord_df$group_lbl, match(ord_df$predictor, predictors)), ]
    pred_order <- ord_df$predictor
  } else {
    pred_order <- predictors
  }

  # complete grid (ensures 0s for missing combos)
  all_df <- expand.grid(
    Outcome   = unique(d$Outcome),
    predictor = predictors,
    stringsAsFactors = FALSE
  )
  d <- merge(all_df, d, by = c("Outcome","predictor"), all.x = TRUE)
  d$pDelta <- ifelse(is.na(d$pDelta), 0, d$pDelta)

  d$predictor <- factor(d$predictor, levels = pred_order)
  if (!is.null(group_map)) {
    d$group_lbl <- factor(unname(group_map[as.character(d$predictor)]), levels = group_levels)
  } else {
    d$group_lbl <- factor(NA_character_)
  }

  # x-axis labels
  xlabs <- if (!is.null(predictor_labels)) {
    lab <- predictor_labels[as.character(levels(d$predictor))]
    lab[is.na(lab)] <- as.character(levels(d$predictor))[is.na(lab)]
    lab
  } else {
    levels(d$predictor)
  }

  # legend visibility
  present_outcomes <- as.character(unique(d$Outcome))
  legend_on <- if (is.null(show_legend)) length(present_outcomes) > 1 else isTRUE(show_legend)

  # thresholds styling
  if (is.null(threshold_colors)) {
    threshold_colors <- rep("firebrick", length(thresholds))
  } else if (length(threshold_colors) == 1L && length(thresholds) > 1L) {
    threshold_colors <- rep(threshold_colors, length(thresholds))
  }

  # ----- plot
  pd <- ggplot2::position_dodge(width = dodge_width)
  sep_x <- seq(0.5, length(levels(d$predictor)) - 0.5, by = 1)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$predictor, y = .data$pDelta, fill = .data$Outcome)) +
    ggplot2::geom_vline(xintercept = sep_x, color = "grey92", linewidth = 0.4) +
    ggplot2::geom_col(position = pd, width = bar_width, color = "white", linewidth = 0.2) +
    {
      layer_list <- list()
      if (length(thresholds)) {
        for (i in seq_along(thresholds)) {
          layer_list[[i]] <- ggplot2::geom_hline(yintercept = thresholds[i],
                                                 linetype = "dashed",
                                                 color = threshold_colors[i],
                                                 linewidth = 0.6)
        }
      }
      layer_list
    } +
    ggplot2::labs(x = NULL, y = "PIP", fill = legend_title) +
    ggplot2::scale_x_discrete(labels = xlabs) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = if (legend_on) "bottom" else "none",
      legend.box = "horizontal",
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
      axis.ticks.x = ggplot2::element_line(color = "grey50"),
      axis.ticks.length.x = grid::unit(3, "pt"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.background   = ggplot2::element_rect(fill = "white", color = NA),
      plot.background    = ggplot2::element_rect(fill = "white", color = NA),
      legend.background  = ggplot2::element_rect(fill = "white", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.12, max(1, max(d$pDelta, na.rm = TRUE))), clip = "off")

  # colors for outcomes
  if (!is.null(colors)) {
    if (is.null(names(colors)) && length(colors) == 1L) {
      colors <- stats::setNames(colors, present_outcomes)
    } else {
      colors <- colors[names(colors) %in% present_outcomes]
    }
    if (length(colors)) {
      p <- p + ggplot2::scale_fill_manual(values = colors, limits = names(colors), drop = FALSE, name = legend_title)
    }
  }

  # ---- optional group brackets
  if (use_groups && !is.null(group_map)) {
    p <- add_group_brackets_scaled_(p, d, text_size = bracket_text_size)
  }

  list(plot = p, data = d)
}

# internal helper: add group brackets using scaled y placement
add_group_brackets_scaled_ <- function(
    p, data_df,
    y_frac = 0.10, tick_frac = 0.04, text_frac = 0.15,
    line_col = "grey35", line_size = 0.6, text_size = bracket_text_size
){
  ymax <- max(data_df$pDelta, na.rm = TRUE); if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  y_min  <- -y_frac  * ymax
  tick_h <-  tick_frac * ymax
  y_txt  <- -text_frac * ymax

  df <- dplyr::distinct(data_df, .data$predictor, .data$group_lbl) |>
    dplyr::arrange(.data$predictor) |>
    dplyr::mutate(x = as.numeric(.data$predictor))

  ranges <- df |>
    dplyr::filter(!is.na(.data$group_lbl)) |>
    dplyr::group_by(.data$group_lbl) |>
    dplyr::summarise(xmin = min(.data$x), xmax = max(.data$x), .groups = "drop")

  p +
    ggplot2::geom_segment(data = ranges,
                          ggplot2::aes(x = .data$xmin - 0.45, xend = .data$xmax + 0.45,
                                       y = y_min, yend = y_min),
                          inherit.aes = FALSE, linewidth = line_size, color = line_col
    ) +
    ggplot2::geom_segment(
      data = ranges,
      ggplot2::aes(x = .data$xmin - 0.45, xend = .data$xmin - 0.45, y = y_min, yend = y_min + tick_h),
      inherit.aes = FALSE, linewidth = line_size, color = line_col
    ) +
    ggplot2::geom_segment(
      data = ranges,
      ggplot2::aes(x = .data$xmax + 0.45, xend = .data$xmax + 0.45, y = y_min, yend = y_min + tick_h),
      inherit.aes = FALSE, linewidth = line_size, color = line_col
    ) +
    ggplot2::annotate(
      "text",
      x = (ranges$xmin + ranges$xmax)/2, y = y_txt,
      label = ranges$group_lbl, size = text_size, color = "grey20"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::expand_limits(y = y_txt) +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 20))
}
