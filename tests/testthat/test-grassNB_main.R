# tests/testthat/test-MUSC_main.R

skip_if_not_installed("spdep")

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

# Tiny nb objects for dimension/class checks
nb        <- spdep::cell2nb(3, 3, type = "queen")
n         <- length(nb)
X         <- cbind(matrix(1, n, 2))   # 2 cols of 1s (shape only)
y         <- rep(1L, n)
grp       <- c(1, 2)
niter     <- 10                        # keep tests fast

# ------------------------------------------------------------------
# Dimension checks
# ------------------------------------------------------------------
test_that("dimension checks fail with clear messages", {
  set.seed(1)

  expect_error(
    silently(MUSC(y = y[-1], X = X, NeighborhoodList = nb, group_ind = grp, niter = niter)),
    "Length of y and number of rows of X do not match",
    fixed = TRUE
  )

  nb_mismatch <- spdep::cell2nb(4, 4, type = "queen") # still class 'nb'
  expect_error(
    silently(MUSC(y = y, X = X, NeighborhoodList = nb_mismatch, group_ind = grp, niter = niter)),
    "Length of y and length of NeighborhoodList do not match",
    fixed = TRUE
  )

  expect_error(
    silently(MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = c(0, 3, 1), niter = niter)),
    "Length of group_ind and number of columns of X do not match",
    fixed = TRUE
  )
})

# ------------------------------------------------------------------
# NeighborhoodList class
# ------------------------------------------------------------------
test_that("NeighborhoodList must be class 'nb'", {
  bad_nb <- replicate(n, integer(0), simplify = FALSE) # not 'nb'
  expect_error(
    silently(MUSC(y = y, X = X, NeighborhoodList = bad_nb, group_ind = grp, niter = niter)),
    "'NeighborhoodList' must be an object of class 'nb' (e.g., from spdep::poly2nb).",
    fixed = TRUE
  )
})

# ------------------------------------------------------------------
# Rate argument pairing & validation
# ------------------------------------------------------------------
test_that("rate arguments must travel together and be valid", {
  set.seed(2)
  pop <- sample(10000:50000, n, TRUE)

  expect_error(
    silently(MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
                  pop_col = pop, scale = NULL)),
    "Both 'pop_col' and 'scale' must be provided together"
  )

  expect_error(
    silently(MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
                  pop_col = NULL, scale = 100000)),
    "Both 'pop_col' and 'scale' must be provided together"
  )

  expect_error(
    silently(MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
                  pop_col = pop, scale = 1.5)),
    "'scale' must be a single positive integer (e.g., 100000).",
    fixed = TRUE
  )

  expect_error(
    silently(MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
                  pop_col = pop[-1], scale = 100000)),
    "Length of 'pop_col' must equal the number of rows of X (and length of y).",
    fixed = TRUE
  )

  bad_pop <- replace(pop, 1, 0)
  expect_error(
    silently(MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
                  pop_col = bad_pop, scale = 100000)),
    "positive and finite"
  )

  expect_silent(
    silently(MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
                  pop_col = pop, scale = 100000, verbose = FALSE))
  )
})

# ------------------------------------------------------------------
# Count-data sanity warnings (mock children so we don't descend)
# ------------------------------------------------------------------
test_that("count-data sanity warnings fire only when verbose", {
  y_nonint <- rep(1.25, n)
  y_unique <- as.numeric(seq_len(n))

  testthat::with_mocked_bindings(
    VS_Group           = function(...) list(mock = "VS_Group"),
    VS_Standard        = function(...) list(mock = "VS_Standard"),
    VS_Group_offset    = function(...) list(mock = "VS_Group_offset"),
    VS_Standard_offset = function(...) list(mock = "VS_Standard_offset"),
    {
      expect_warning(
        silently(MUSC(y = y_nonint, X = X, NeighborhoodList = nb, group_ind = grp,
                      niter = niter, verbose = TRUE)),
        "Outcome 'y' should be nonnegative integers for counts; found non-integer or negative values.",
        fixed = TRUE
      )

      expect_warning(
        silently(MUSC(y = y_unique, X = X, NeighborhoodList = nb, group_ind = grp,
                      niter = niter, verbose = TRUE)),
        "The data is possibly continuous; please double-check that the outcome is a discrete count.",
        fixed = TRUE
      )

      expect_silent(
        silently(MUSC(y = y_nonint, X = X, NeighborhoodList = nb, group_ind = grp,
                      niter = niter, verbose = FALSE))
      )
      expect_silent(
        silently(MUSC(y = y_unique, X = X, NeighborhoodList = nb, group_ind = grp,
                      niter = niter, verbose = FALSE))
      )
    }
  )
})


# skip_if_not_installed("spdep")
#
# nb  <- spdep::cell2nb(3, 3, type = "queen")
# n   <- length(nb)
# X   <- cbind(matrix(1, n, 2))   # 2 cols of 1s (ok for shape checks)
# y   <- rep(1L, n)
# grp <- c(1, 2)
# niter <- 10000
#
# test_that("dimension checks fail with clear messages", {
#   set.seed(1)
#
#   expect_error(
#     MUSC(y = y[-1], X = X, NeighborhoodList = nb, group_ind = grp, niter = niter),
#     "Length of y and number of rows of X do not match", fixed = TRUE
#   )
#
#   nb_mismatch <- spdep::cell2nb(4, 4, type = "queen") # still class 'nb'
#   expect_error(
#     MUSC(y = y, X = X, NeighborhoodList = nb_mismatch, group_ind = grp, niter = niter),
#     "Length of y and length of NeighborhoodList do not match", fixed = TRUE
#   )
#
#   expect_error(
#     MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = c(0, 3, 1), niter = niter),
#     "Length of group_ind and number of columns of X do not match", fixed = TRUE
#   )
# })
#
# test_that("NeighborhoodList must be class 'nb'", {
#   bad_nb <- replicate(n, integer(0), simplify = FALSE) # not 'nb'
#   expect_error(
#     invisible(capture.output(
#       MUSC(y = y, X = X, NeighborhoodList = bad_nb, group_ind = grp, niter = niter)
#     )),
#     "'NeighborhoodList' must be an object of class 'nb' (e.g., from spdep::poly2nb).",
#     fixed = TRUE
#   )
#
#
# })
#
# test_that("rate arguments must travel together and be valid", {
#   set.seed(2)
#   pop <- sample(10000:50000, n, TRUE)
#
#   expect_error(
#     invisible(capture.output(
#       MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
#            pop_col = pop, scale = NULL)
#     )),
#     "Both 'pop_col' and 'scale' must be provided together",                # substring
#   )
#   expect_error(
#     MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
#          pop_col = NULL, scale = 100000),
#     "Both 'pop_col' and 'scale' must be provided together",
#   )
#   expect_error(
#     MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
#          pop_col = pop, scale = 1.5),
#     "'scale' must be a single positive integer (e.g., 100000).", fixed = TRUE
#   )
#
#   expect_error(
#     MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
#          pop_col = pop[-1], scale = 100000),
#     "Length of 'pop_col' must equal the number of rows of X (and length of y).",
#     fixed = TRUE
#   )
#
#   bad_pop <- replace(pop, 1, 0)
#   expect_error(
#     MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
#          pop_col = bad_pop, scale = 100000),
#     "positive and finite"   # substring is fine
#   )
#
#   expect_silent(
#     MUSC(y = y, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
#          pop_col = pop, scale = 100000, verbose = FALSE)
#   )
# })
#
# test_that("count-data sanity warnings fire only when verbose", {
#   y_nonint <- rep(1.25, n)
#   expect_error(
#     invisible(capture.output(
#       MUSC(y = y_nonint, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
#            verbose = TRUE)
#     )),
#     "Outcome 'y' should be nonnegative integers for counts; found non-integer or negative values.",
#     fixed = TRUE
#   )
#
#   y_unique <- as.numeric(seq_len(n))
#   expect_warning(
#     MUSC(y = y_unique, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
#          verbose = TRUE),
#     "The data is possibly continuous; please double-check that the outcome is a discrete count.",
#     fixed = TRUE
#   )
#
#   expect_silent(
#     MUSC(y = y_nonint, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
#          verbose = FALSE)
#   )
#   expect_silent(
#     MUSC(y = y_unique, X = X, NeighborhoodList = nb, group_ind = grp, niter = niter,
#          verbose = FALSE)
#   )
# })
#
#
