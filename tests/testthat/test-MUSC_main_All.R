# --- Routing & arg forwarding for all 4 branches -------------------
test_that("routes to correct child and forwards args", {
  nb  <- spdep::cell2nb(2, 2, type = "queen")
  X   <- matrix(rnorm(length(nb) * 2), ncol = 2)
  y   <- seq_along(nb)
  grp <- c(1, 2)
  pop <- rep(10000, length(nb))

  called <- new.env(parent = emptyenv())
  reset_called <- function() {
    rm(list = ls(envir = called), envir = called)
    called$which <- NULL
  }

  make_mock <- function(which_name) {
    function(K, y, NeighborhoodList, which.prior, niter, verbose, pop_col = NULL, scale = NULL, ...) {
      called$which         <- which_name
      called$K             <- K
      called$y             <- y
      called$NeighborhoodList <- NeighborhoodList
      called$which.prior   <- which.prior
      called$niter         <- niter
      called$verbose       <- verbose
      called$pop_col       <- pop_col
      called$scale         <- scale
      if (!missing(...)) called$dots <- list(...)
      list(ok = TRUE)
    }
  }

  testthat::with_mocked_bindings(
    VS_Standard        = make_mock("VS_Standard"),
    VS_Standard_offset = make_mock("VS_Standard_offset"),
    VS_Group           = function(K, y, group_ind, NeighborhoodList, which.prior, niter, verbose, pop_col = NULL, scale = NULL, ...) {
      # identical to make_mock but keeps group_ind
      called$which         <- "VS_Group"
      called$K             <- K
      called$y             <- y
      called$group_ind     <- group_ind
      called$NeighborhoodList <- NeighborhoodList
      called$which.prior   <- which.prior
      called$niter         <- niter
      called$verbose       <- verbose
      called$pop_col       <- pop_col
      called$scale         <- scale
      if (!missing(...)) called$dots <- list(...)
      list(ok = TRUE)
    },
    VS_Group_offset    = function(K, y, group_ind, NeighborhoodList, which.prior, niter, verbose, pop_col = NULL, scale = NULL, ...) {
      called$which         <- "VS_Group_offset"
      called$K             <- K
      called$y             <- y
      called$group_ind     <- group_ind
      called$NeighborhoodList <- NeighborhoodList
      called$which.prior   <- which.prior
      called$niter         <- niter
      called$verbose       <- verbose
      called$pop_col       <- pop_col
      called$scale         <- scale
      if (!missing(...)) called$dots <- list(...)
      list(ok = TRUE)
    },
    {
      # 1) standard, no offset
      reset_called()
      silently(MUSC(X = X, y = y, NeighborhoodList = nb, group_ind = NULL, niter = 5, which.prior = "HS", verbose = FALSE))
      expect_equal(called$which, "VS_Standard")
      expect_identical(called$NeighborhoodList, nb)
      expect_equal(ncol(called$K), ncol(X) + 1)
      expect_true(all(called$K[, 1] == 1))

      # 2) standard + offset
      reset_called()
      silently(MUSC(X = X, y = y, NeighborhoodList = nb, group_ind = NULL, niter = 5, which.prior = "SS",
                    pop_col = pop, scale = 100000, verbose = FALSE))
      expect_equal(called$which, "VS_Standard_offset")
      expect_equal(called$which.prior, "SS")
      expect_identical(called$pop_col, pop)
      expect_identical(called$scale, 100000)

      # 3) group, no offset
      reset_called()
      silently(MUSC(X = X, y = y, NeighborhoodList = nb, group_ind = grp, niter = 5, which.prior = "HS", verbose = FALSE))
      expect_equal(called$which, "VS_Group")
      expect_identical(called$group_ind, grp)

      # 4) group + offset
      reset_called()
      silently(MUSC(X = X, y = y, NeighborhoodList = nb, group_ind = grp, niter = 5, which.prior = "HS",
                    pop_col = pop, scale = 100000, verbose = FALSE))
      expect_equal(called$which, "VS_Group_offset")
      expect_identical(called$scale, 100000)
    }
  )
})


# --- Coercions: y matrix first column; X vector -> matrix -----------
test_that("coercions: y matrix first col; X vector -> matrix; K has intercept", {
  nb <- spdep::cell2nb(2, 3, type="queen"); n <- length(nb)
  y_mat <- cbind(1:n, 999)                 # second col should be ignored
  X_vec <- rnorm(n)                        # becomes column matrix
  grp   <- rep(1L, 1)                      # P=1

  called <- new.env(); called$K <- NULL; called$y <- NULL
  testthat::with_mocked_bindings(
    VS_Group = function(K, y, ...) { called$K <<- K; called$y <<- y; list(ok=TRUE) },
    {
      silently(MUSC(X=X_vec, y=y_mat, NeighborhoodList=nb, group_ind=grp, niter=3, verbose=FALSE))
      expect_true(is.matrix(called$K))
      expect_equal(ncol(called$K), 2)      # intercept + 1 feature
      expect_true(all(called$K[,1]==1))
      expect_equal(called$y, y_mat[,1])    # used first column only
    }
  )
})

# --- which.prior matching & errors ----------------------------------
test_that("which.prior must be HS or SS", {
  nb <- spdep::cell2nb(2,2, type="queen")
  X  <- matrix(0, nrow=length(nb), ncol=1)
  y  <- rep(1L, length(nb))
  expect_error(
    silently(MUSC(X=X, y=y, NeighborhoodList=nb, group_ind=NULL, which.prior="BAD", niter=1, verbose=FALSE)),
    "'arg' should be one of",
    fixed = TRUE
  )
  expect_silent(silently(MUSC(X=X, y=y, NeighborhoodList=nb, group_ind=NULL, which.prior="HS", niter=1, verbose=FALSE)))
  expect_silent(silently(MUSC(X=X, y=y, NeighborhoodList=nb, group_ind=NULL, which.prior="SS", niter=1, verbose=FALSE)))
})

# --- Warnings only when verbose ------------------------------------
test_that("count-data warnings only when verbose", {
  nb <- spdep::cell2nb(2,2, type="queen"); n <- length(nb)
  X  <- matrix(0, nrow=n, ncol=1)
  y_nonint <- rep(1.25, n)
  y_unique <- as.numeric(seq_len(n))

  testthat::with_mocked_bindings(
    VS_Group           = function(...) list(mock="VS_Group"),
    VS_Standard        = function(...) list(mock="VS_Standard"),
    VS_Group_offset    = function(...) list(mock="VS_Group_offset"),
    VS_Standard_offset = function(...) list(mock="VS_Standard_offset"),
    {
      expect_warning(silently(MUSC(X=X, y=y_nonint, NeighborhoodList=nb, group_ind=NULL, niter=1, verbose=TRUE)),
                     "nonnegative integers for counts", fixed = FALSE)
      expect_warning(silently(MUSC(X=X, y=y_unique, NeighborhoodList=nb, group_ind=NULL, niter=1, verbose=TRUE)),
                     "possibly continuous", fixed = FALSE)

      expect_silent(silently(MUSC(X=X, y=y_nonint, NeighborhoodList=nb, group_ind=NULL, niter=1, verbose=FALSE)))
      expect_silent(silently(MUSC(X=X, y=y_unique, NeighborhoodList=nb, group_ind=NULL, niter=1, verbose=FALSE)))
    }
  )
})

# --- Optional: niter sanity (add validation in MUSC, then test) -----
# If you add something like:
#   if (!is.numeric(niter) || length(niter)!=1L || niter<1 || !is.finite(niter)) stop("niter must be a positive scalar.")
# then test:
test_that("niter must be a positive scalar (if validated)", {
  nb <- spdep::cell2nb(2,2, type="queen"); n <- length(nb)
  X  <- matrix(0, nrow=n, ncol=1); y <- rep(1L, n)
  expect_error(silently(MUSC(X=X, y=y, NeighborhoodList=nb, group_ind=NULL, niter=0, verbose=FALSE)),
               "niter.*positive", ignore.case = TRUE)
  expect_error(silently(MUSC(X=X, y=y, NeighborhoodList=nb, group_ind=NULL, niter=c(10,20), verbose=FALSE)),
               "niter", ignore.case = TRUE)
})
