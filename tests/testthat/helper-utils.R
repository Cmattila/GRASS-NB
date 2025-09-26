# Swallow stdout/stderr from cat()/print() so testthat doesn't flag "produced output".
silently <- function(expr) invisible(capture.output(force(expr)))
