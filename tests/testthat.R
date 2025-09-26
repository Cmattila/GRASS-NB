# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(MUSC)

test_check("MUSC")

# devtools::check()


# devtools::load_all()
# devtools::test(filter = "MUSC_main") # [ FAIL 0 | WARN 0 | SKIP 0 | PASS 39 ]
# devtools::test(filter = "MUSC_main_All") # [ FAIL 0 | WARN 0 | SKIP 0 | PASS 25 ]
#
# devtools::document()
