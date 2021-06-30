library(devtools)
packageVersion("devtools")
library(tidyverse)
library(fs)
create_package(getwd())
use_git()
dir_info(all = TRUE, regexp = "^[.]git$") %>%
  select(path, type)

(a <- factor(c("character", "hits", "your", "eyeballs")))
(b <- factor(c("but", "integer", "where it", "counts")))
c(a, b)
factor(c(as.character(a), as.character(b)))
fbind <- function(a, b) {
  factor(c(as.character(a), as.character(b)))
}
use_r("fbind")
load_all()
fbind(a, b)
exists("fbind", where = globalenv(), inherits = FALSE)
use_mit_license("Mao Yinan")
document()
check()
install()

use_testthat()
use_test("fbind")
test()

use_package("forcats")
use_r("fcount")
load_all()
fcount(iris$Species)

document()
use_readme_rmd()

check()
install()
