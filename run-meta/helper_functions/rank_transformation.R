rank_within_sample <- function(mat) {
  apply(mat, 2, function(x) rank(x, ties.method = "average"))
}
