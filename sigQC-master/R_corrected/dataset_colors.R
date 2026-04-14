#' Generate N maximally-distinguishable colors for dataset identification.
#'
#' First three colors are anchored at red, blue, green (ColorBrewer Set1).
#' Additional colors are placed in HCL space at hues that maximize the minimum
#' angular distance from all previously placed hues.
#'
#' @param n Number of colors needed (>= 1)
#' @return Character vector of n hex color strings
dataset_colors <- function(n) {
  anchor_hex  <- c("#E41A1C", "#377EB8", "#4DAF4A")  # red, blue, green
  if (n <= 3) return(anchor_hex[seq_len(n)])

  anchor_hues <- c(12, 255, 135)
  placed <- anchor_hues
  for (j in seq_len(n - 3)) {
    candidates <- seq(0, 359, by = 1)
    min_dists <- vapply(candidates, function(h) {
      d <- abs(h - placed)
      min(pmin(d, 360 - d))
    }, numeric(1))
    placed <- c(placed, candidates[which.max(min_dists)])
  }
  extra <- grDevices::hcl(h = placed[4:n], c = 80, l = 60)
  c(anchor_hex, extra)
}
