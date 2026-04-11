# make_radar_chart_loc.R — REFACTORED
#
# Creates the final summary radar chart. Split into compute_radar() + plot_radar().

# ============================================================================
# compute_radar: Pure computation — no side effects
# ============================================================================
# Returns a list with:
#   $radar_plot_mat — matrix of radar chart values (rows = sig-dataset combos, cols = metrics)
#   $areas          — area ratios for each sig-dataset combination
#   $legend_labels  — formatted legend labels with area ratios
#   $radarplot_rownames — row names for the output table
compute_radar <- function(radar_plot_values, names_sigs, names_datasets) {
  radar_plot_mat <- c()

  all_metrics <- c('sd_median_ratio', 'abs_skewness_ratio', 'prop_top_10_percent',
                   'prop_top_25_percent', 'prop_top_50_percent', 'coeff_of_var_ratio',
                   'med_prop_na', 'med_prop_above_med', 'autocor_median',
                   'rho_mean_med', 'rho_pca1_med', 'rho_mean_pca1',
                   'prop_pca1_var', 'standardization_comp')

  # Fill in missing metrics with 0
  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      diff_names <- base::setdiff(all_metrics,
                                   names(radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]))
      if (length(diff_names) > 0) {
        for (j in seq_along(diff_names)) {
          radar_plot_values[[names_sigs[k]]][[names_datasets[i]]][diff_names[j]] <- 0
        }
      }
    }
  }

  # Flatten nested list into matrix
  for (k in seq_along(names_sigs)) {
    lapply(radar_plot_values[[names_sigs[k]]],
           function(x) radar_plot_mat <<- rbind(radar_plot_mat, x))
  }

  # Add max/min rows for fmsb::radarchart
  radar_plot_mat <- rbind(rep(1, length(radar_plot_values[[1]][[1]])),
                          rep(0, length(radar_plot_values[[1]][[1]])),
                          radar_plot_mat)
  radar_plot_mat <- abs(radar_plot_mat)

  # Compute area ratios
  areas <- c()
  for (i in 3:dim(radar_plot_mat)[1]) {
    areas <- c(areas, sum(sapply(seq_along(radar_plot_mat[i, ]), function(x) {
      if (x < length(radar_plot_mat[i, ])) {
        radar_plot_mat[i, x] * radar_plot_mat[i, x + 1]
      } else {
        radar_plot_mat[i, x] * radar_plot_mat[i, 1]
      }
    })))
  }
  areas <- areas / dim(radar_plot_mat)[2]

  # Build legend labels with area ratios
  legend_labels <- c()
  count <- 1
  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      legend_labels <- c(legend_labels,
                         paste0(names_datasets[i], ' ', names_sigs[k],
                                ' (', format(areas[count], digits = 2), ')'))
      count <- count + 1
    }
  }

  # Build row names for output table
  radarplot_rownames <- c()
  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      radarplot_rownames <- c(radarplot_rownames,
                               paste0(gsub(' ', '.', names_datasets[i]), '_',
                                      gsub(' ', '.', names_sigs[k])))
    }
  }

  # Strip max/min rows for the output matrix
  out_mat <- radar_plot_mat[3:(dim(radar_plot_mat)[1]), ]
  out_mat <- as.matrix(out_mat)
  if (length(radarplot_rownames) == 1) {
    new_colnames <- rownames(out_mat)
    dim(out_mat) <- c(1, length(out_mat))
    colnames(out_mat) <- new_colnames
  }
  row.names(out_mat) <- radarplot_rownames

  list(
    radar_plot_mat = radar_plot_mat,
    output_table = out_mat,
    areas = areas,
    legend_labels = legend_labels,
    radarplot_rownames = radarplot_rownames
  )
}

# ============================================================================
# plot_radar: Side effects — PDF + table output
# ============================================================================
plot_radar <- function(compute_result, names_sigs, names_datasets,
                       out_dir, showResults = FALSE) {
  radar_plot_mat <- compute_result$radar_plot_mat
  areas <- compute_result$areas
  legend_labels <- compute_result$legend_labels

  colours_array <- grDevices::rainbow(length(names_datasets))
  legend_cols <- c()
  legend_lty <- c()
  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      legend_cols <- c(legend_cols, colours_array[i])
      legend_lty <- c(legend_lty, k)
    }
  }

  max_title_length <- max(nchar(outer(names_datasets, names_sigs, paste)))

  # Calculate legend box size
  graphics::plot.new()
  l <- graphics::legend(0.05, 0, legend = legend_labels, seg.len = 2, title = "Datasets",
                        lty = legend_lty, bty = "n", lwd = 1, col = legend_cols,
                        cex = min(0.8, 3 * 10 / max_title_length), plot = FALSE)

  radius <- 1.5
  if (graphics::grconvertX(l$rect$h, from = "ndc", to = "device") > radius) {
    x_dist <- radius
  } else {
    x_dist <- sqrt((radius^2) - ((radius - graphics::grconvertX(l$rect$h, from = "ndc", to = "device"))^2))
  }
  padding <- (graphics::grconvertX(x_dist, from = "user", to = "device") +
               graphics::grconvertX(l$rect$w, from = "ndc", to = "device")) -
    graphics::grconvertX(radius, from = "user", to = "device")
  padding <- padding * (padding > 0)
  padding <- graphics::grconvertX(padding, from = 'device', to = "inches")

  if (showResults) {
    grDevices::dev.new()
  } else {
    grDevices::pdf(file.path(out_dir, 'sig_radarplot.pdf'), width = 10, height = 10)
  }
  orig_par <- graphics::par('mai')
  graphics::par(xpd = TRUE, mai = (graphics::par('mai') + c(0, 0, 0, padding)))

  fmsb::radarchart(as.data.frame(radar_plot_mat),
                   maxmin = TRUE, axistype = 1,
                   cglcol = 'grey', axislabcol = 'black',
                   caxislabels = seq(0, 1, length.out = 5),
                   cglty = 1, cglwd = 1, calcex = 0.5,
                   vlabels = c('Relative Med. SD', 'Skewness',
                               expression(sigma["" >= "10%"]), expression(sigma["" >= "25%"]),
                               expression(sigma["" >= "50%"]), 'Coef. of Var.',
                               'Non-NA Prop.', 'Prop. Expressed',
                               'Intra-sig. Corr.', expression(rho["Mean,Med"]),
                               expression(rho["PCA1,Med"]), expression(rho["Mean,PCA1"]),
                               expression(sigma["PCA1"]), expression(rho["Med,Z-Med"])),
                   vlcex = 0.6, title = 'Signature Summary',
                   pty = 16, plty = legend_lty, pcol = legend_cols, plwd = 2)

  sorted_idx <- order(-areas)
  graphics::legend(x_dist, 1.25, xpd = NA,
                   legend = legend_labels[sorted_idx], seg.len = 2, title = "Datasets",
                   lty = legend_lty[sorted_idx], bty = "n", lwd = 1,
                   col = legend_cols[sorted_idx],
                   cex = min(0.8, 3 * 10 / max_title_length))

  if (showResults) {
    grDevices::dev.copy(grDevices::pdf, file.path(out_dir, 'sig_radarplot.pdf'),
                        width = 10, height = 10)
  }
  if (grDevices::dev.cur() != 1) grDevices::dev.off()

  # Write radar chart table
  if (!dir.exists(file.path(out_dir, 'radarchart_table'))) {
    dir.create(file.path(out_dir, 'radarchart_table'))
  }
  utils::write.table(compute_result$output_table,
                     file = file.path(out_dir, 'radarchart_table', 'radarchart_table.txt'),
                     quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE)

  graphics::par(mar = orig_par)
}

# ============================================================================
# make_radar_chart_loc: Original API wrapper — backward-compatible
# ============================================================================
make_radar_chart_loc <- function(radar_plot_values, showResults = FALSE, names_sigs,
                                  names_datasets, out_dir = '~', file) {

  result <- compute_radar(radar_plot_values, names_sigs, names_datasets)
  plot_radar(result, names_sigs, names_datasets, out_dir, showResults)

  cat('Radar chart made successfully.\n', file = file)
}
