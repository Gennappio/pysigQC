# grab_grob.R
#
# Captures the current R graphics device content as a grid grob (graphical object).
# Used by eval_struct_loc to capture individual biclustering heatmaps into a list
# that can later be arranged in a grid layout by draw.heatmaps().
# Relies on gridGraphics::grid.echo() to convert base graphics to grid graphics.
# @keywords grab_grob

grab_grob <- function(){
  # require(gridGraphics)
  # require(grid)
  # In non-interactive Rscript runs gridGraphics::grid.echo() emits the warning
  # "No graphics to replay" and grid::grid.grab() then returns NULL. A NULL grob
  # propagates into draw.heatmaps() -> grid::editGrob(NULL, ...) which raises
  # 'no applicable method for validGrob applied to NULL', aborting eval_*_loc()
  # before the radar values are stored. Coerce any failure mode to a nullGrob()
  # so the surrounding lapply finishes and the caller's metric assignments run.
  g <- tryCatch({
    gridGraphics::grid.echo()
    grid::grid.grab()
  }, error = function(e) NULL)
  if (is.null(g) || !grid::is.grob(g)) grid::nullGrob() else g
}
