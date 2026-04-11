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
  gridGraphics::grid.echo()
  grid::grid.grab()
}
