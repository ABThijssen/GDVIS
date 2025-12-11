#' Plot GDVIS plot
#'
#' This function plots a GDVIS visualization based on GDVIS parameters.

#' @param input_triangle_parameters the path to the RData GDVIS output, which ends with 3D.triangle_parameters.RData
#' @param show.rendering option to show the rendering of the plot in an external window, default is TRUE
#' @param show.names option to show the names of the groups, default is TRUE
#' @param x_upper change axes to make the plots all the same size or make the legend appear outside of the plot (only for 2D plot)
#' @param x_lower change axes to make the plots all the same size or make the legend appear outside of the plot (only for 2D plot)
#' @param y_upper change axes to make the plots all the same size or make the legend appear outside of the plot (only for 2D plot)
#' @param y_lower change axes to make the plots all the same size or make the legend appear outside of the plot (only for 2D plot)
#' @export
# Function to plot GDVIS triangles
GDVIS_plot <- function(input_triangle_parameters, show.rendering = TRUE, show.names = TRUE, x_lower = NULL, x_upper = NULL, y_lower = NULL, y_upper = NULL) {

# Create a new environment
temp_triangle_env <- new.env()

# Set to NULL before overwriting the correct one
temp_triangle_env$triangle.output.list <- NULL
temp_triangle_env$double.triangle.output.list <- NULL
temp_triangle_env$CD.triangle.output.list <- NULL

# Load the parameters file into this new environment
load(input_triangle_parameters, envir = temp_triangle_env)

# Extract the list
if (is.null(temp_triangle_env$triangle.output.list) == TRUE & is.null(temp_triangle_env$double.triangle.output.list) == TRUE & is.null(temp_triangle_env$CD.triangle.output.list) == TRUE) { temp_triangle_env$triangle.output.list }
if (is.null(temp_triangle_env$triangle.output.list) == FALSE) { data.path <- temp_triangle_env$triangle.output.list }
if (is.null(temp_triangle_env$double.triangle.output.list) == FALSE) { data.path <- temp_triangle_env$double.triangle.output.list }
if (is.null(temp_triangle_env$CD.triangle.output.list) == FALSE) { data.path <- temp_triangle_env$CD.triangle.output.list }

# Set all plot types to FALSE
temp_triangle_env$plot_3D <- FALSE
temp_triangle_env$plot_2D.2D <- FALSE
temp_triangle_env$plot_CD <- FALSE

# Extract the data (to check the plot type)
list2env(data.path, envir = temp_triangle_env)

# Redirect to correct plot function

## 2D
if (temp_triangle_env$plot_3D == FALSE & temp_triangle_env$plot_CD == FALSE & temp_triangle_env$plot_2D.2D == FALSE) {
  return(GDVIS::GDVIS_plot_2D(input_triangle_parameters, x_lower = x_lower, x_upper = x_upper, y_lower = y_lower, y_upper = y_upper))
}

## 3D
if (temp_triangle_env$plot_3D == TRUE) {
  return(GDVIS::GDVIS_plot_3D(input_triangle_parameters, show.rendering = show.rendering, show.names = show.names))
}

## 2D.2D
if (temp_triangle_env$plot_2D.2D == TRUE) {
  return(GDVIS::GDVIS_plot_2D_2D(input_triangle_parameters, show.rendering = show.rendering, show.names = show.names))
}

## CD
if (temp_triangle_env$plot_CD == TRUE) {
  return(GDVIS::GDVIS_plot_CD(input_triangle_parameters, show.rendering = show.rendering, show.names = show.names))
}


} # end GDVIS_plot
