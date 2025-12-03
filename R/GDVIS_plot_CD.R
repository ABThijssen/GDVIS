#' Plot GDVIS cross-disorder plot
#'
#' This function plots a GDVIS visualization based on GDVIS parameters.

#' @param input_triangle_parameters the path to the RData GDVIS output, which ends with CD.triangle_parameters.RData
#' @param show.rendering option to show the rendering of the plot in an external window, default is TRUE
#' @param show.names option to show the names of the groups, default is TRUE
#' @param webversion needs to stay on FALSE
#' @export


# Function to plot 3D triangles
GDVIS_plot_CD <- function(input_triangle_parameters, show.rendering = TRUE, show.names = TRUE, webversion = FALSE) {



#### Set up the environment and load data -----------------------------------------------------------



  # Create a new environment
  temp_triangle_env <- new.env()

  # Load the parameters file into this new environment
  load(input_triangle_parameters, envir = temp_triangle_env)

  # Extract the list from the environment
  triangle_list <- temp_triangle_env$CD.triangle.output.list

  # Show message
  cli::cli_h1(paste0("Running ", cli::col_cyan("GDVIS plot"), " on ", cli::col_cyan(triangle_list$name_trait1), " with ", cli::col_cyan(triangle_list$name_trait2), " with ", cli::col_cyan(triangle_list$name_trait2)))

  # Directly access the objects
  with(temp_triangle_env, {

    # Extract the variables of the list
    list2env(CD.triangle.output.list, envir = temp_triangle_env)



#### Set up for plot CD triangle --------------------------------------------------------------------------------


    coord_trait1 <- c(x.trait1, y.trait1, z.trait1)
    coord_trait2 <- c(x.trait2, y.trait2, z.trait2)
    coord_trait3 <- c(x.trait3, y.trait3, z.trait3)
    coord_trait1.con <- c(x.trait1.con, y.trait1.con, z.trait1.con)
    coord_trait2.con <- c(x.trait2.con, y.trait2.con, z.trait2.con)
    coord_trait3.con <- c(x.trait3.con, y.trait3.con, z.trait3.con)

# Calculate midpoints for arrow
coord_arrow.trait1_line <- (coord_trait1 - coord_trait1.con) / 2
coord_arrow.trait2_line <- (coord_trait2 - coord_trait2.con) / 2
coord_arrow.trait3_line <-  (coord_trait3 - coord_trait3.con) / 2

# Calculate dir vectors of line
dir.vect_trait1_line <- coord_trait1 - coord_trait1.con
dir.vect_trait2_line <- coord_trait2 - coord_trait2.con
dir.vect_trait3_line <- coord_trait3 - coord_trait3.con

# Function to create arrowheads
add_arrowhead <- function(fig, dirvect, midpoint,color) {

  # Get rotation matrix
  cos30 <- cos(pi/6)  # cos(30 degrees)
  sin30 <- sin(pi/6)  # sin(30 degrees)
  rotation_matrix_right <- matrix(c(cos30, -sin30, 0,  sin30, cos30, 0,   0, 0, 1  ), nrow = 3, byrow = TRUE)
  rotation_matrix_left <- matrix(c( cos30, sin30, 0,-sin30, cos30, 0,   0, 0, 1   ), nrow = 3, byrow = TRUE)
  rotation_matrix_up <- matrix(c(1, 0, 0, 0, cos30, -sin30,0, sin30,  cos30 ), nrow = 3, byrow = TRUE)
  rotation_matrix_down <- matrix(c( 1, 0, 0,0, cos30, sin30, 0, -sin30, cos30), nrow = 3, byrow = TRUE)

  # Get rotated dir vectors
  rotated_dirvect_right <- rotation_matrix_right %*% dirvect
  rotated_dirvect_left <- rotation_matrix_left %*% dirvect
  rotated_dirvect_up <- rotation_matrix_up %*% dirvect
  rotated_dirvect_down <- rotation_matrix_down %*% dirvect
  # Normalize the dirvectors
  dir_vector_right_normalized <- rotated_dirvect_right / sqrt(sum(rotated_dirvect_right^2))
  dir_vector_left_normalized <- rotated_dirvect_left / sqrt(sum(rotated_dirvect_left^2))
  dir_vector_up_normalized <- rotated_dirvect_up / sqrt(sum(rotated_dirvect_up^2))
  dir_vector_down_normalized <- rotated_dirvect_down / sqrt(sum(rotated_dirvect_down^2))
  # Get new length of the dirvectors
  scaled_vector_right <- dir_vector_right_normalized * 0.0075
  scaled_vector_left <- dir_vector_left_normalized * 0.0075
  scaled_vector_up <- dir_vector_up_normalized * 0.0075
  scaled_vector_down <- dir_vector_down_normalized * 0.0075
  # Get coordinates
  arrow_start_right <- midpoint - scaled_vector_right
  arrow_start_left <- midpoint - scaled_vector_left
  arrow_start_up <- midpoint - scaled_vector_up
  arrow_start_down <- midpoint - scaled_vector_down

  # Adding the arrowhead lines to the plot
  fig <- fig %>%
    plotly::add_trace(x = c(midpoint[1], arrow_start_right[1]),
              y = c(midpoint[2], arrow_start_right[2]),
              z = c(midpoint[3], arrow_start_right[3]),
              type = 'scatter3d', mode = 'lines', line = list(color = color, width = 3), showlegend = FALSE, hoverinfo = 'none') %>%
    plotly::add_trace(x = c(midpoint[1], arrow_start_left[1]),
              y = c(midpoint[2], arrow_start_left[2]),
              z = c(midpoint[3], arrow_start_left[3]),
              type = 'scatter3d', mode = 'lines', line = list(color = color, width = 3), showlegend = FALSE, hoverinfo = 'none') %>%
    plotly::add_trace(x = c(midpoint[1], arrow_start_up[1]),
              y = c(midpoint[2], arrow_start_up[2]),
              z = c(midpoint[3], arrow_start_up[3]),
              type = 'scatter3d', mode = 'lines', line = list(color = color, width = 3), showlegend = FALSE, hoverinfo = 'none') %>%
    plotly::add_trace(x = c(midpoint[1], arrow_start_down[1]),
              y = c(midpoint[2], arrow_start_down[2]),
              z = c(midpoint[3], arrow_start_down[3]),
              type = 'scatter3d', mode = 'lines', line = list(color = color, width = 3), showlegend = FALSE, hoverinfo = 'none') %>%
    plotly::add_trace(x = c(arrow_start_left[1], arrow_start_right[1]),
              y = c(arrow_start_left[2], arrow_start_right[2]),
              z = c(arrow_start_left[3], arrow_start_right[3]),
              type = 'scatter3d', mode = 'lines', line = list(color = color, width = 3), showlegend = FALSE, hoverinfo = 'none') %>%
    plotly::add_trace(x = c(arrow_start_up[1], arrow_start_down[1]),
              y = c(arrow_start_up[2], arrow_start_down[2]),
              z = c(arrow_start_up[3], arrow_start_down[3]),
              type = 'scatter3d', mode = 'lines', line = list(color = color, width = 3), showlegend = FALSE, hoverinfo = 'none') %>%
    plotly::add_trace(x = c(arrow_start_right[1], arrow_start_left[1], midpoint[1]),
              y = c(arrow_start_right[2], arrow_start_left[2], midpoint[2]),
              z = c(arrow_start_right[3], arrow_start_left[3], midpoint[3]),
              type = 'mesh3d',
              i = c(0), j = c(1), k = c(2),
              vertexcolor = matrix(rep(col2rgb(color)/255, 3), nrow = 3, byrow = TRUE),  # Apply the exact color to each vertex
              opacity = 1,
              showlegend = FALSE, showscale = FALSE,
              hoverinfo = 'none') %>%
    plotly::add_trace(x = c(arrow_start_up[1], arrow_start_down[1], midpoint[1]),
              y = c(arrow_start_up[2], arrow_start_down[2], midpoint[2]),
              z = c(arrow_start_up[3], arrow_start_down[3], midpoint[3]),
              type = 'mesh3d',
              i = c(0), j = c(1), k = c(2),
              vertexcolor = matrix(rep(col2rgb(color)/255, 3), nrow = 3, byrow = TRUE),  # Apply the exact color to each vertex
              opacity = 1,
              showlegend = FALSE, showscale = FALSE,
              hoverinfo = 'none')

  return(fig)    }



#### Create CD triangle --------------------------------------------------------------------------------


# Initialize a 3D scatter plot
fig <- plotly::plot_ly() %>%

  # trait1 line
  plotly::add_trace(x = c(x.trait1.con, x.trait1), y = c(y.trait1.con, y.trait1), z = c(z.trait1.con, z.trait1),
            type = 'scatter3d', mode = 'lines',
            line = list(color = "#083F80", width = 6),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # trait2 line
  plotly::add_trace(x = c(x.trait2.con, x.trait2), y = c(y.trait2.con, y.trait2), z = c(z.trait2.con, z.trait2),
            type = 'scatter3d', mode = 'lines',
            line = list(color = "#3696D2", width = 6),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # trait3 line
  plotly::add_trace(x = c(x.trait3.con, x.trait3), y = c(y.trait3.con, y.trait3), z = c(z.trait3.con, z.trait3),
            type = 'scatter3d', mode = 'lines',
            line = list(color = "#B4D8F0", width = 6),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # trait1-trait2 line
  plotly::add_trace(x = c(x.trait1, x.trait2), y = c(y.trait1, y.trait2), z = c(z.trait1, z.trait2),
            type = 'scatter3d', mode = 'lines',
            line = list(color = "#ED7D31", width = 6),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # trait1-trait3 line
  plotly::add_trace(x = c(x.trait1, x.trait3), y = c(y.trait1, y.trait3), z = c(z.trait1, z.trait3),
            type = 'scatter3d', mode = 'lines',
            line = list(color = "#00B050", width = 6),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # trait2-trait3 line
  plotly::add_trace(x = c(x.trait2, x.trait3), y = c(y.trait2, y.trait3), z = c(z.trait2, z.trait3),
            type = 'scatter3d', mode = 'lines',
            line = list(color = "#7030A0", width = 6),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # # trait1.con line to trait2.con
  # plotly::add_trace(x = c(x.trait1.con, x.trait2.con), y = c(y.trait1.con, y.trait2.con), z = c(z.trait1.con, z.trait2.con),
  #           type = 'scatter3d', mode = 'lines',
  #           line = list(color = "black", width = 2),
  #           showlegend = FALSE,
  #           hoverinfo = 'none') %>%
  #
  # # trait1.con line to trait3.con
  # plotly::add_trace(x = c(x.trait1.con, x.trait3.con), y = c(y.trait1.con, y.trait3.con), z = c(z.trait1.con, z.trait3.con),
  #           type = 'scatter3d', mode = 'lines',
  #           line = list(color = "black", width = 2),
  #           showlegend = FALSE,
  #           hoverinfo = 'none') %>%
  #
  # # trait2.con line to trait3.con
  # plotly::add_trace(x = c(x.trait2.con, x.trait3.con), y = c(y.trait2.con, y.trait3.con), z = c(z.trait2.con, z.trait3.con),
  #           type = 'scatter3d', mode = 'lines',
  #           line = list(color = "black", width = 2),
  #           showlegend = FALSE,
  #           hoverinfo = 'none') %>%

  # Add popmean
  plotly::add_trace(x = x.popmean, y = y.popmean, z = z.popmean,
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 5, color = 'black', symbol = "cross"),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # Add trait1.con point
  plotly::add_trace(x = x.trait1.con,
            y = y.trait1.con,
            z = z.trait1.con,
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 4, color = 'black', line = list(color = "black", width = 3), symbol = "circle"),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # Add trait1 point
  plotly::add_trace(x = x.trait1,
            y = y.trait1,
            z = z.trait1,
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 4, color = 'black', line = list(color = "black", width = 3), symbol = "circle"),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # Add trait2.con point
  plotly::add_trace(x = x.trait2.con,
            y = y.trait2.con,
            z = z.trait2.con,
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 4, color = 'black', line = list(color = "black", width = 3), symbol = "circle"),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # Add trait2 point
  plotly::add_trace(x = x.trait2,
            y = y.trait2,
            z = z.trait2,
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 4, color = 'black', line = list(color = "black", width = 3), symbol = "circle"),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # Add trait3.con point
  plotly::add_trace(x = x.trait3.con,
            y = y.trait3.con,
            z = z.trait3.con,
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 4, color = 'black', line = list(color = "black", width = 3), symbol = "circle"),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # Add trait3 point
  plotly::add_trace(x = x.trait3,
            y = y.trait3,
            z = z.trait3,
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 4, color = 'black', line = list(color = "black", width = 3), symbol = "circle"),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # Clean up plot
  plotly::layout(scene = list(
    xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, showspikes = FALSE, title = ''),
    yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, showspikes = FALSE, title = ''),
    zaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, showspikes = FALSE, title = ''),
    aspectmode = "data"),
    hovermode = 'none') %>%

  # Add text annotations
  plotly::add_trace(x = c(x.trait1, x.trait2, x.trait3, x.trait1.con, x.trait2.con, x.trait3.con),
            y = c(y.trait1, y.trait2, y.trait3, y.trait1.con, y.trait2.con, y.trait3.con),
            z = c(z.trait1, z.trait2, z.trait3, z.trait1.con, z.trait2.con, z.trait3.con),
            type = "scatter3d", mode = "text",
            text = c(name_trait1, name_trait2, name_trait3, name_trait1.con, name_trait2.con, name_trait3.con),
            textposition = c("bottom right", "bottom right", "bottom right", "bottom right", "bottom right", "bottom right"),
            textfont = list(size = 12),
            showlegend = FALSE,
            hoverinfo = 'none') %>%

  # Make the cases triangle gray
  plotly::add_trace(x = c(x.trait1, x.trait2, x.trait3),
            y = c(y.trait1, y.trait2, y.trait3),
            z = c(z.trait1, z.trait2, z.trait3),
            type = 'mesh3d',
            i = c(0),  j = c(1), k = c(2),
            colorscale = list(c(0, "lightgray"), c(1, "lightgray")), # Define a gray colorscale
            intensity = rep(1, 3), # Use intensity to apply the colorscale uniformly
            opacity = 0.5,
            showlegend = FALSE, showscale = FALSE,
            hoverinfo = 'none')

# Add arrowhead to the plot
fig <- add_arrowhead(fig, dir.vect_trait1_line, coord_arrow.trait1_line, color = "#083F80")
fig <- add_arrowhead(fig, dir.vect_trait2_line, coord_arrow.trait2_line, color = "#3696D2")
fig <- add_arrowhead(fig, dir.vect_trait3_line, coord_arrow.trait3_line, color = "#B4D8F0")

if (webversion == FALSE) {
  # Save
  save_path <- paste0(folder_location, "/", name_trait1, "_",name_trait2,"_" ,name_trait3, "_CD.plot.html") %>% gsub("^/", "", .)
  htmlwidgets::saveWidget(fig, save_path)

  # Return message
  cli::cli_alert_info(paste0(" Plot saved as ", folder_location,"/",name_trait1, "_",name_trait2,"_" ,name_trait3, "_CD.plot.html"))
}



#### Create legend --------------------------------------------------------------------------------


### Round values
round_three_decimals <- function(objects) {
  for (i in objects) {assign(i, formatC(get(i,envir = temp_triangle_env), digits = 3, format = "f"),envir = temp_triangle_env)  }}
round_no_decimals <- function(objects) {
  for (i in objects) {assign(i, formatC(get(i,envir = temp_triangle_env), digits = 0, format = "f"),envir = temp_triangle_env)  }}
round_three_decimals(c("h2_trait1", "h2_trait2", "h2_trait3", "h2_trait1.trait2", "h2_trait1.trait3", "h2_trait2.trait3", "rg_trait1.cases_trait2.cases_vs_trait1.cases_trait3.cases","rg_trait1.cases_trait2.cases_vs_trait2.cases_trait3.cases","rg_trait1.cases_trait3.cases_vs_trait2.cases_trait3.cases"))
round_no_decimals(c("a.deg_trait1_trait2", "a.deg_trait1_trait3", "a.deg_trait2_trait3", "a.deg_trait1.cases_trait2.cases_vs_trait1.cases_trait3.cases","a.deg_trait1.cases_trait2.cases_vs_trait2.cases_trait3.cases", "a.deg_trait1.cases_trait3.cases_vs_trait2.cases_trait3.cases"))

# Make legend
legend <- ggplot() +

  # Plot pop mean
  geom_point(ggplot2::aes(x = 0.1, y = -0.125), shape = 3, size = 4, stroke = 2) +
  geom_text(ggplot2::aes(x = 0.2, y = -0.125), label = "Population mean", color = "black", hjust = 0, size = 5) +

  # Plot angles with degrees
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(0.125, 0.095)), color = "#00B050", linewidth = 2) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(0.125, 0.175)), color = "#ED7D31", linewidth = 2) +
  geom_text(ggplot2::aes(x = 0.3, y = 0.14), label = paste0(rg_trait1.cases_trait2.cases_vs_trait1.cases_trait3.cases, " (", a.deg_trait1.cases_trait2.cases_vs_trait1.cases_trait3.cases, "°)"), color = "black", hjust = 0, size = 6) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(0.26, 0.23)), color = "#7030A0", linewidth = 2) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(0.26, 0.31)), color = "#ED7D31", linewidth = 2) +
  geom_text(ggplot2::aes(x = 0.3, y = 0.27), label = paste0(rg_trait1.cases_trait2.cases_vs_trait2.cases_trait3.cases, " (", a.deg_trait1.cases_trait2.cases_vs_trait2.cases_trait3.cases, "°)"), color = "black", hjust = 0, size = 6) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(0.385, 0.355)), color = "#7030A0", linewidth = 2) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(0.385, 0.435)), color = "#00B050", linewidth = 2) +
  geom_text(ggplot2::aes(x = 0.3, y = 0.395), label = paste0(rg_trait1.cases_trait3.cases_vs_trait2.cases_trait3.cases, " (", a.deg_trait1.cases_trait3.cases_vs_trait2.cases_trait3.cases, "°)"), color = "black", hjust = 0, size = 6) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(0.51, 0.48)), color = "#083F80", linewidth = 2) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(0.51, 0.56)), color = "#3696D2", linewidth = 2) +
  geom_text(ggplot2::aes(x = 0.3, y = 0.52), label = paste0(rg_trait1_trait2, " (", a.deg_trait1_trait2, "°)"), color = "black", hjust = 0, size = 6) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(0.635, 0.605)), color = "#083F80", linewidth = 2) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(0.635, 0.685)), color = "#B4D8F0", linewidth = 2) +
  geom_text(ggplot2::aes(x = 0.3, y = 0.645), label = paste0(rg_trait1_trait3, " (", a.deg_trait1_trait3, "°)"), color = "black", hjust = 0, size = 6) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(0.765, 0.735)), color = "#B4D8F0", linewidth = 2) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(0.765, 0.815)), color = "#3696D2", linewidth = 2) +
  geom_text(ggplot2::aes(x = 0.3, y = 0.775), label = paste0(rg_trait2_trait3, " (", a.deg_trait2_trait3, "°)"), color = "black", hjust = 0, size = 6) +
  annotate("text", x = 0, y = 0.925, label = "italic(r)[g] ~ '(degrees)'", color = "black", hjust = 0, size = 6.5, parse = TRUE) +


  # plot h2 lines
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(1.25, 1.25)), color = "#7030A0", linewidth = 2) +
  geom_text(ggplot2::aes(x = 0.3, y = 1.26), label = h2_trait2.trait3, color = "black", hjust = 0, size = 6) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(1.375, 1.375)), color = "#00B050", linewidth = 2) +
  geom_text(ggplot2::aes(x = 0.3, y = 1.385), label = h2_trait1.trait3, color = "black", hjust = 0, size = 6) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(1.5, 1.5)), color = "#ED7D31", linewidth = 2) +
  geom_text(ggplot2::aes(x = 0.3, y = 1.51), label = h2_trait1.trait2, color = "black", hjust = 0, size = 6) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(1.625, 1.625)), color = "#B4D8F0", linewidth = 2) +
  geom_text(ggplot2::aes(x = 0.3, y = 1.635), label = h2_trait3, color = "black", hjust = 0, size = 6) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(1.75, 1.75)), color = "#3696D2", linewidth = 2) +
  geom_text(ggplot2::aes(x = 0.3, y = 1.76), label = h2_trait2, color = "black", hjust = 0, size = 6) +
  geom_line(ggplot2::aes(x = c(0.1, 0.25), y = c(1.875, 1.875)), color = "#083F80", linewidth = 2) +
  geom_text(ggplot2::aes(x = 0.3, y = 1.885), label = h2_trait1, color = "black", hjust = 0, size = 6) +

  annotate("text", x = 0, y = 2.035, label = "Heritability", color = "black", hjust = 0, size = 6.5) +

  # Clean background
  theme_minimal() +
  theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),   # Remove axes titles
        axis.text.x = ggplot2::element_blank(),  axis.text.y = ggplot2::element_blank(),    # Remove axes text
        axis.ticks = ggplot2::element_blank(),                                     # Remove axes ticks
        panel.grid.major = ggplot2::element_blank(),                               # Remove major gridlines
        panel.grid.minor = ggplot2::element_blank(),                               # Remove minor gridlines
        panel.background = ggplot2::element_rect(fill = "white", colour = NA)) +   # Set background color to white
  xlim(0, 1.6) +
  # Scale axes
  coord_fixed(ratio = 1)


### Save legend
if(webversion == TRUE) {
  suppressWarnings({ ggsave(paste0(folder_location, "/", name_trait1, "_",name_trait2,"_" ,name_trait3,"_CD.legend.png"), plot = legend, width = 4, height = 8, dpi = 300) }) }

# Return message
cli::cli_alert_info(paste0(" Legend saved as ", folder_location,"/",name_trait1, "_",name_trait2,"_" ,name_trait3,"_CD.legend.png"))
cli::cli_alert_success("GDVIS plot succesfully finished")
if (show.rendering == TRUE & webversion == FALSE) {  cli::cli_alert_info("Showing plot in external window ") }


#### Make shiny object to have plot and legend together --------------------------------------------------------------------------------


if (webversion == FALSE) {
# Define UI for the shiny app
ui <- shiny::fluidPage(fluidRow(
  shiny::column(8, plotly::plotlyOutput("plotly_plot", height = "600px")),   # 8/12 width for the plotly plot
  shiny::column(4, shiny::plotOutput("ggplot_legend", height = "600px"))))    # 4/12 width for the ggplot legend

# Define server logic for the shiny app
server <- function(input, output, session) {

  # Render the plotly plot (fig)
  output$plotly_plot <- plotly::renderPlotly({  fig  })

  # Render the ggplot legend
  output$ggplot_legend <-  suppressWarnings({ shiny::renderPlot({  legend  })   })  }

# Run the shiny application
shiny_app <- shiny::shinyApp(ui = ui, server = server)

# Remove the folder if it already already exists (shiny doesn't overwrite)
save_path <- paste0(folder_location, "/",name_trait1, "_",name_trait2,"_" ,name_trait3, "_CD.plot_files")
if (dir.exists(save_path)) { unlink(save_path, recursive = T)   }

# Save the Shiny app script as a .R file
save_path_shiny <- paste0(folder_location, "/", name_trait1, "_",name_trait2,"_" ,name_trait3, "_shiny_app.rds") %>% gsub("^/", "", .)
saveRDS(shiny_app, file = save_path_shiny)

# Save the UI and server components separately
save_path_ui <- paste0(folder_location, "/", name_trait1, "_",name_trait2,"_" ,name_trait3, "_ui.rds") %>% gsub("^/", "", .)
saveRDS(ui, file = save_path_ui)
save_path_server <- paste0(folder_location, "/", name_trait1, "_",name_trait2,"_" ,name_trait3, "_server.rds")  %>% gsub("^/", "", .)
saveRDS(server, file = save_path_server)

# Show if asked
if (show.rendering == TRUE) { assign("shiny_app", shiny_app, envir = .GlobalEnv); shiny::runApp(get("shiny_app", envir = .GlobalEnv)) }}


#### For website plot --------------------------------------------------------------------------------------------------------

# # Define UI (only the inner UI element you want to inject via renderUI)
# cd_ui_static <- fluidRow(
#   column(8, plotlyOutput("plotly_plot", height = "600px")),
#   column(4, plotOutput("ggplot_legend", height = "600px")))
#
# # Define server (must match exact input/output names used above)
# cd_server_static <- function(input, output, session) {
#   output$plotly_plot <- renderPlotly({ fig })
#   output$ggplot_legend <- renderPlot({
#     req(current_page() == "show_cd_plot")
#     legend
#   })
# }
#
# # Save UI and server for later use in the main app
# saveRDS(cd_ui_static, file = file.path(folder_location, paste0(name_trait1, "_", name_trait2, "_", name_trait3, "_ui.rds")))
# saveRDS(cd_server_static, file = file.path(folder_location, paste0(name_trait1, "_", name_trait2, "_", name_trait3, "_server.rds")))

# For user website plot
# Return plot objects
if (webversion == TRUE) {

  # Return plot objects
  return(list(
    fig = fig,
    legend = legend
  ))
}


  }) # End with env

} # End function
