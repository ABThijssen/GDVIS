#' Plot GDVIS 3D plot
#'
#' This function plots a GDVIS visualization based on GDVIS parameters.

#' @param input_triangle_parameters the path to the RData GDVIS output, which ends with 3D.triangle_parameters.RData
#' @param show.rendering option to show the rendering of the plot in an external window, default is TRUE
#' @param show.names option to show the names of the groups, default is TRUE
#' @param webversion needs to stay on FALSE
#' @export
# Function to plot 3D triangles
GDVIS_plot_3D <- function(input_triangle_parameters, show.rendering = TRUE, show.names = TRUE, webversion = FALSE) {



#### Set up the environment and load data -----------------------------------------------------------



    # Create a new environment
    temp_triangle_env <- new.env()

    # Load the parameters file into this new environment
    load(input_triangle_parameters, envir = temp_triangle_env)

    # Extract the list from the environment
    triangle_list <- temp_triangle_env$triangle.output.list

    # Show message
    cli::cli_h1(paste0("Running ", cli::col_cyan("GDVIS plot"), " on ", cli::col_cyan(triangle_list$plot_title), " with ", cli::col_cyan(triangle_list$name_ext)))

    # Directly access the objects
    with(temp_triangle_env, {

    # Extract the variables of the list
    list2env(triangle.output.list, envir = temp_triangle_env)



#### Set up for plot 3D triangle --------------------------------------------------------------------------------



      # Function to create a gradient line with distinct segments for continuous traits
      add_explicit_color_line <- function(fig, x, y, z, colors, width) {
        n_segments <- length(colors) - 1
        for (i in seq(1, n_segments)) {
          fig <- fig %>%
            plotly::add_trace(
              x = x[c(i, i + 1)], y = y[c(i, i + 1)], z = z[c(i, i + 1)],
              type = 'scatter3d', mode = 'lines',
              line = list(color = colors[i], width = width),
              showlegend = FALSE,
              hoverinfo = 'none'
            )
        }
        return(fig)
      }

     # Generate the gradient colors for continuous traits
      gradient_colors <- c("#E8F4FB", "#DCEEFF", "#CFE7FB", "#B4D8F0",
                           "#A5D0EA", "#91C7E7", "#77BBE2", "#6AB9DF", "#56ADD8",
                           "#4EA5D7", "#3696D2", "#2A8FCB", "#1481C4", "#2075B0",
                           "#1063A3", "#105CA0", "#083F80", "#08306B", "#052954")


    # Names based on whether it is a continuous trait or not
    name_extcon <- ifelse(pop.prev_ext == 0.5, paste0(name_ext, " -"), "controls")
    name_ext <- ifelse(pop.prev_ext == 0.5, paste0(name_ext, " +"), name_ext)
    h2_ext_name <- ifelse(pop.prev_ext == 0.5, "  (h2)", "  (h2.obs.50/50)")

    ### Put coordinates in vectors
    coord_con                       <- c(x.con, y.con, z.con)
    coord_sub1                      <- c(x.sub1, y.sub1,z.sub1)
    coord_sub2                      <- c(x.sub2,y.sub2,z.sub2)
    coord_allcases                  <- c(x.allcases, y.allcases, z.allcases)
    coord_popmean                   <- c(x.popmean, y.popmean, z.popmean)
    coord_ext                       <- c(x.ext,y.ext,z.ext)
    coord_extcontrol                <- c(x.extcontrol,y.extcontrol,z.extcontrol)

    # Calculate midpoints for arrow
    coord_arrow.sub1_line <- coord_con + (coord_sub1 - coord_con) / 2
    coord_arrow.sub2_line <- coord_con + (coord_sub2 - coord_con) / 2
    coord_arrow.allcases_line <- coord_con + (coord_allcases - coord_con) / 2
    coord_arrow.subsub_line <- coord_sub2 + (coord_sub1 - coord_sub2) / 2
    coord_arrow.ext_line <- coord_extcontrol + (coord_ext - coord_extcontrol) / 2
    coord_arrow.exttop_line <- coord_popmean + (coord_ext - coord_popmean) / 2
    coord_arrow.extbottom_line <- coord_popmean + (coord_extcontrol - coord_popmean) / 2

    # Calculate dir vectors of line
    dir.vect_sub1_line <- coord_sub1 - coord_con
    dir.vect_sub2_line <- coord_sub2 - coord_con
    dir.vect_subsub_line <- coord_sub1 - coord_sub2
    dir.vect_allcases_line <- coord_allcases - coord_con
    dir.vect_ext_line <- coord_ext - coord_extcontrol

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




#### Create 3D triangle --------------------------------------------------------------------------------



    # Initialize a 3D scatter plot
    fig <- plotly::plot_ly() %>%

      # Subsub line
      plotly::add_trace(x = c(x.sub1, x.sub2), y = c(y.sub1, y.sub2), z = c(z.sub1, z.sub2),
          type = 'scatter3d', mode = 'lines',
          line = list(color = "#7030A0", width = 6),
          showlegend = FALSE,
          hoverinfo = 'none') %>%

      # sub2_con line
      plotly::add_trace(x = c(x.sub2, x.con), y = c(y.sub2, y.con), z = c(z.sub2, z.con),
          type = 'scatter3d', mode = 'lines',
          line = list(color = "#00B050", width = 6),
          showlegend = FALSE,
          hoverinfo = 'none') %>%

      # sub1_con line
      plotly::add_trace(x = c(x.con, x.sub1), y = c(y.con, y.sub1), z = c(z.con, z.sub1),
              type = 'scatter3d', mode = 'lines',
              line = list(color = "#ED7D31", width = 6),
              showlegend = FALSE,
              hoverinfo = 'none') %>%

      # Allcases line
      plotly::add_trace(x = c(x.allcases, x.con), y = c(y.con, y.allcases), z = c(z.con, z.allcases),
          type = 'scatter3d', mode = 'lines',
          line = list(color = 'gray', width = 3),
          showlegend = FALSE,
          hoverinfo = 'none')

    # External trait line conditional on whether it is a continuous trait or not
    if (pop.prev_ext == 0.5) {

      fig <- fig %>%
        # External trait from popmean
        padd_explicit_color_line(
          x = seq( x.ext, x.popmean,length.out = length(gradient_colors)),
          y = seq( y.ext,y.popmean, length.out = length(gradient_colors)),
          z = seq( z.ext,z.popmean, length.out = length(gradient_colors)),
          colors = gradient_colors, width = 6 )  %>%

        # External trait to controls
        add_explicit_color_line(
          x = seq( x.extcontrol, x.popmean,length.out = length(gradient_colors)),
          y = seq(y.extcontrol,y.popmean,  length.out = length(gradient_colors)),
          z = seq( z.extcontrol, z.popmean,length.out = length(gradient_colors)),
          colors = gradient_colors, width = 6  )

       } else {

      fig <- fig %>%

        # External trait from popmean line
        plotly::add_trace(x = c(x.popmean,x.ext), y = c(y.popmean, y.ext), z = c(z.popmean, z.ext),
                  type = 'scatter3d', mode = 'lines',
                  line = list(color = "#052954", width = 6),
                  showlegend = FALSE,
                  hoverinfo = 'none') %>%

        # External trait to controls line
        plotly::add_trace(x = c(x.popmean,x.extcontrol), y = c(y.popmean, y.extcontrol), z = c(z.popmean, z.extcontrol),
                  type = 'scatter3d', mode = 'lines',
                  line = list(color = "#052954", width = 6),
                  showlegend = FALSE,
                  hoverinfo = 'none')        }

       fig <- fig %>%
       # Add popmean
         plotly::add_trace(x = x.popmean, y = y.popmean, z = z.popmean,
                          plotly = 'scatter3d', mode = 'markers',
          marker = list(size = 5, color = 'black', symbol = "cross"),
          showlegend = FALSE,
          hoverinfo = 'none') %>%

       # Add allcases point
         plotly::add_trace(x = x.allcases,
          y = y.allcases,
          z = z.allcases,
          type = 'scatter3d', mode = 'markers',
          marker = list(size = 4, color = 'black', symbol = "square"),
          showlegend = FALSE,
          hoverinfo = 'none') %>%

        # Add sub1 point
        plotly::add_trace(x = x.con,
                y = y.con,
                z = z.con,
                type = 'scatter3d', mode = 'markers',
                marker = list(size = 4, color = 'black', line = list(color = "black", width = 3), symbol = "circle"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%

        # Add sub1 point
        plotly::add_trace(x = x.sub1,
                y = y.sub1,
                z = z.sub1,
                type = 'scatter3d', mode = 'markers',
                marker = list(size = 4, color = 'black', line = list(color = "black", width = 3), symbol = "circle"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%

        # Add sub2 point
        plotly::add_trace(x = x.sub2,
                y = y.sub2,
                z = z.sub2,
                type = 'scatter3d', mode = 'markers',
                marker = list(size = 4, color = 'black', line = list(color = "black", width = 3), symbol = "circle"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%

       # Make the bottom triangle gray
       plotly::add_trace(x = c(x.sub1, x.sub2, x.con),
          y = c(y.sub1, y.sub2, y.con),
          z = c(z.sub1, z.sub2, z.con),
          type = 'mesh3d',
          i = c(0),  j = c(1), k = c(2),
          colorscale = list(c(0, "lightgray"), c(1, "lightgray")), # Define a gray colorscale
          intensity = rep(1, 3), # Use intensity to apply the colorscale uniformly
          opacity = 0.5,
          showlegend = FALSE, showscale = FALSE,
          hoverinfo = 'none')  %>%

        # Clean up plot
        plotly::layout(scene = list(
          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, showspikes = FALSE, title = ''),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, showspikes = FALSE, title = ''),
          zaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, showspikes = FALSE, title = ''),
          aspectmode = "data"),
          hovermode = 'none')

       # Add text annotations
       if (show.names == TRUE) {
       fig <- fig %>%
         plotly::add_trace(x = c(x.sub1, x.sub2, x.con, x.ext, x.extcontrol),
                   y = c(y.sub1, y.sub2, y.con, y.ext, y.extcontrol),
                   z = c(z.sub1, z.sub2, z.con, z.ext, z.extcontrol),
                   type = "scatter3d", mode = "text",
                   text = c(name_sub1, name_sub2, name_con, name_ext, name_extcon),
                   textposition = c("bottom right", "bottom right", "bottom right", "bottom right", "bottom right"),
                   textfont = list(size = 12),
                   showlegend = FALSE,
                   hoverinfo = 'none') }

      # Add arrowhead to the plot
      fig <- add_arrowhead(fig, dir.vect_sub1_line, coord_arrow.sub1_line,  color = "#ED7D31")
      fig <- add_arrowhead(fig, dir.vect_sub2_line, coord_arrow.sub2_line, color = "#00B050")
      fig <- add_arrowhead(fig, dir.vect_subsub_line, coord_arrow.subsub_line,color = "#7030A0")
      fig <- add_arrowhead(fig, dir.vect_allcases_line, coord_arrow.allcases_line,color = "gray")
      if (pop.prev_ext == 0.5) {
        fig <- add_arrowhead(fig, dir.vect_ext_line, coord_arrow.exttop_line,  color = "#77BBE2")
        fig <- add_arrowhead(fig, dir.vect_ext_line, coord_arrow.extbottom_line,color = "#77BBE2") } else {
        fig <- add_arrowhead(fig, dir.vect_ext_line, coord_arrow.ext_line, color = "#2075B0")    }


      # Save
      htmlwidgets::saveWidget(fig, paste0(folder_location, "/", filename3D,"_3D.plot.html"))

      # Return message
      cli::cli_alert_info(paste0(" Plot saved as ", folder_location,"/",filename3D,"_3D.plot.html"))

      # Box for the figure
      # fig <- fig %>%
      #   add_trace(x = c(-0.2, 0.2), y = c(0.35, 0.35), z = c(0.2, 0.2),
      #             type = 'scatter3d', mode = 'lines',
      #             showlegend = FALSE,
      #             line = list(color = "black", width = 5)) %>%
      #   add_trace(x = c(-0.2, 0.2), y = c(0.35, 0.35), z = c(-0.2, -0.2),
      #             type = 'scatter3d', mode = 'lines',
      #             line = list(color = "black", width = 5),
      #             showlegend = FALSE) %>%
      #   add_trace(x = c(-0.2, 0.2), y = c(-0.05,-0.05), z = c(0.2, 0.2),
      #             type = 'scatter3d', mode = 'lines',
      #             showlegend = FALSE,
      #             line = list(color = "black", width = 5)) %>%
      #   add_trace(x = c(-0.2, 0.2), y = c(-0.05,-0.05), z = c(-0.2, -0.2),
      #             type = 'scatter3d', mode = 'lines',
      #             showlegend = FALSE,
      #             line = list(color = "black", width = 5)) %>%
      #   add_trace(x = c(-0.2,-0.2), y = c(-0.05, 0.35), z = c(0.2,0.2),
      #             type = 'scatter3d', mode = 'lines',
      #             showlegend = FALSE,
      #             line = list(color = "black", width = 5)) %>%
      #   add_trace(x = c(-0.2,-0.2), y = c(-0.05, 0.35), z = c(-0.2,-0.2),
      #             type = 'scatter3d', mode = 'lines',
      #             showlegend = FALSE,
      #             line = list(color = "black", width = 5)) %>%
      #   add_trace(x = c(0.2,0.2), y = c(-0.05, 0.35), z = c(0.2,0.2),
      #             type = 'scatter3d', mode = 'lines',
      #             showlegend = FALSE,
      #             line = list(color = "black", width = 5)) %>%
      #   add_trace(x = c(0.2,0.2), y = c(-0.05, 0.35), z = c(-0.2,-0.2),
      #             type = 'scatter3d', mode = 'lines',
      #             showlegend = FALSE,
      #             line = list(color = "black", width = 5)) %>%
      #   add_trace(x = c(-0.2,-0.2), y = c(0.35, 0.35), z = c(-0.2, 0.2),
      #             type = 'scatter3d', mode = 'lines',
      #             showlegend = FALSE,
      #             line = list(color = "black", width = 5)) %>%
      #   add_trace(x = c(-0.2, -0.2), y = c(-0.05, -0.05), z = c(-0.2, 0.2),
      #             type = 'scatter3d', mode = 'lines',
      #             showlegend = FALSE,
      #             line = list(color = "black", width = 5)) %>%
      #   add_trace(x = c(0.2,0.2), y = c(0.35, 0.35), z = c(-0.2, 0.2),
      #             type = 'scatter3d', mode = 'lines',
      #             showlegend = FALSE,
      #             line = list(color = "black", width = 5)) %>%
      #   add_trace(x = c(0.2, 0.2), y = c(-0.05, -0.05), z = c(-0.2, 0.2),
      #             type = 'scatter3d', mode = 'lines',
      #             showlegend = FALSE,
      #             line = list(color = "black", width = 5))



#### Create legend --------------------------------------------------------------------------------



      ### Round values
      round_three_decimals <- function(objects) {
        for (i in objects) {assign(i, formatC(get(i,envir = temp_triangle_env), digits = 3, format = "f"),envir = temp_triangle_env)  }}
      round_no_decimals <- function(objects) {
        for (i in objects) {assign(i, formatC(get(i,envir = temp_triangle_env), digits = 0, format = "f"),envir = temp_triangle_env)  }}
      round_three_decimals(c("h2_ext","h2_allcases.con","h2_sub1.sub2","h2_sub1.con", "h2_sub2.con", "rg_sub1.con_ext", "rg_sub1.sub2_ext", "rg_sub2.con_ext", "rg_sub1.con_sub1.sub2", "rg_sub2.con_sub1.sub2", "rg_sub1.con_sub2.con", "rg_allcases.con_ext"))
      round_no_decimals(c("a.deg_sub1.con_ext","a.deg_sub1.sub2_ext","a.deg_sub2.con_ext", "a.deg_sub1.con_sub1.sub2", "a.deg_sub2.con_sub1.sub2","a.deg_sub1.con_sub2.con", "a.deg_allcases.con_ext"))

      # Make legend
      legend <- ggplot2::ggplot() +

        # Plot allcases and pop mean
        ggplot2::geom_point(aes(x = 0.1, y = -0.25), shape = 15, size = 4.5, color = "black") +
        ggplot2::geom_text(aes(x = 0.2, y = -0.25), label = name_allcases, color = "black", hjust = 0, size = 5) +
        ggplot2::geom_point(aes(x = 0.1, y = -0.125), shape = 3, size = 4, stroke = 2) +
        ggplot2::geom_text(aes(x = 0.2, y = -0.125), label = "Population mean", color = "black", hjust = 0, size = 5) +

        # Plot angles with degrees
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.125, 0.095)), color = "#00B050", linewidth = 2) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.125, 0.175)), color = "#ED7D31", linewidth = 2) +
        ggplot2::geom_text(aes(x = 0.3, y = 0.14), label = paste0(rg_sub1.con_sub2.con, " (", a.deg_sub1.con_sub2.con, "°)"), color = "black", hjust = 0, size = 6) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.26, 0.23)), color = "#7030A0", linewidth = 2) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.26, 0.31)), color = "#ED7D31", linewidth = 2) +
        ggplot2::geom_text(aes(x = 0.3, y = 0.27), label = paste0(rg_sub1.con_sub1.sub2, " (", a.deg_sub1.con_sub1.sub2, "°)"), color = "black", hjust = 0, size = 6) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.385, 0.355)), color = "#7030A0", linewidth = 2) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.385, 0.435)), color = "#00B050", linewidth = 2) +
        ggplot2::geom_text(aes(x = 0.3, y = 0.395), label = paste0(rg_sub2.con_sub1.sub2, " (", a.deg_sub2.con_sub1.sub2, "°)"), color = "black", hjust = 0, size = 6) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.51, 0.48)), color = "gray", linewidth = 2) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.51, 0.56)), color = "#2075B0", linewidth = 2) +
        ggplot2::geom_text(aes(x = 0.3, y = 0.52), label = paste0(rg_allcases.con_ext, " (", a.deg_allcases.con_ext, "°)"), color = "black", hjust = 0, size = 6) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.635, 0.605)), color = "#00B050", linewidth = 2) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.635, 0.685)), color = "#2075B0", linewidth = 2) +
        ggplot2::geom_text(aes(x = 0.3, y = 0.645), label = paste0(rg_sub2.con_ext, " (", a.deg_sub2.con_ext, "°)"), color = "black", hjust = 0, size = 6) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.765, 0.735)), color = "#ED7D31", linewidth = 2) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.765, 0.815)), color = "#2075B0", linewidth = 2) +
        ggplot2::geom_text(aes(x = 0.3, y = 0.775), label = paste0(rg_sub1.con_ext, " (", a.deg_sub1.con_ext, "°)"), color = "black", hjust = 0, size = 6) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.89, 0.86)), color = "#7030A0", linewidth = 2) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(0.89, 0.94)), color = "#2075B0", linewidth = 2) +
        ggplot2::geom_text(aes(x = 0.3, y = 0.9), label = paste0(rg_sub1.sub2_ext, " (", a.deg_sub1.sub2_ext, "°)"), color = "black", hjust = 0, size = 6) +
        ggplot2::annotate("text", x = 0, y = 1.05, label = "italic(r)[g] ~ '(degrees)'", color = "black", hjust = 0, size = 6.5, parse = TRUE) +

        # plot h2 lines
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(1.25, 1.25)), color = "#2075B0", linewidth = 2) +
        ggplot2::geom_text(aes(x = 0.3, y = 1.26), label = h2_ext, color = "black", hjust = 0, size = 6) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(1.375, 1.375)), color = "gray", linewidth = 2) +
        ggplot2::geom_text(aes(x = 0.3, y = 1.385), label = h2_allcases.con, color = "black", hjust = 0, size = 6) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(1.5, 1.5)), color = "#7030A0", linewidth = 2) +
        ggplot2::geom_text(aes(x = 0.3, y = 1.51), label = h2_sub1.sub2, color = "black", hjust = 0, size = 6) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(1.625, 1.625)), color = "#ED7D31", linewidth = 2) +
        ggplot2::geom_text(aes(x = 0.3, y = 1.635), label = h2_sub1.con, color = "black", hjust = 0, size = 6) +
        ggplot2::geom_line(aes(x = c(0.1, 0.25), y = c(1.75, 1.75)), color = "#00B050", linewidth = 2) +
        ggplot2::geom_text(aes(x = 0.3, y = 1.76), label = h2_sub2.con, color = "black", hjust = 0, size = 6) +
        #annotate("text", x = 0, y = 1.91, label = "italic(h)^2 ~ italic(obs)[50/50]", color = "black", hjust = 0, size = 6.5, parse = TRUE) +
        ggplot2::annotate("text", x = 0, y = 1.91, label = "Heritability", color = "black", hjust = 0, size = 6.5) +

        # Clean background
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.title.x = element_blank(), axis.title.y = element_blank(),   # Remove axes titles
              axis.text.x = element_blank(),  axis.text.y = element_blank(),    # Remove axes text
              axis.ticks = element_blank(),                                     # Remove axes ticks
              panel.grid.major = element_blank(),                               # Remove major gridlines
              panel.grid.minor = element_blank(),                               # Remove minor gridlines
              panel.background = element_rect(fill = "white", colour = NA)) +   # Set background color to white
        ggplot2::xlim(0, 1.6) +
        # Scale axes
        ggplot2::coord_fixed(ratio = 1)


      ### Save legend
      suppressWarnings({ ggplot2::ggsave(paste0(folder_location, "/", filename3D,"_3D.legend.png"), plot = legend, width = 4, height = 8, dpi = 300) })


      # Return message
      cli::cli_alert_info(paste0(" Legend saved as ", folder_location,"/",filename3D,"_3D.legend.png"))
      cli::cli_alert_success("GDVIS plot succesfully finished")
      if (show.rendering == TRUE) {  cli::cli_alert_info("Showing plot in external window ") }



#### Make shiny object to have plot and legend together --------------------------------------------------------------------------------



      if (webversion == FALSE) {

      # Define UI for the shiny app
      ui <- shiny::fluidPage(shiny::fluidRow(
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
      save_path <- paste0(folder_location, "/",plot_title, ".", name_ext, "_3D.plot_files")
      if (dir.exists(save_path)) { unlink(save_path, recursive = T)   }

      # Save the Shiny app script as a .R file
      saveRDS(shiny_app, file = paste0(folder_location, "/", filename3D, "_shiny_app.rds"))

      # Save the UI and server components separately
      saveRDS(ui, file = paste0(folder_location, "/", filename3D, "_ui.rds"))
      saveRDS(server, file = paste0(folder_location, "/", filename3D, "_server.rds"))

      # Show if asked
      if (show.rendering == TRUE) { assign("shiny_app", shiny_app, envir = .GlobalEnv); runApp(get("shiny_app", envir = .GlobalEnv)) }
      }

      if (webversion == TRUE) {
        # Return plot objects
        return(list(
          fig = fig,
          legend = legend  ))
      }


    }) # End with env

}  # End function
