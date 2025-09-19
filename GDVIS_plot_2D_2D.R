library(tidyverse)
library(plotly)
library(data.table)
library(htmlwidgets)
library(shiny)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(cli)


# Function to plot 2D_2D triangles 
GDVIS_plot_2D_2D <- function(input_triangle_parameters, show.rendering = TRUE, webversion = FALSE) {
  
  
  
#### Set up the environment and load data -----------------------------------------------------------
  
  
    # Create a new environment 
    temp_triangle_env <- new.env()
  
    # Load the parameters file into this new environment
    load(input_triangle_parameters, envir = temp_triangle_env)
    
    # Extract the list from the environment
    triangle_list <- temp_triangle_env$double.triangle.output.list
    
    # Show message  
    cli::cli_h1(paste0("Running ", cli::col_cyan("GDVIS plot 2D_2D"), " on ", cli::col_cyan(triangle_list$triangle1.filename), " with " , cli::col_cyan(triangle_list$triangle2.filename))) 
  
    # Directly access the objects
    with(temp_triangle_env, {
    
    # Extract the variables of the list
    list2env(double.triangle.output.list, envir = temp_triangle_env)
    
    
    
    
#### Calculate extra points --------------------------------------------------------------------------------     
    
    
    
    # Put coordinates in vectors
    triangle1.coord_sub1 <- c(triangle1.x.sub1, triangle1.y.sub1, triangle1.z.sub1)
    triangle1.coord_sub2 <- c(triangle1.x.sub2, triangle1.y.sub2, triangle1.z.sub2)
    triangle2.coord_sub1 <- c(triangle2.x.sub1, triangle2.y.sub1, triangle2.z.sub1)
    triangle2.coord_sub2 <- c(triangle2.x.sub2, triangle2.y.sub2, triangle2.z.sub2)
    coord_con <- c(triangle1.x.con, triangle1.y.con, triangle1.z.con)
    coord_allcases <- c(triangle1.x.allcases, triangle1.y.allcases, triangle1.z.allcases)
    
    # Calculate midpoints for arrow 
    triangle1.coord_arrow.sub1_line <- coord_con + (triangle1.coord_sub1 - coord_con) / 2
    triangle1.coord_arrow.sub2_line <- coord_con + (triangle1.coord_sub2 - coord_con) / 2
    triangle1.coord_arrow.subsub_line <- triangle1.coord_sub2 + (triangle1.coord_sub1 - triangle1.coord_sub2) / 2
    triangle1.coord_arrow.allcases_line <- coord_con + (coord_allcases - coord_con) / 2
    triangle2.coord_arrow.sub1_line <- coord_con + (triangle2.coord_sub1 - coord_con) / 2
    triangle2.coord_arrow.sub2_line <- coord_con + (triangle2.coord_sub2 - coord_con) / 2
    triangle2.coord_arrow.subsub_line <- triangle2.coord_sub2 + (triangle2.coord_sub1 - triangle2.coord_sub2) / 2
    
    # Calculate dir vectors of line
    triangle1.dir.vect_sub1_line <- triangle1.coord_sub1 - coord_con
    triangle1.dir.vect_sub2_line <- triangle1.coord_sub2 - coord_con
    triangle1.dir.vect_subsub_line <- triangle1.coord_sub1 - triangle1.coord_sub2
    triangle1.dir.vect_allcases_line <- coord_allcases - coord_con
    triangle2.dir.vect_sub1_line <- triangle2.coord_sub1 - coord_con
    triangle2.dir.vect_sub2_line <- triangle2.coord_sub2 - coord_con
    triangle2.dir.vect_subsub_line <- triangle2.coord_sub1 - triangle2.coord_sub2
    
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
        add_trace(x = c(midpoint[1], arrow_start_right[1]), 
                  y = c(midpoint[2], arrow_start_right[2]), 
                  z = c(midpoint[3], arrow_start_right[3]),
                  type = 'scatter3d', mode = 'lines', line = list(color = color, width = 3), showlegend = FALSE) %>%
        add_trace(x = c(midpoint[1], arrow_start_left[1]), 
                  y = c(midpoint[2], arrow_start_left[2]), 
                  z = c(midpoint[3], arrow_start_left[3]),
                  type = 'scatter3d', mode = 'lines', line = list(color = color, width = 3), showlegend = FALSE) %>%
        add_trace(x = c(midpoint[1], arrow_start_up[1]), 
                  y = c(midpoint[2], arrow_start_up[2]), 
                  z = c(midpoint[3], arrow_start_up[3]),
                  type = 'scatter3d', mode = 'lines', line = list(color = color, width = 3), showlegend = FALSE) %>%
        add_trace(x = c(midpoint[1], arrow_start_down[1]), 
                  y = c(midpoint[2], arrow_start_down[2]), 
                  z = c(midpoint[3], arrow_start_down[3]),
                  type = 'scatter3d', mode = 'lines', line = list(color = color, width = 3), showlegend = FALSE) %>%
        add_trace(x = c(arrow_start_left[1], arrow_start_right[1]), 
                  y = c(arrow_start_left[2], arrow_start_right[2]), 
                  z = c(arrow_start_left[3], arrow_start_right[3]),
                  type = 'scatter3d', mode = 'lines', line = list(color = color, width = 3), showlegend = FALSE) %>%
        add_trace(x = c(arrow_start_up[1], arrow_start_down[1]), 
                  y = c(arrow_start_up[2], arrow_start_down[2]), 
                  z = c(arrow_start_up[3], arrow_start_down[3]),
                  type = 'scatter3d', mode = 'lines', line = list(color = color, width = 3), showlegend = FALSE) %>%
        add_trace(x = c(arrow_start_right[1], arrow_start_left[1], midpoint[1]),
                  y = c(arrow_start_right[2], arrow_start_left[2], midpoint[2]),
                  z = c(arrow_start_right[3], arrow_start_left[3], midpoint[3]),
                  type = 'mesh3d',
                  i = c(0), j = c(1), k = c(2),  
                  vertexcolor = matrix(rep(col2rgb(color)/255, 3), nrow = 3, byrow = TRUE),  # Apply the exact color to each vertex
                  opacity = 1,
                  showlegend = FALSE, showscale = FALSE,
                  hoverinfo = 'none') %>%
        add_trace(x = c(arrow_start_up[1], arrow_start_down[1], midpoint[1]),
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
    fig <- plot_ly() %>%
      
      # Subsub line 
      add_trace(x = c(triangle1.x.sub1, triangle1.x.sub2), y = c(triangle1.y.sub1, triangle1.y.sub2), z = c(triangle1.z.sub1, triangle1.z.sub2),
                type = 'scatter3d', mode = 'lines',
                line = list(color = "#7030A0", width = 9),
                name = paste0(triangle1.filename, ": ", triangle1.h2_sub1.sub2,"  (h2.obs.50/50)"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # Subsub line 
      add_trace(x = c(triangle2.x.sub1, triangle2.x.sub2), y = c(triangle2.y.sub1, triangle2.y.sub2), z = c(triangle2.z.sub1, triangle2.z.sub2),
                type = 'scatter3d', mode = 'lines',
                line = list(color = "#7030A0", width = 9, dash = "dashdot"),
                name = paste0(triangle2.filename, ": ", triangle2.h2_sub1.sub2,"  (h2.obs.50/50)"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # sub2_con line
      add_trace(x = c(triangle1.x.sub2, triangle1.x.con), y = c(triangle1.y.sub2, triangle1.y.con), z = c(triangle1.z.sub2, triangle1.z.con),
                type = 'scatter3d', mode = 'lines',
                line = list(color = "#00B050", width = 9),
                name = paste0(triangle1.filename, ": ", triangle1.h2_sub2.con,"  (h2.obs.50/50)"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # sub2_con line
      add_trace(x = c(triangle2.x.sub2, triangle2.x.con), y = c(triangle2.y.sub2, triangle2.y.con), z = c(triangle2.z.sub2, triangle2.z.con),
                type = 'scatter3d', mode = 'lines',
                line = list(color ="#00B050", width = 9, dash = "dashdot"),
                name = paste0(triangle2.filename, ": ", triangle2.h2_sub2.con,"  (h2.obs.50/50)"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # sub1_con line
      add_trace(x = c(triangle1.x.con, triangle1.x.sub1), y = c(triangle1.y.con, triangle1.y.sub1), z = c(triangle1.z.con, triangle1.z.sub1),
                type = 'scatter3d', mode = 'lines',
                line = list(color = "#ED7D31", width = 9),
                name = paste0(triangle1.filename, ": ", triangle1.h2_sub1.con,"  (h2.obs.50/50)"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # sub1_con line
      add_trace(x = c(triangle2.x.con, triangle2.x.sub1), y = c(triangle2.y.con, triangle2.y.sub1), z = c(triangle2.z.con, triangle2.z.sub1),
                type = 'scatter3d', mode = 'lines',
                line = list(color = "#ED7D31", width = 9, dash = "dashdot"),
                name = paste0(triangle2.filename, ": ", triangle2.h2_sub1.con,"  (h2.obs.50/50)"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # Allcases line
      add_trace(x = c(triangle1.x.allcases, triangle1.x.con), y = c(triangle1.y.con, triangle1.y.allcases), z = c(triangle1.z.con, triangle1.z.allcases),
                type = 'scatter3d', mode = 'lines',
                line = list(color = 'gray', width = 6),
                name = paste0(triangle2.h2_allcases.con,"  (h2.obs.50/50)"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # Add popmean point
      add_trace(x = triangle1.x.popmean, y = triangle1.y.popmean, z = triangle1.z.popmean,
                type = 'scatter3d', mode = 'markers',
                marker = list(size = 6, color = 'black', symbol = "cross"),
                showlegend = FALSE,
                name = "pop mean",
                hoverinfo = 'none') %>%
      
      # Add allcases point
      add_trace(x = triangle1.x.allcases, 
                y = triangle1.y.allcases, 
                z = triangle1.z.allcases,
                type = 'scatter3d', mode = 'markers',
                marker = list(size = 5, color = 'black', symbol = "square"),
                showlegend = FALSE,
                name = triangle1.name_allcases,
                hoverinfo = 'none') %>%
      
      # Add con point
      add_trace(x = triangle1.x.con, 
                y = triangle1.y.con, 
                z = triangle1.z.con,
                type = 'scatter3d', mode = 'markers',
                marker = list(size = 6, color = 'black', line = list(color = "black", width = 3), symbol = "circle"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # Add sub1 point
      add_trace(x = triangle1.x.sub1, 
                y = triangle1.y.sub1, 
                z = triangle1.z.sub1,
                type = 'scatter3d', mode = 'markers',
                marker = list(size = 6, color = 'white', line = list(color = "black", width = 3), symbol = "circle"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # Add sub1 point
      add_trace(x = triangle2.x.sub1, 
                y = triangle2.y.sub1, 
                z = triangle2.z.sub1,
                type = 'scatter3d', mode = 'markers',
                marker = list(size = 6, color = 'white', line = list(color = "black", width = 3), symbol = "circle"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # Add sub2 point
      add_trace(x = triangle2.x.sub2, 
                y = triangle2.y.sub2, 
                z = triangle2.z.sub2,
                type = 'scatter3d', mode = 'markers',
                marker = list(size = 6, color = 'darkgray', line = list(color = "black", width = 3), symbol = "circle"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # Add sub2 point
      add_trace(x = triangle1.x.sub2, 
                y = triangle1.y.sub2, 
                z = triangle1.z.sub2,
                type = 'scatter3d', mode = 'markers',
                marker = list(size = 6, color = 'darkgray', line = list(color = "black", width = 3), symbol = "circle"),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # Add text annotations 
      add_trace(x = c(triangle1.x.sub1, triangle1.x.sub2, triangle1.x.con, triangle2.x.sub1, triangle2.x.sub2, triangle2.x.con),
                y = c(triangle1.y.sub1, triangle1.y.sub2, triangle1.y.con, triangle2.y.sub1, triangle2.y.sub2, triangle2.y.con),
                z = c(triangle1.z.sub1, triangle1.z.sub2, triangle1.z.con, triangle2.z.sub1, triangle2.z.sub2, triangle2.z.con),
                type = "scatter3d", mode = "text",
                text = c(triangle1.name_sub1, triangle1.name_sub2, triangle1.name_con, triangle2.name_sub1, triangle2.name_sub2, triangle2.name_con),
                textposition = c("bottom right", "bottom right", "bottom right","bottom right", "bottom right", "bottom right"),
                textfont = list(size = 12),
                showlegend = FALSE,
                hoverinfo = 'none') %>%
      
      # Make the bottom triangle gray  
      add_trace(x = c(triangle1.x.sub1, triangle1.x.sub2, triangle1.x.con),
                y = c(triangle1.y.sub1, triangle1.y.sub2, triangle1.y.con),
                z = c(triangle1.z.sub1, triangle1.z.sub2, triangle1.z.con),
                type = 'mesh3d',
                i = c(0),  # Indices of x starting from 0
                j = c(1),  # Indices of y
                k = c(2),  # Indices of z
                colorscale = list(c(0, "lightgray"), c(1, "lightgray")), # Define a gray colorscale
                intensity = rep(1, 3), # Use intensity to apply the colorscale uniformly
                opacity = 0.5,
                showlegend = FALSE, showscale = FALSE,
                hoverinfo = 'none')  %>%
      
      # Make the bottom triangle gray  
      add_trace(x = c(triangle2.x.sub1, triangle2.x.sub2, triangle2.x.con),
                y = c(triangle2.y.sub1, triangle2.y.sub2, triangle2.y.con),
                z = c(triangle2.z.sub1, triangle2.z.sub2, triangle2.z.con),
                type = 'mesh3d',
                i = c(0),  # Indices of x starting from 0
                j = c(1),  # Indices of y
                k = c(2),  # Indices of z
                colorscale = list(c(0, "lightgray"), c(1, "lightgray")), # Define a gray colorscale
                intensity = rep(1, 3), # Use intensity to apply the colorscale uniformly
                opacity = 0.5,
                showlegend = FALSE, showscale = FALSE,
                hoverinfo = 'none')  %>%
      
      # Clean up plot
      layout(scene = list(
        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, showspikes = FALSE, title = ''), 
        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, showspikes = FALSE, title = ''),
        zaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, showspikes = FALSE, title = ''),
        aspectmode = "data"),
        hovermode = 'none') 
    
    # Add arrowhead to the plot
    fig <- add_arrowhead(fig, triangle1.dir.vect_sub1_line, triangle1.coord_arrow.sub1_line,  color = "#ED7D31")
    fig <- add_arrowhead(fig, triangle1.dir.vect_sub2_line, triangle1.coord_arrow.sub2_line, color = "#00B050")
    fig <- add_arrowhead(fig, triangle1.dir.vect_subsub_line, triangle1.coord_arrow.subsub_line,color = "#7030A0")
    fig <- add_arrowhead(fig, triangle1.dir.vect_allcases_line, triangle1.coord_arrow.allcases_line,color = "gray")
    fig <- add_arrowhead(fig, triangle2.dir.vect_sub1_line, triangle2.coord_arrow.sub1_line,  color = "#ED7D31")
    fig <- add_arrowhead(fig, triangle2.dir.vect_sub2_line, triangle2.coord_arrow.sub2_line, color = "#00B050")
    fig <- add_arrowhead(fig, triangle2.dir.vect_subsub_line, triangle2.coord_arrow.subsub_line,color = "#7030A0")
    
    # # Box for the figure
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
    
    
    # Save
    htmlwidgets::saveWidget(fig, paste0(triangle1.folder_location, "/", triangle1.filename, ".vs.",triangle2.filename,".html"))
    
    # Return message
    cli::cli_alert_info(paste0(" Plot saved as ", triangle1.folder_location, "/", triangle1.filename, ".vs.",triangle2.filename,".html"))
    
    
    
#### Create legend --------------------------------------------------------------------------------            
    
    
    
    ### Round values
    round_three_decimals <- function(objects) {
      for (i in objects) {assign(i, formatC(get(i,envir = temp_triangle_env), digits = 3, format = "f"),envir = temp_triangle_env)  }}
    round_no_decimals <- function(objects) {
      for (i in objects) {assign(i, formatC(get(i,envir = temp_triangle_env), digits = 0, format = "f"),envir = temp_triangle_env)  }}

    round_three_decimals(c("triangle1.h2_sub2.con", "triangle1.h2_sub1.con", "triangle1.h2_sub1.sub2", "triangle2.h2_sub2.con", "triangle2.h2_sub1.con", "triangle2.h2_sub1.sub2", "triangle2.h2_allcases.con",
                           "triangle1.rg_sub1.con_sub1.sub2", "triangle1.rg_sub2.con_sub1.sub2", "triangle1.rg_sub1.con_sub2.con",
                           "triangle2.rg_sub2.con_sub1.sub2", "triangle2.rg_sub1.con_sub2.con","triangle2.rg_sub1.con_sub1.sub2",
                           "rg_rsub1.rsub2_sub1.sub2","rg_rsub1.con_sub1.con","rg_rsub2.con_sub2.con","rg_rsub1.con_sub2.con","rg_rsub2.con_sub1.con"))
    round_no_decimals(c("triangle2.a.deg_sub1.con_sub1.sub2","triangle2.a.deg_sub1.con_sub2.con","triangle2.a.deg_sub2.con_sub1.sub2","triangle1.a.deg_sub1.con_sub2.con","triangle1.a.deg_sub2.con_sub1.sub2","triangle1.a.deg_sub1.con_sub1.sub2",
                        "a.deg_rsub2.con_sub1.con", "a.deg_rsub1.con_sub2.con", "a.deg_rsub2.con_sub2.con","a.deg_rsub1.con_sub1.con","a.deg_rsub1.rsub2_sub1.sub2"))
  
    # Make legend
    legend <- ggplot() +
      
      # Plot allcases and pop mean
      geom_point(aes(x = 0, y = -1.435), shape = 15, size = 4.5, color = "black") +
      geom_text(aes(x = 0.115, y = -1.42), label = triangle1.name_allcases, color = "black", hjust = 0, size = 5) +
      geom_point(aes(x = 0, y = -1.31), shape = 4, size = 4.5, stroke = 2) +
      geom_text(aes(x = 0.115, y = -1.31), label = "Population mean", color = "black", hjust = 0, size = 5) +
      
      # Plot angles
      # Triangle relations
      geom_line(aes(x = c(0.1, 0.23), y = c(-1.08, -1.11)), color = "#ED7D31", linewidth = 1.5) +
      geom_line(aes(x = c(0.1, 0.23), y = c(-1.08, -1.03)), color = "#00B050", linewidth = 1.5, linetype = "11") +
      geom_text(aes(x = 0.3, y = -1.07), label = paste0(rg_rsub2.con_sub1.con, " (", a.deg_rsub2.con_sub1.con, "°)"), color = "black", hjust = 0, size = 6) +
      geom_line(aes(x = c(0.1, 0.23), y = c(-0.955, -0.985)), color = "#00B050", linewidth = 1.5) +
      geom_line(aes(x = c(0.1, 0.23), y = c(-0.955, -0.91)), color = "#ED7D31", linewidth = 1.5, linetype = "11") +
      geom_text(aes(x = 0.3, y = -0.945), label = paste0(rg_rsub1.con_sub2.con, " (", a.deg_rsub1.con_sub2.con, "°)"), color = "black", hjust = 0, size = 6) +
      geom_line(aes(x = c(0.1, 0.23), y = c(-0.83, -0.86)), color = "#00B050", linewidth = 1.5) +
      geom_line(aes(x = c(0.1, 0.23), y = c(-0.83, -0.78)), color = "#00B050", linewidth = 1.5, linetype = "11") +
      geom_text(aes(x = 0.3, y = -0.815), label = paste0(rg_rsub2.con_sub2.con, " (", a.deg_rsub2.con_sub2.con, "°)"), color = "black", hjust = 0, size = 6) +
      geom_line(aes(x = c(0.1, 0.23), y = c(-0.705, -0.66)), color = "#ED7D31", linewidth = 1.5) +
      geom_line(aes(x = c(0.1, 0.23), y = c(-0.705, -0.74)), color = "#ED7D31", linewidth = 1.5, linetype = "11") +
      geom_text(aes(x = 0.3, y = -0.695), label = paste0(rg_rsub1.con_sub1.con, " (", a.deg_rsub1.con_sub1.con, "°)"), color = "black", hjust = 0, size = 6) +
      geom_line(aes(x = c(0.1, 0.23), y = c(-0.58, -0.61)), color = "#7030A0", linewidth = 1.5) +
      geom_line(aes(x = c(0.1, 0.23), y = c(-0.58, -0.53)), color = "#7030A0", linewidth = 1.5, linetype = "11") +
      geom_text(aes(x = 0.3, y = -0.57), label = paste0(rg_rsub1.rsub2_sub1.sub2, " (", a.deg_rsub1.rsub2_sub1.sub2, "°)"), color = "black", hjust = 0, size = 6) +
      annotate("text", x = 0, y = -0.445, label = "Combined triangles", color = "black", hjust = 0, size = 5) +
      
      # Triangle 2
      geom_line(aes(x = c(0.1, 0.23), y = c(-0.25, -0.28)), color = "#7030A0", linewidth = 1.5, linetype = "11") +
      geom_line(aes(x = c(0.1, 0.23), y = c(-0.25, -0.2)), color = "#ED7D31", linewidth = 1.5, linetype = "11") +
      geom_text(aes(x = 0.3, y = -0.24), label = paste0(triangle2.rg_sub1.con_sub1.sub2, " (", triangle2.a.deg_sub1.con_sub1.sub2, "°)"), color = "black", hjust = 0, size = 6) +
      geom_line(aes(x = c(0.1, 0.23), y = c(-0.125, -0.155)), color = "#ED7D31", linewidth = 1.5, linetype = "11") +
      geom_line(aes(x = c(0.1, 0.23), y = c(-0.125, -0.07)), color = "#00B050", linewidth = 1.5, linetype = "11") +
      geom_text(aes(x = 0.3, y = -0.11), label = paste0(triangle2.rg_sub1.con_sub2.con, " (", triangle2.a.deg_sub1.con_sub2.con, "°)"), color = "black", hjust = 0, size = 6) +
      geom_line(aes(x = c(0.1, 0.23), y = c(0.005, -0.025)), color = "#7030A0", linewidth = 1.5, linetype = "11") +
      geom_line(aes(x = c(0.1, 0.23), y = c(0.005, 0.055)), color = "#00B050", linewidth = 1.5, linetype = "11") +
      geom_text(aes(x = 0.3, y = 0.02), label = paste0(triangle2.rg_sub2.con_sub1.sub2, " (", triangle2.a.deg_sub2.con_sub1.sub2, "°)"), color = "black", hjust = 0, size = 6) +
      annotate("text", x = 0, y = 0.14, label = triangle2.filename, color = "black", hjust = 0, size = 5) +
      
      # Triangle 1
      geom_line(aes(x = c(0.1, 0.23), y = c(0.34, 0.31)), color = "#ED7D31", linewidth = 1.5) +
      geom_line(aes(x = c(0.1, 0.23), y = c(0.34, 0.39)), color = "#00B050", linewidth = 1.5) +
      geom_text(aes(x = 0.3, y = 0.35), label = paste0(triangle1.rg_sub1.con_sub2.con, " (", triangle1.a.deg_sub1.con_sub2.con, "°)"), color = "black", hjust = 0, size = 6) +
      geom_line(aes(x = c(0.1, 0.23), y = c(0.47, 0.44)), color = "#7030A0", linewidth = 1.5) +
      geom_line(aes(x = c(0.1, 0.23), y = c(0.47, 0.52)), color = "#00B050", linewidth = 1.5) +
      geom_text(aes(x = 0.3, y = 0.475), label = paste0(triangle1.rg_sub2.con_sub1.sub2, " (", triangle1.a.deg_sub2.con_sub1.sub2, "°)"), color = "black", hjust = 0, size = 6) +
      geom_line(aes(x = c(0.1, 0.23), y = c(0.595, 0.565)), color = "#7030A0", linewidth = 1.5) +
      geom_line(aes(x = c(0.1, 0.23), y = c(0.595, 0.645)), color = "#ED7D31", linewidth = 1.5) +
      geom_text(aes(x = 0.3, y = 0.595), label = paste0(triangle1.rg_sub1.con_sub1.sub2, " (", triangle1.a.deg_sub1.con_sub1.sub2, "°)"), color = "black", hjust = 0, size = 6) +
      annotate("text", x = 0, y = 0.725, label = triangle1.filename, color = "black", hjust = 0, size = 5) +
      annotate("text", x = -0.1, y = 0.875, label = "italic(r)[g] ~ '(degrees)'", color = "black", hjust = 0, size = 6.5, parse = TRUE) +

      
      # plot h2 lines
      # Allcases
      geom_line(aes(x = c(0.1, 0.25), y = c(1.13, 1.13)), color = "gray", linewidth = 2) +
      geom_text(aes(x = 0.3, y = 1.14), label = triangle2.h2_allcases.con, color = "black", hjust = 0, size = 6) +
      annotate("text", x = 0, y = 1.26, label = triangle1.name_allcases, color = "black", hjust = 0, size = 5) +
      
      # Triangle 2
      geom_line(aes(x = c(0.1, 0.25), y = c(1.4, 1.4)), color = "#7030A0", linewidth = 2, linetype = "11") +
      geom_text(aes(x = 0.3, y = 1.41), label = triangle2.h2_sub1.sub2, color = "black", hjust = 0, size = 6) +
      geom_line(aes(x = c(0.1, 0.25), y = c(1.525, 1.525)), color = "#ED7D31", linewidth = 2, linetype = "11") +
      geom_text(aes(x = 0.3, y = 1.535), label = triangle2.h2_sub1.con, color = "black", hjust = 0, size = 6) +
      geom_line(aes(x = c(0.1, 0.25), y = c(1.65, 1.65)), color = "#00B050", linewidth = 2, linetype = "11") +
      geom_text(aes(x = 0.3, y = 1.66), label = triangle2.h2_sub2.con, color = "black", hjust = 0, size = 6) +
      annotate("text", x = 0, y = 1.78, label = triangle2.filename, color = "black", hjust = 0, size = 5) +
      
      # Triangle 1
      geom_line(aes(x = c(0.1, 0.25), y = c(1.92)), color = "#7030A0", linewidth = 2) +
      geom_text(aes(x = 0.3, y = 1.93), label = triangle1.h2_sub1.sub2, color = "black", hjust = 0, size = 6) +
      geom_line(aes(x = c(0.1, 0.25), y = c(2.045, 2.045)), color = "#ED7D31", linewidth = 2) +
      geom_text(aes(x = 0.3, y = 2.055), label = triangle1.h2_sub1.con, color = "black", hjust = 0, size = 6) +
      geom_line(aes(x = c(0.1, 0.25), y = c(2.17, 2.17)), color = "#00B050", linewidth = 2) +
      geom_text(aes(x = 0.3, y = 2.18), label = triangle1.h2_sub2.con, color = "black", hjust = 0, size = 6) +
      annotate("text", x = 0, y = 2.3, label = triangle1.filename, color = "black", hjust = 0, size = 5) +
      annotate("text", x = -0.1, y = 2.48, label = "Heritability", color = "black", hjust = 0, size = 6.5) +
      
      # Clean background
      theme_minimal() +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(),   # Remove axes titles
            axis.text.x = element_blank(),  axis.text.y = element_blank(),    # Remove axes text
            axis.ticks = element_blank(),                                     # Remove axes ticks
            panel.grid.major = element_blank(),                               # Remove major gridlines
            panel.grid.minor = element_blank(),                               # Remove minor gridlines
            panel.background = element_rect(fill = "white", colour = NA)) +   # Set background color to white
      xlim(-0.1, 1.6) +
      # Scale axes
      coord_fixed(ratio = 1)
    
    
    ### Save legend
    ggsave(paste0(triangle1.folder_location, "/", triangle1.filename, ".vs.", triangle2.filename, "2D_2Dlegend.png"), plot = legend, width = 4, height = 8, dpi = 300) 
    
    # Return message
    cli::cli_alert_info(paste0(" Legend saved as ", triangle1.folder_location, "/", triangle1.filename, ".vs.", triangle2.filename, "2D_2Dlegend.png"))
    cli::cli_alert_success("GDVIS plot succesfully finished")
    if (show.rendering == TRUE) {  cli::cli_alert_info("Showing plot in external window ") }
    
    
    
  #### Make shiny object to have plot and legend together --------------------------------------------------------------------------------     
    
    
    if (webversion == FALSE) {
    
    # Define UI for the shiny app
    ui <- fluidPage(fluidRow(
      column(8, plotlyOutput("plotly_plot", height = "600px")),   # 8/12 width for the plotly plot
      column(4, plotOutput("ggplot_legend", height = "600px"))))    # 4/12 width for the ggplot legend
    
    # Define server logic for the shiny app
    server <- function(input, output, session) {
      
      # Render the plotly plot (fig)
      output$plotly_plot <- renderPlotly({  fig  })
      
      # Render the ggplot legend
      output$ggplot_legend <-  suppressWarnings({ renderPlot({  legend  })   })  }
    
    # Run the shiny application
    shiny_app <- shinyApp(ui = ui, server = server)
    
    # Remove the folder if it already already exists (shiny doesn't overwrite)
    save_path <- paste0(triangle1.folder_location, "/", triangle1.filename, ".vs.", triangle2.filename, "_files")
    if (dir.exists(save_path)) { 
      Sys.chmod(save_path, mode = "777", use_umask = FALSE)  # Grant full permissions
      unlink(save_path, recursive = TRUE, force = TRUE)      # Force delete
      }
    
    # Save the Shiny app script as a .R file
    saveRDS(shiny_app, file = paste0(triangle1.folder_location, "/", triangle1.filename, ".vs.", triangle2.filename, ".double_triangles_shiny_app.rds"))
    
    # Save the UI and server components separately
    saveRDS(ui, file = paste0(triangle1.folder_location, "/", triangle1.filename, ".vs.", triangle2.filename, ".double_triangles_ui.rds"))
    saveRDS(server, file = paste0(triangle1.folder_location, "/", triangle1.filename, ".vs.", triangle2.filename, ".double_triangles_server.rds"))
    
    # Show if asked, otherwise print that the plots are saved
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