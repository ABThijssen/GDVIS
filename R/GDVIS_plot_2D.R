#' Plot GDVIS 2D plot
#'
#' This function plots a GDVIS visualization based on GDVIS parameters.

#' @param input_triangle_parameters the path to the RData GDVIS output, which ends with 2D.triangle_parameters.RData
#' @param x_lower set lower value for the x axis
#' @param x_upper set upper value for the x axis
#' @param y_lower set lower value for the y asis
#' @param y_upper set upper value for the y axis
#' @export
# Function to plot 2D triangles
GDVIS_plot_2D <- function(input_triangle_parameters, x_lower = NULL, x_upper = NULL, y_lower = NULL, y_upper = NULL ) {

#### Set up the environment and load data -----------------------------------------------------------


    # Create a new environment
    temp_triangle_env <- new.env()

    # Load the parameters file into this new environment
    load(input_triangle_parameters, envir = temp_triangle_env)
    temp_triangle_env$x_lower <- x_lower
    temp_triangle_env$x_upper <- x_upper
    temp_triangle_env$y_lower <- y_lower
    temp_triangle_env$y_upper <- y_upper

    # Extract the list from the environment
    triangle_list <- temp_triangle_env$triangle.output.list

    # Show message
    cli::cli_h1(paste0("Running ", cli::col_cyan("GDVIS plot"), " on ", cli::col_cyan(triangle_list$plot_title)))

    # Directly access the objects
    with(temp_triangle_env, {

    # Extract the variables of the list
    list2env(triangle.output.list, envir = temp_triangle_env)





#### Create dataframes --------------------------------------------------------------------------------


    # Create line and points data
    data <- data.frame(
      x = c(x.con, x.con, x.con, x.sub1, x.allcases, x.sub2, x.sub1, x.sub2, x.con, x.allcases, x.sub1, x.sub2),
      y = c(y.con, y.con, y.con, y.sub1, y.allcases, y.sub2, y.sub1, y.sub2, y.con, y.allcases, y.sub1, y.sub2),
      line_category = c("allcases", "sub2", "sub1", "sub1", "allcases", "sub2", "sub1-sub2", "sub1-sub2", NA, NA, NA, NA),
      label = c(rep(NA, 8), name_con, name_allcases, name_sub1,name_sub2),
      label_color = c(rep(NA, 8), "name_con", NA, "name_sub1","name_sub2"),
      type = c(rep("line", 8),"point", NA, rep("point", 2)),
      linetype = c("dashed", rep("solid", 3), "dashed", rep("solid", 3),rep(NA, 4)),
      vjust = c(rep(NA, 8), 2, 0, 0, 0),
      hjust = c(rep(NA, 8), 0.5, 0, 0, 0),
      nudge_x = c(rep(NA, 8), 0, 0.01, 0.005, 0.015),
      nudge_y = c(rep(NA, 8), 0, 0.01, 0.015, -0),
      stringsAsFactors = FALSE) %>%
      dplyr::mutate(line_category = factor(line_category, levels = c("sub1", "sub2", "sub1-sub2", "allcases")), linetype = factor(linetype, levels = c("solid", "dashed")))
    data_angle <- data %>%
      dplyr::filter(!is.na(label_color))
    data_annot <- data %>%
      dplyr::filter(is.na(type)) %>%
      dplyr::select(x, y, label) %>%
      dplyr::mutate(label_type = "allcases") %>%
      tibble::add_row(x = x.popmean, y = y.popmean, label = "population mean", label_type = "population mean")

    # Adjust for rg = 1
    sub1_line_width <- if(rg_sub1.con_sub2.con == 1) {2.5} else {1}
    sub2_line_width <- if(rg_sub1.con_sub2.con == 1) {1.3} else {1}
    subsub_line_width <- if(rg_sub1.con_sub2.con == 1) {1.3} else {1}
    allcases_line_width <- if(rg_sub1.con_sub2.con == 1) {0.3} else {0.7}

    # Adjust for nonsignificant
    subsub_line <- if (filter1 == 1) {"dotted"} else {"solid"}
    sub1_line <- if (filter2 == 1) {"dotted"} else {"solid"}
    sub2_line <- if (filter3 == 1) {"dotted"} else {"solid"}
    h2_sub2.con <- if(filter3 == 1) {paste0(h2_sub2.con, " (!)")} else {h2_sub2.con}
    h2_sub1.con <- if(filter2 == 1) {paste0(h2_sub1.con, " (!)")} else {h2_sub1.con}
    h2_sub1.sub2 <- if (filter1 == 1) {paste0(h2_sub1.sub2, " (!)")} else {h2_sub1.sub2}

    ### Put coordinates in vectors
    coord_con                       <- c(x.con, y.con, z.con)
    coord_sub1                      <- c(x.sub1, y.sub1,z.sub1)
    coord_sub2                      <- c(x.sub2,y.sub2,z.sub2)
    coord_allcases                  <- c(x.allcases, y.allcases, z.allcases)
    coord_popmean                   <- c(x.popmean, y.popmean, z.popmean)

    # Calculate midpoints for arrow end
    arrow.sub1_line.end <- coord_con + (coord_sub1 - coord_con) / 2
    arrow.sub2_line.end <- coord_con + (coord_sub2 - coord_con) / 2
    arrow.allcases_line.end <- coord_con + (coord_allcases - coord_con) / 2
    arrow.subsub_line.end <- coord_sub2 + (coord_sub1 - coord_sub2) / 2

    # Start arrow
    arrow.sub1_line.start <- coord_con + 0.95*((coord_sub1 - coord_con) / 2)
    arrow.sub2_line.start <- coord_con + 0.95*((coord_sub2 - coord_con) / 2)
    arrow.allcases_line.start <- coord_con + 0.95*((coord_allcases - coord_con) / 2)
    arrow.subsub_line.start <- coord_sub2 + 0.95*((coord_sub1 - coord_sub2) / 2)

    ### Round values
    round_three_decimals <- function(objects) {
      for (i in objects) {assign(i, formatC(get(i,envir = temp_triangle_env), digits = 3, format = "f"),envir = temp_triangle_env)  }}
    round_no_decimals <- function(objects) {
      for (i in objects) {assign(i, formatC(get(i,envir = temp_triangle_env), digits = 0, format = "f"),envir = temp_triangle_env)  }}
    round_three_decimals(c("rg_sub1.con_sub2.con","rg_sub1.con_sub1.sub2","rg_sub2.con_sub1.sub2","h2_sub2.con", "h2_sub1.con","h2_sub1.sub2","h2_allcases.con"))
    round_no_decimals(c("a.deg_sub1.con_sub2.con","a.deg_sub1.con_sub1.sub2","a.deg_sub2.con_sub1.sub2"))



#### plot 2D triangle --------------------------------------------------------------------------------



    plot <- ggplot2::ggplot() +

      # Plot lines
      ggplot2::geom_line(ggplot2::aes(x = x, y = y, group = line_category, color = line_category, linewidth = line_category, linetype = line_category), data = subset(data, type == "line")) +
      ggplot2::scale_color_manual(values = c("sub2" = "#00B050", "sub1" = "#ED7D31", "sub1-sub2" = "#7030A0", "allcases" = "lightgray"),name = "Heritability",
                         # bquote(italic(h)^2 ~ obs[50/50])
                         labels = c("sub2" = h2_sub2.con, "sub1" = h2_sub1.con,"sub1-sub2" = h2_sub1.sub2, "allcases" = h2_allcases.con)) +

      ggplot2::scale_linewidth_manual(values = c("sub2" = sub2_line_width, "sub1" = sub1_line_width, "sub1-sub2" = subsub_line_width, "allcases" = allcases_line_width), guide = "none") +
      ggplot2::scale_linetype_manual(values = c("sub2" = sub2_line, "sub1" = sub1_line, "sub1-sub2" = subsub_line, "allcases" = "solid"), guide = "none") +

      # Plot the subgroup points
      ggplot2::geom_point(ggplot2::aes(x = x, y = y, fill = label_color), size = 2.5, shape = 21, color = "black", data = data_angle) +
      ggplot2::scale_fill_manual(name = expression(italic(r)[g]~"(degrees)"),
                        values = c("name_con" = "black", "name_sub1" = "white",  name_sub2 = "darkgray"),
                        labels = c("name_con" = paste0(rg_sub1.con_sub2.con, " (", a.deg_sub1.con_sub2.con, "°)"),
                                   "name_sub1" = paste0(rg_sub1.con_sub1.sub2, " (", a.deg_sub1.con_sub1.sub2, "°)"),
                                   "name_sub2" = paste0(rg_sub2.con_sub1.sub2," (", a.deg_sub2.con_sub1.sub2, "°)"))) +
      ggplot2::geom_text(data = data_angle,
                ggplot2::aes(x = x, y = y, label = label, vjust = vjust, hjust = hjust),
                nudge_x = data_angle$nudge_x, nudge_y = data_angle$nudge_y) +

      # Plot the other points
      ggplot2::geom_point(ggplot2::aes(x = x, y = y, shape = label_type, size = label_type), data = data_annot) +
      ggplot2::scale_shape_manual(name = "", values = c("allcases" = 15, "population mean" = 8),
                         labels = c("allcases" = name_allcases, "population mean" = "population mean")) +
      ggplot2::scale_size_manual(values = c("allcases" = 1.5, "population mean" = 1), guide = "none") +

      # Clean background
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),   # Remove axes titles
            axis.text.x = ggplot2::element_blank(),  axis.text.y = ggplot2::element_blank(),    # Remove axes text
            axis.ticks = ggplot2::element_blank(),                                     # Remove axes ticks
            panel.grid.major = ggplot2::element_blank(),                               # Remove major gridlines
            panel.grid.minor = ggplot2::element_blank(),                               # Remove minor gridlines
            panel.background = ggplot2::element_rect(fill = "white", colour = NA)) +   # Set background color to white


      # Set title (set title like this so it is always centered)
      ggplot2::annotate("text", x = 0, y = y.con -0.04, label = plot_title, vjust = 2.5, hjust = 0.5, size = 5, fontface = "bold") +

      # Fix legends
      ggplot2::theme(
            #legend.position.inside = c(0.15, 0.5),                                    # Position the legend at the bottom left
            legend.position = c(0.15, 0.5),
            legend.title = ggplot2::element_text(size = 10),                           # Adjusts the size of the legend title
            legend.text = ggplot2::element_text(size = 8),                             # Adjusts the size of the legend text (labels)
            legend.key.size = grid::unit(0.5, "cm")) +                              # Adjusts the size of the keys in the legend
      ggplot2::guides(color = ggplot2::guide_legend(order = 1), fill = ggplot2::guide_legend(order = 2), shape = ggplot2::guide_legend(order = 3)) +

      # Scale axes
      ggplot2::coord_fixed(ratio = 1)

      # Plot the arrows
      if(rg_sub1.con_sub2.con < 1) {
        plot <- plot +
          ggplot2::geom_segment(ggplot2::aes(x = arrow.subsub_line.start[1], y = arrow.subsub_line.start[2], xend = arrow.subsub_line.end[1], yend = arrow.subsub_line.end[2]), color = "#7030A0", linewidth = 0.5, arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "closed")) +
          ggplot2::geom_segment(ggplot2::aes(x = arrow.allcases_line.start[1], y = arrow.allcases_line.start[2], xend = arrow.allcases_line.end[1], yend = arrow.allcases_line.end[2]), color = "lightgray", linewidth = 0.5, arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "closed")) +
          ggplot2::geom_segment(ggplot2::aes(x = arrow.sub1_line.start[1], y = arrow.sub1_line.start[2], xend = arrow.sub1_line.end[1], yend = arrow.sub1_line.end[2]), color = "#ED7D31", linewidth = 0.5, arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "closed")) +
          ggplot2::geom_segment(ggplot2::aes(x = arrow.sub2_line.start[1], y = arrow.sub2_line.start[2], xend = arrow.sub2_line.end[1], yend = arrow.sub2_line.end[2]), color = "#00B050", linewidth = 0.5, arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "closed"))
        }

      # Conditionally set x and y axis limits
      if (!is.null(x_lower) && !is.null(x_upper)) {
        plot <- plot + ggplot2::scale_x_continuous(limits = c(x_lower, x_upper))  }
      if (!is.null(y_lower) && !is.null(y_upper)) {
        plot <- plot + ggplot2::scale_y_continuous(limits = c(y_lower, y_upper))  }

      # Save plot
      ggplot2::ggsave(plot, file = paste0(folder_location,"/",filename, "_2D.triangle_plot.png"), height = 6, width = 12)

      # Return message
      cli::cli_alert_info(paste0(" Plot saved as ", folder_location,"/",filename,"_2D.triangle_plot.png"))
      cli::cli_alert_success("GDVIS plot succesfully finished")

      # Return the plot object
      return(plot)


  }) # End with env

  }  # End function



