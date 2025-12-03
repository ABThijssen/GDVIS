#' Calculate GDVIS parameters
#'
#' This function calculates GDVIS parameters. Input is based on LDSC estimates. h2 needs to be on the observed scale with 50:50 ascertainment.
#' For 2D the input is:
#' h2_sub1.con: the heritability of the GWAS of the subtype1-cases versus controls.
#' h2_se_sub1.con: the standard error of the heritability
#' name_sub1: the name of subtype1-cases
#' N_sub1: the number of individuals in the subtype1-cases group
#' h2_sub2.con: the heritability of the GWAS of the subtype2-cases versus controls.
#' h2_se_sub2.con: the standard error of the heritability
#' name_sub2: the name of subtype2-cases
#' N_sub1: the number of individuals in the subtype2-cases group
#' name_allcases: the name of the pool of the subtype1-cases and subtype2-cases group
#' name_con: the name of the controls
#' rg_sub1.con_sub2.con: the genetic correlation between the GWAS subtype1-cases versus controls and the GWAS subtype2-cases versus controls
#' rg_se_sub1.con_sub2.con: the standard error of the genetic correlation
#' folder_location: the folder that you want the data to be stored at
#' filename: the name of the file for saving
#' pop.prev_case: the population prevalence of the cases
#' There are optional variables that you can also add, this will allow GDVIS to check its calculated values against LDSC values, the result of which can be found in the logfile.
#' However, this is not necessary and you can create the input list without these variables.
#' optional_LDSC_rg_allcases.con_sub1.con: the genetic correlation between the GWAS all cases versus controls and the GWAS subtype1-cases versus controls
#' optional_LDSC_rg_se_allcases.con_sub1.con: the standard error of the genetic correlation
#' optional_LDSC_rg_allcases.con_sub2.con: the genetic correlation between the GWAS all cases versus controls and the GWAS subtype2-cases versus controls
#' optional_LDSC_rg_se_allcases.con_sub2.con: the standard error of the genetic correlation
#' optional_LDSC_h2_sub1.sub2: the heritability of the GWAS of subtype1-cases versus subtype2-cases
#' optional_LDSC_h2_se_sub1.sub2: the standard error of the heritability
#' optional_LDSC_rg_sub1.con_sub1.sub2: the genetic correlation of the GWAS subtype1-cases versus controls and the GWAS subtype1-cases versus subtype2-cases
#' optional_LDSC_rg_se_sub1.con_sub1.sub2: the standard error of the genetic correlation
#' optional_LDSC_rg_sub2.con_sub1.sub2: the genetic correlation of the GWAS subtype2-cases versus controls and the GWAS subtype1-cases versus subtype2-cases
#' optional_LDSC_rg_se_sub2.con_sub1.sub2: the standard error of the genetic correlation
#' optional_LDSC_h2_allcases.con: the heritability of the GWAS of allversus controls
#' optional_LDSC_h2_se_allcases.con: the standard error of the heritability
#'
#' If you want to plot a subtype with an external trait, you need this additional input:
#' plot_3D: this tells GDVIS that you want to calculate a subtype with an external trait and should be set to TRUE
#' h2_ext: the heritabitlity of the external trait
#' h2_se_ext: the standard error of the heritabilty
#' rg_sub1.con_ext: the genetic correlation of the GWAS subtype1-cases versus controls and the GWAS of the external trait
#' rg_se_sub1.con_ext: the standard error of the genetic correlation
#' rg_sub2.con_ext: the genetic correlation of the GWAS subtype2-cases versus controls and the GWAS of the external trait
#' rg_se_sub2.con_ext: the standard error of the genetic correlation
#' name_ext: the name of the external trait
#' pop.prev_ext: the population prevalence of the external trait, if the trait is continuous, set the prevalence to 0.5
#' There is some optional input here as well, will allow GDVIS to check its calculated values against LDSC values:
#' optional_LDSC_rg_sub1.sub2_ext: the genetic correlation of the GWAS subtype1-cases versus subtype2-cases and the GWAS of the external trait
#' optional_LDSC_rg_se_sub1.sub2_ext: the standard error of the genetic correlation
#' optional_LDSC_rg_allcases.con_ext: the genetic correlation of the GWAS all cases and the GWAS of the external trait
#' optional_LDSC_rg_se_allcases.con_ext: the standard error of the genetic correlation
#' If you want to compare two different subtypes, you need the following:
#' all input for both subtypes as described above, but staarting with the triangle1./triangle2. e.g. triangle1.h2_sub1.con and triangle2.h2_sub2.con
#' plot_2D.2D: this tells GDVIS to run in the subtype vs subtype mode
#' rg_triangle1.sub1.sub2_triangle2.sub1.sub2: the genetic correlation of the GWAS subtype1-cases versus subtype2-cases from subtype definition A and the GWAS subtype1-cases versus subtype2-cases from subtype definition B
#' rg_se_triangle1.sub1.sub2_triangle2.sub1.sub2: the standard error the the genetic correlation
#' rg_triangle1.sub1.con_triangle2.sub1.con: the genetic correlation of the GWAS subtype1-cases versus controls from subtype definition A and the GWAS subtype1-cases versus controls from subtype definition B
#' rg_se_triangle1.sub1.con_triangle2.sub1.con: the standard error the the genetic correlation
#' rg_triangle1.sub2.con_triangle2.sub1.con: the genetic correlation of the GWAS subtype12-cases versus controls from subtype definition A and the GWAS subtype1-cases versus controls from subtype definition B
#' rg_se_triangle1.sub2.con_triangle2.sub1.con: the standard error the the genetic correlation
#' rg_triangle1.sub1.con_triangle2.sub2.con: the genetic correlation of the GWAS subtype1-cases versus controls from subtype definition A and the GWAS subtype2-cases versus controls from subtype definition B
#' rg_se_triangle1.sub1.con_triangle2.sub2.con: the standard error the the genetic correlation
#' rg_triangle1.sub2.con_triangle2.sub2.con: the genetic correlation of the GWAS subtype2-cases versus controls from subtype definition A and the GWAS subtype2-cases versus controls from subtype definition B
#' rg_se_triangle1.sub2.con_triangle2.sub2.con: the standard error the the genetic correlation
#' @param triangle.input.list input list with all the parameters
#' @param log_fun needs to stay on message, important for webversion
#' @param webversion needs to stay on FALSE
#' @export
# Function to calculate parameters needed for triangle plotting
GDVIS_calc <- function(triangle.input.list, log_fun = message, webversion = FALSE) {


#### Set up the environment and load + check data ------------------------------------------------------------------------------------


  ## Add plot_3D, plot_2D.2D and plot_CD to the input list if it is not there
  if (!"plot_3D" %in% names(triangle.input.list))  {triangle.input.list$plot_3D <- FALSE }
  if (!"plot_2D.2D" %in% names(triangle.input.list))  {triangle.input.list$plot_2D.2D <- FALSE }
  if (!"plot_CD" %in% names(triangle.input.list))  {triangle.input.list$plot_CD <- FALSE }

  # Add folder_location if not in input list
  if (!"folder_location" %in% names(triangle.input.list)) {triangle.input.list$folder_location <- ""}
  if (!"triangle1.folder_location" %in% names(triangle.input.list)) {triangle.input.list$triangle1.folder_location <- ""}
  if (!"triangle2.folder_location" %in% names(triangle.input.list)) {triangle.input.list$triangle2.folder_location <- ""}

  ## Print the subgroup to the console
  if (triangle.input.list$plot_3D == TRUE & webversion == TRUE) {
    log_fun(paste0("Running GDVIS calc on ", triangle.input.list$plot_title, " with ",  triangle.input.list$name_ext), type = "info") }
  if (triangle.input.list$plot_3D == TRUE & webversion == FALSE) {
    cli::cli_h1(paste0("Running ", cli::col_cyan("GDVIS calc"), " on ", cli::col_cyan(triangle.input.list$plot_title), " with ",  cli::col_cyan(triangle.input.list$name_ext)))  }
  if (triangle.input.list$plot_2D.2D == TRUE  & webversion == FALSE) {
    cli::cli_h1(paste0("Running ", cli::col_cyan("GDVIS calc")," in 2D.2D mode on ", cli::col_cyan(triangle.input.list$triangle1.plot_title),  " and ",cli::col_cyan(triangle.input.list$triangle2.plot_title))) }
  if (triangle.input.list$plot_2D.2D == FALSE & triangle.input.list$plot_3D == FALSE  & webversion == FALSE & triangle.input.list$plot_CD == FALSE) {
    cli::cli_h1(paste0("Running ", cli::col_cyan("GDVIS calc"), " on ", cli::col_cyan(triangle.input.list$plot_title))) }
  if (triangle.input.list$plot_2D.2D == FALSE & triangle.input.list$plot_3D == FALSE  & webversion == TRUE) {
    log_fun(paste0("Running GDVIS calc on ", triangle.input.list$plot_title), type = "info") }
  if (triangle.input.list$plot_CD == TRUE & webversion == FALSE) {
    cli::cli_h1(paste0("Running ", cli::col_cyan("GDVIS calc")," in CD mode on ", cli::col_cyan(triangle.input.list$name_trait1), ", ", cli::col_cyan(triangle.input.list$name_trait2)," and ",cli::col_cyan(triangle.input.list$name_trait3))) }

  # Set required names
  if (triangle.input.list$plot_3D == TRUE) {
    required_pars <- c(
      "h2_sub1.con", "h2_se_sub1.con", "name_sub1", "N_sub1",
      "h2_sub2.con", "h2_se_sub2.con", "name_sub2", "N_sub2",
      "name_allcases", "name_con", "rg_sub1.con_sub2.con", "rg_se_sub1.con_sub2.con",
      "plot_title", "folder_location", "filename", "pop.prev_case",
      "h2_ext" , "h2_se_ext" ,"rg_sub1.con_ext", "rg_se_sub1.con_ext", "rg_sub2.con_ext",
      "rg_se_sub2.con_ext", "name_ext", "pop.prev_ext")

  } else if (triangle.input.list$plot_CD == TRUE) {
    required_pars <- c(
      "name_trait1",  "name_trait2", "name_trait3",  "name_trait1.con", "name_trait2.con", "name_trait3.con", "rg_trait1_trait2",  "rg_trait1_trait3",  "rg_trait2_trait3",
      "h2_trait1","h2_trait2","h2_trait3","pop.prev_trait1","pop.prev_trait2","pop.prev_trait3", "folder_location" )

  } else if (triangle.input.list$plot_2D.2D == TRUE) {
    required_pars <- c(
      "triangle1.h2_sub1.con","triangle1.h2_se_sub1.con", "triangle1.name_sub1","triangle1.N_sub1","triangle1.h2_sub2.con","triangle1.h2_se_sub2.con","triangle1.name_sub2",
      "triangle1.N_sub2", "triangle1.name_allcases","triangle1.name_con", "triangle1.rg_sub1.con_sub2.con", "triangle1.rg_se_sub1.con_sub2.con", #"triangle1.h2_allcases.con","triangle1.h2_se_allcases.con" ,
      "triangle1.plot_title", "triangle1.folder_location","triangle1.filename","triangle1.pop.prev_case",
      "triangle2.h2_sub1.con","triangle2.h2_se_sub1.con", "triangle2.name_sub1","triangle2.N_sub1","triangle2.h2_sub2.con","triangle2.h2_se_sub2.con","triangle2.name_sub2",
      "triangle2.N_sub2", "triangle2.name_allcases","triangle2.name_con", "triangle2.rg_sub1.con_sub2.con", "triangle2.rg_se_sub1.con_sub2.con", #"triangle2.h2_allcases.con","triangle2.h2_se_allcases.con" ,
      "triangle2.plot_title", "triangle2.folder_location","triangle2.filename","triangle2.pop.prev_case",
      "rg_triangle1.sub1.sub2_triangle2.sub1.sub2","rg_se_triangle1.sub1.sub2_triangle2.sub1.sub2","rg_triangle1.sub1.con_triangle2.sub1.con","rg_se_triangle1.sub1.con_triangle2.sub1.con","rg_triangle1.sub2.con_triangle2.sub1.con" ,
      "rg_se_triangle1.sub2.con_triangle2.sub1.con","rg_triangle1.sub1.con_triangle2.sub2.con","rg_se_triangle1.sub1.con_triangle2.sub2.con","rg_triangle1.sub2.con_triangle2.sub2.con","rg_se_triangle1.sub2.con_triangle2.sub2.con"  )
  } else {
    required_pars <- c(
      "h2_sub1.con", "h2_se_sub1.con", "name_sub1", "N_sub1",
      "h2_sub2.con", "h2_se_sub2.con", "name_sub2", "N_sub2",
      "name_allcases", "name_con", "rg_sub1.con_sub2.con", "rg_se_sub1.con_sub2.con",
      "plot_title", "folder_location", "filename", "pop.prev_case" )
    }

  # Check whether all required input is given
  if (webversion == FALSE) {
    missing_pars <- required_pars[!required_pars %in% names(triangle.input.list)]
    if(length(missing_pars) > 0) {
      missing_pars_str <- paste(missing_pars, collapse = ", ")
      error_message <- sprintf("The following parameters are missing: %s", missing_pars_str)
      cli::cli_alert_danger(error_message)
      invokeRestart("abort")  # Stops execution without an error message
    }}
  if (webversion == TRUE) {
    missing_pars <- required_pars[
      !required_pars %in% names(triangle.input.list) |
        sapply(triangle.input.list[required_pars], function(x) is.null(x) || any(is.na(x)))
    ]

    if(length(missing_pars) > 0) {
      missing_pars_str <- paste(missing_pars, collapse = ", ")
      error_message <- sprintf("The following parameters are missing or invalid (if you want to run in 3D, 2D.2D or CD mode make sure to specify that): %s", missing_pars_str)
      log_fun(error_message, type = "danger")
      stop(error_message)
    }
  }

    # Create a new environment
    temp_triangle_env <- new.env()

    # Directly access the objects
    with(temp_triangle_env, {

    ## Set default optional data
    optional_variables <- c("optional_LDSC_h2_sub1.sub2", "optional_LDSC_h2_se_sub1.sub2",  "optional_LDSC_h2_allcases.con", "optional_LDSC_h2_se_allcases.con",
      "optional_LDSC_rg_allcases.con_sub1.con", "optional_LDSC_rg_se_allcases.con_sub1.con", "optional_LDSC_rg_allcases.con_sub2.con", "optional_LDSC_rg_se_allcases.con_sub2.con",
      "optional_LDSC_rg_sub1.con_sub1.sub2", "optional_LDSC_rg_se_sub1.con_sub1.sub2", "optional_LDSC_rg_sub2.con_sub1.sub2", "optional_LDSC_rg_se_sub2.con_sub1.sub2",
      "plot_3D",  "filename3D", "optional_LDSC_rg_sub1.sub2_ext", "optional_LDSC_rg_se_sub1.sub2_ext",
      "optional_LDSC_rg_allcases.con_ext", "optional_LDSC_rg_se_allcases.con_ext", "check5_difference_rg_sub1.sub2_ext_and_optional_LDSC_rg_sub1.sub2_ext", "check5.propSE",
      "check5a.internal.propSE","check5b.internal.propSE","check6_difference_rg_allcases.con_ext_and_optional_LDSC_rg_allcases.con_ext","check6.propSE","x.ext","y.ext","z.ext",
      "x.extcontrol","y.extcontrol","z.extcontrol","coord_ext","a.deg_sub1.con_ext","a.deg_sub2.con_ext","a.deg_sub1.sub2_ext", "rg_sub1.sub2_ext",
      "rg_allcases_con_ext","rg_se_allcases_con_ext","rg_se_sub1_sub2_ext", "rg_sub1.con_ext","rg_se_sub1.con_ext", "rg_sub2.con_ext","rg_se_sub2.con_ext","rg_se_sub1.con_ext","rg_se_sub2.con_ext", "rg_sub1.con_ext.NR", "rg_sub2.con_ext.NR", "name_ext", "pop.prev_ext", "h2_ext", "h2_se_ext", "rg_allcases.con_ext",
      "internal_check_rg_sub1.con_ext", "internal_check_rg_sub2.con_ext", "rg_sub1.con_sub2.con_original", "optional_LDSC_rg_allcases.con_sub1.con_original", "optional_LDSC_rg_allcases.con_sub2.con_original",
      "optional_LDSC_rg_sub1.con_sub1.sub2_original", "optional_LDSC_rg_sub2.con_sub1.sub2_original",
      "optional_LDSC_rg_sub1.sub2_ext_original", "optional_LDSC_rg_allcases.con_ext_original", "rg_sub1.sub2_ext_original", "internal_check_a_sub1.con_ext", "internal_check_a_sub2.con_ext","a_sub1.con_ext", "a_sub2.con_ext",
      "rg_allcases_con_ext_original", "rg_sub1.con_ext_original", "rg_sub2.con_ext_original","rg_allcases.con_ext_original","a.deg_allcases.con_ext","message_triangle_not_possible",
      "a.rad_sub1.sub2_ext", "a.deg_sub1.sub2_ext","a.rad_allcases.con_ext","a.deg_allcases.con_ext", "coordinates_ext", "coordinates_extcontrol",
      "checkmode",
      "triangle1.optional_LDSC_rg_allcases.con_sub1.con","triangle1.optional_LDSC_rg_se_allcases.con_sub1.con","triangle1.optional_LDSC_rg_allcases.con_sub2.con","triangle1.optional_LDSC_rg_se_allcases.con_sub2.con","triangle1.optional_LDSC_h2_sub1.sub2",
      "triangle1.optional_LDSC_h2_se_sub1.sub2","triangle1.optional_LDSC_rg_sub1.con_sub1.sub2","triangle1.optional_LDSC_rg_se_sub1.con_sub1.sub2","triangle1.optional_LDSC_rg_sub2.con_sub1.sub2", "triangle1.optional_LDSC_rg_se_sub2.con_sub1.sub2",
      "triangle2.optional_LDSC_rg_allcases.con_sub1.con","triangle2.optional_LDSC_rg_se_allcases.con_sub1.con","triangle2.optional_LDSC_rg_allcases.con_sub2.con","triangle2.optional_LDSC_rg_se_allcases.con_sub2.con","triangle2.optional_LDSC_h2_sub1.sub2",
      "triangle2.optional_LDSC_h2_se_sub1.sub2","triangle2.optional_LDSC_rg_sub1.con_sub1.sub2","triangle2.optional_LDSC_rg_se_sub1.con_sub1.sub2","triangle2.optional_LDSC_rg_sub2.con_sub1.sub2", "triangle2.optional_LDSC_rg_se_sub2.con_sub1.sub",
      "different_h2", "message_triangle_not_possible_2D.2D")


    # Assign NA to all the optional variables
    for (i in optional_variables) {  assign(i, NA)  }

    # Set to false unless explicitly set to true in input.list
    plot_3D <- FALSE
    plot_2D.2D <- FALSE
    plot_CD <- FALSE

    # Set save data to true as default
    save_data <- TRUE

    # Extract the variables of the list and potentially overwrite default variables
    list2env(triangle.input.list, envir = temp_triangle_env)



#### Create functions --------------------------------------------------------------------------------



    # Function to get line length from two points (Pythagoras: A^2 + B^2 = C^2)
    f.length_from_coords <- function(coord.a, coord.b) {  sqrt(sum((coord.a-coord.b)^2))     }

    # Function to get angle from line lengths  (Cosine law: angle_c = acos( (A^2 + B^2 - C^2) / 2*A*B))
    f.angle_from_lengths <- function(length.a,length.b,length.c) {   acos(((length.a^2 + length.b^2 - length.c^2)) / (2 * length.a * length.b))    }

    # Function to get angle from coordinates of two lines (First get the dir vectors, ngle_c = acos( (A^2 + B^2 - C^2) / 2*A*B),Let everything go the the origin so that the reference is the same)
    f.angle_from_coords <- function(coord.start.line.a, coord.end.line.a, coord.start.line.b, coord.end.line.b) {
      dir.vect.a <- coord.start.line.a - coord.end.line.a
      dir.vect.b <- coord.start.line.b - coord.end.line.b
      acos(  { (f.length_from_coords(c(0,0,0),dir.vect.a)^2) + (f.length_from_coords(c(0,0,0),dir.vect.b)^2) - (f.length_from_coords(dir.vect.a,dir.vect.b)^2) } / {  2 * f.length_from_coords(c(0,0,0),dir.vect.a) * f.length_from_coords(c(0,0,0),dir.vect.b)  }  )   }

    # Function to get line length from two lines and their angle (Cosine law: Length C = A^2 + B^2 - 2*A*B*cos(angle_c))
    f.length_from_angle_and_lines <- function(length.a, length.b, angle.a.b) {
      sqrt(length.a^2 + length.b^2 - 2 * length.a * length.b * cos(angle.a.b)) }

    # Function to extend the line by a certain distance
    f.extend_line <- function(x.coord, y.coord, z.coord, distance) {
      dirvect <- c(0 - x.coord, 0 - y.coord, 0 - z.coord) # Get directional vector for ext
      unitvect <- dirvect / sqrt(sum(dirvect^2))#  Normalize to get unit vector
      coords <- unitvect * distance }# Scale to length

    # Some simple oneliner functions
    rad2deg  <- function(radians)   {  radians * (180 / pi)   }  # Function to convert radians to degrees
    angle2rg <- function(angle)     {  cos(angle)             }  # Function to convert angle in rad to rg
    rg2angle <- function(rg)        {  acos(rg)               }  # Function to convert rg to angle
    d2h2     <- function(distance)  {  distance^2             }  # Function to convert distance to h2
    h22d     <- function(h2)        {  sqrt(h2)               }  # Function to convert h2 to d

    # Function to test whether the two smallest lines together are bigger than the largest line
    f.triangle_possible <- function(length1, length2, length3) {
      sum_smallest_two <- sum(sort(c(length1, length2, length3))[1:2])
      longest <- sort(c(length1, length2, length3))[3]
      if (sum_smallest_two < longest) {
         return(FALSE)  # Return FALSE to indicate invalid triangle
      } else {
        return(TRUE)   # Return TRUE if triangle is possible
      }
    }


## Start code chunk that depends on 2D.2D

    run_code_chunk <- function(triangle = NULL) {

      # This part will only rename the triangles if 2D.2D setting is on
      if (!is.null(triangle)) {

        # Update all variables based on the triangle
        h2_sub1.con                                   <- get(paste0(triangle, ".h2_sub1.con"))
        h2_se_sub1.con                                <- get(paste0(triangle, ".h2_se_sub1.con"))
        h2_sub2.con                                   <- get(paste0(triangle, ".h2_sub2.con"))
        h2_se_sub2.con                                <- get(paste0(triangle, ".h2_se_sub2.con"))
        h2_allcases.con                               <- get(paste0(triangle, ".optional_LDSC_h2_allcases.con")) # This is not an optional input for 2D.2D
        h2_se_allcases.con                            <- get(paste0(triangle, ".optional_LDSC_h2_se_allcases.con"))
        optional_LDSC_h2_sub1.sub2                    <- get(paste0(triangle, ".optional_LDSC_h2_sub1.sub2"))
        optional_LDSC_h2_se_sub1.sub2                 <- get(paste0(triangle, ".optional_LDSC_h2_se_sub1.sub2"))
        rg_sub1.con_sub2.con                          <- get(paste0(triangle, ".rg_sub1.con_sub2.con"))
        rg_se_sub1.con_sub2.con                       <- get(paste0(triangle, ".rg_se_sub1.con_sub2.con"))
        optional_LDSC_rg_sub1.con_sub1.sub2           <- get(paste0(triangle, ".optional_LDSC_rg_sub1.con_sub1.sub2"))
        optional_LDSC_rg_se_sub1.con_sub1.sub2        <- get(paste0(triangle, ".optional_LDSC_rg_se_sub1.con_sub1.sub2"))
        optional_LDSC_rg_sub2.con_sub1.sub2           <- get(paste0(triangle, ".optional_LDSC_rg_sub2.con_sub1.sub2"))
        optional_LDSC_rg_se_sub2.con_sub1.sub2        <- get(paste0(triangle, ".optional_LDSC_rg_se_sub2.con_sub1.sub2"))
        optional_LDSC_rg_allcases.con_sub1.con        <- get(paste0(triangle, ".optional_LDSC_rg_allcases.con_sub1.con"))
        optional_LDSC_rg_se_allcases.con_sub1.con     <- get(paste0(triangle, ".optional_LDSC_rg_se_allcases.con_sub1.con"))
        optional_LDSC_rg_allcases.con_sub2.con        <- get(paste0(triangle, ".optional_LDSC_rg_allcases.con_sub2.con"))
        optional_LDSC_rg_se_allcases.con_sub2.con     <- get(paste0(triangle, ".optional_LDSC_rg_se_allcases.con_sub2.con"))
        name_con                                      <- get(paste0(triangle, ".name_con"))
        name_allcases                                 <- get(paste0(triangle, ".name_allcases"))
        name_sub1                                     <- get(paste0(triangle, ".name_sub1"))
        name_sub2                                     <- get(paste0(triangle, ".name_sub2"))
        N_sub1                                        <- get(paste0(triangle, ".N_sub1"))
        N_sub2                                        <- get(paste0(triangle, ".N_sub2"))
        filename                                      <- get(paste0(triangle, ".filename"))
        folder_location                               <- get(paste0(triangle, ".folder_location"))
        plot_title                                    <- get(paste0(triangle, ".plot_title"))
        pop.prev_case                                 <- get(paste0(triangle, ".pop.prev_case"))

        # Print message start 2D
        if (webversion == FALSE) {
          cli::cli_inform(paste0("Starting 2D calculations for ", cli::col_cyan(plot_title))) }
        if (webversion == TRUE) {
          log_fun(paste0("Starting 2D calculations for ", plot_title), type = "info")   }

      } # End renaming triangles



#### Check if values exceed [-1,1] and restrain if necessary -----------------------------------------------------------------


    # Get the variables that need to be checked
    rg_variables <- c("rg_sub1.con_sub2.con", "optional_LDSC_rg_allcases.con_sub1.con", "optional_LDSC_rg_allcases.con_sub2.con",
                      "optional_LDSC_rg_sub1.con_sub1.sub2", "optional_LDSC_rg_sub2.con_sub1.sub2",
                      "optional_LDSC_rg_sub1.sub2_ext", "optional_LDSC_rg_allcases.con_ext", "rg_sub1.sub2_ext",
                      "rg_allcases_con_ext", "rg_sub1.con_ext", "rg_sub2.con_ext","rg_allcases.con_ext" )

    # Initialize log variable
    log_messages <- " The following values were updated to perform calculations:"

    # Loop through the variables and clip to [-1,1] if needed
    for (var in rg_variables) {
      original <- get(var)  # Retrieve the original values
      if (is.na(original)) {next } # Skip to the next variable in the loop
      updated <- pmin(pmax(original, -1), 1)  # Perform the update if needed

      # Compare and log a message if values are updated
      if (!(original == updated)) {
        log_messages <- c(log_messages, sprintf(paste("Value", get(var), "updated to", updated, "for %s"), var))

         assign(var, as.numeric(updated))  # Save the updated values back to the variable
         assign(paste0(var, "_original"), original) # Keep track of the original variable
    }}

    # Print log messages to the console
    if (length(log_messages) > 1 & webversion == FALSE) { cli::cli_alert_info(paste(log_messages, collapse = "\n")) }
    if (length(log_messages) > 1 & webversion == TRUE) { log_fun(paste(log_messages, collapse = "\n"), type = "info") }



#### Calculate filters -----------------------------------------------------------------------------------------------



    # Filter 1: rg significantly different from 1 and heritabilities nonsignificantly different from each other? (indicative of nonsignificant subgroups)
    if ((rg_sub1.con_sub2.con + 1.96 * rg_se_sub1.con_sub2.con > 1) &  ## rg not significantly different from 1
        ( (h2_sub1.con >= h2_sub2.con - 1.96 * h2_se_sub2.con & h2_sub1.con <= h2_sub2.con + 1.96 * h2_se_sub2.con) |  # h2 sub 1 within CI sub2
          (h2_sub2.con >= h2_sub1.con - 1.96 * h2_se_sub1.con & h2_sub2.con <= h2_sub1.con + 1.96 * h2_se_sub1.con) )   # h2 sub2 within CI sub1
    )
      {
      filter1 <- 1

  save_data <- FALSE
  filter1_error_message <- "Subtypes are not significantly different from each other, no triangle will be created"
   if (webversion == FALSE) { cli::cli_alert_danger(cli::col_red(filter1_error_message)) }
  if (webversion == TRUE) { log_fun(filter1_error_message, type = "danger")  }
   log_fun(filter1_error_message, type = "danger")
  stop(filter1_error_message)
    } else {  filter1 <- 0    }
    outcome_filter_1a <- rg_sub1.con_sub2.con + 1.96 * rg_se_sub1.con_sub2.con                                                       # Save specific outcomes for logfile
    interval_filter_1b <- paste0("(",round(h2_sub2.con - 1.96 * h2_se_sub2.con,3)," - ", round(h2_sub2.con + 1.96 * h2_se_sub2.con,3), ")")
    outcome_filter_1b <- h2_sub1.con >= h2_sub2.con - 1.96 * h2_se_sub2.con & h2_sub1.con <= h2_sub2.con + 1.96 * h2_se_sub2.con     # Save specific outcomes for logfile
    outcome_filter_1c <- h2_sub2.con >= h2_sub1.con - 1.96 * h2_se_sub1.con & h2_sub2.con <= h2_sub1.con + 1.96 * h2_se_sub1.con     # Save specific outcomes for logfile
    interval_filter_1c <- paste0("(",round(h2_sub1.con - 1.96 * h2_se_sub1.con,3)," - ", round(h2_sub1.con + 1.96 * h2_se_sub1.con,3), ")")

    # Filter 2: h2 sub1.con signficant?
    if (h2_sub1.con - 1.96 * h2_se_sub1.con < 0) {
      filter2 <- 1
      save_data <- FALSE
      filter2_error_message <- paste0("h2 ", name_sub1, " is not significantly different from 0, no triangle will be created")
      if (webversion == FALSE) {  cli::cli_alert_danger(cli::col_red(paste0("h2 ", name_sub1, " is not significantly different from 0, no triangle will be created"))) }
      if (webversion == TRUE) { log_fun(filter2_error_message, type = "danger") }
      stop(filter2_error_message)
    } else {  filter2 <- 0}
    outcome_filter_2 <- h2_sub1.con - 1.96 * h2_se_sub1.con

    # Filter 3: h2 sub2.con signficant?
    if (h2_sub2.con - 1.96 * h2_se_sub2.con < 0) {
      filter3 <- 1
      save_data <- FALSE
      filter3_error_message <- paste0("h2 ", name_sub2, " is not significantly different from 0, no triangle will be created")
      if (webversion == FALSE) { cli::cli_alert_danger(cli::col_red(paste0("h2 ", name_sub2, " is not significantly different from 0, no triangle will be created"))) }
      if (webversion == TRUE) { log_fun(filter3_error_message, type = "danger") }
      stop(filter3_error_message)
    }     else {filter3 <- 0}
    outcome_filter_3 <- h2_sub2.con - 1.96 * h2_se_sub2.con                                                                            # Save specific outcomes for logfile



#### Calculate 2D parameters ------------------------------------------------------------------------



    ### Calculate proportion parameters
    prop_sub1                       <- N_sub1 / (N_sub1 + N_sub2)
    k                               <- pop.prev_case[1]

    ### Calculate angles and line lengths of triangle

    # Step 1: Convert input parameters to line lengths and angles
    a.rad_sub1.con_sub2.con         <- rg2angle(rg_sub1.con_sub2.con)
    d_sub1.con                      <- h22d(h2_sub1.con)
    d_sub2.con                      <- h22d(h2_sub2.con)
    a.deg_sub1.con_sub2.con         <- rad2deg(a.rad_sub1.con_sub2.con)

    # Step 2: Calculate subsub line length
    d_sub1.sub2                     <- f.length_from_angle_and_lines(d_sub1.con,d_sub2.con,a.rad_sub1.con_sub2.con)
    h2_sub1.sub2                    <- d2h2(d_sub1.sub2)

    # Step 3: Calculate angles with subsub line
    if (rg_sub1.con_sub2.con == 1) {
     a.rad_sub2.con_sub1.sub2        <- 0
     a.rad_sub1.con_sub1.sub2        <- 0 } else {
    a.rad_sub2.con_sub1.sub2        <- f.angle_from_lengths(d_sub1.sub2,d_sub2.con,d_sub1.con)
    a.rad_sub1.con_sub1.sub2        <- f.angle_from_lengths(d_sub1.sub2,d_sub1.con,d_sub2.con) }
    rg_sub2.con_sub1.sub2           <- angle2rg(a.rad_sub2.con_sub1.sub2)
    rg_sub1.con_sub1.sub2           <- angle2rg(a.rad_sub1.con_sub1.sub2)
    a.deg_sub2.con_sub1.sub2        <- rad2deg(a.rad_sub2.con_sub1.sub2)
    a.deg_sub1.con_sub1.sub2        <- rad2deg(a.rad_sub1.con_sub1.sub2)

    # Step 4: Calculate lengths of line parts of subsub (distance is proportional to subtype N)
    d_allcases.sub2                 <- prop_sub1 * d_sub1.sub2
    d_allcases.sub1                 <- (1-prop_sub1) * d_sub1.sub2

    # Step 5: Calculate allcases line length
    d_allcases.con                  <- f.length_from_angle_and_lines(d_sub2.con,d_allcases.sub2,a.rad_sub2.con_sub1.sub2)
    h2_allcases.con                 <- d2h2(d_allcases.con)

    # Step 6: Calculate angles allcases line with subtypes
    if (rg_sub1.con_sub2.con == 1) {
      a.rad_allcases.con_sub2.con     <- 0
      a.rad_allcases.con_sub1.con     <- 0 } else {
    a.rad_allcases.con_sub2.con     <- f.angle_from_lengths(d_allcases.con,d_sub2.con,d_allcases.sub2)
    a.rad_allcases.con_sub1.con     <- f.angle_from_lengths(d_allcases.con,d_sub1.con,d_allcases.sub1) }
    rg_allcases.con_sub2.con        <- angle2rg(a.rad_allcases.con_sub2.con)
    rg_allcases.con_sub1.con        <- angle2rg(a.rad_allcases.con_sub1.con)
    a.deg_allcases.con_sub2.con     <- rad2deg(a.rad_allcases.con_sub2.con)
    a.deg_allcases.con_sub1.con     <- rad2deg(a.rad_allcases.con_sub1.con)

    ### Calculate coordinates of subgroups
    x.con                           <- 0
    y.con                           <- -k * d_allcases.con
    z.con                           <- 0
    x.sub2                          <- d_sub2.con * sin(a.rad_allcases.con_sub2.con) + x.con
    y.sub2                          <- d_sub2.con * cos(a.rad_allcases.con_sub2.con) + y.con
    z.sub2                          <- 0 + z.con
    x.sub1                          <- -d_sub1.con * sin(a.rad_allcases.con_sub1.con) + x.con
    y.sub1                          <- d_sub1.con * cos(a.rad_allcases.con_sub1.con) + y.con
    z.sub1                          <- 0 + z.con
    x.allcases                      <- 0
    y.allcases                      <- y.sub2 - (y.sub1 - y.sub2) / (x.sub1 - x.sub2) * x.sub2       # Where sub1.sub2 crosses the y-axis
    z.allcases                      <- 0
    x.popmean                       <- 0
    y.popmean                       <- 0
    z.popmean                       <- 0

    ### Put coordinates in vectors
    coord_con                       <- c(x.con, y.con, z.con)
    coord_sub1                      <- c(x.sub1, y.sub1,z.sub1)
    coord_sub2                      <- c(x.sub2,y.sub2,z.sub2)
    coord_allcases                  <- c(x.allcases, y.allcases, z.allcases)
    coord_popmean                   <- c(x.popmean, y.popmean, z.popmean)

    ### Update angles using coordinate system, so that directionality is taken into account
    if (rg_sub1.con_sub2.con != 1) {
    a.rad_sub2.con_sub1.sub2        <- f.angle_from_coords(coord_con,coord_sub2,coord_sub2,coord_sub1)
    a.rad_sub1.con_sub1.sub2        <- f.angle_from_coords(coord_sub1,coord_con,coord_sub1,coord_sub2)
    rg_sub2.con_sub1.sub2           <- angle2rg(a.rad_sub2.con_sub1.sub2)
    rg_sub1.con_sub1.sub2           <- angle2rg(a.rad_sub1.con_sub1.sub2)
    a.deg_sub2.con_sub1.sub2        <- rad2deg(a.rad_sub2.con_sub1.sub2)
    a.deg_sub1.con_sub1.sub2        <- rad2deg(a.rad_sub1.con_sub1.sub2) }

    # Change coordinates of subgroups if rg = 1
    if (rg_sub1.con_sub2.con == 1) {

      # First need to redo some calculations while forcing the angles to be zero, as the lines will all be drawn on x = 0

      # Step 3: Calculate angles with subsub line --> force to be zero
      a.rad_sub2.con_sub1.sub2        <- 0
      a.rad_sub1.con_sub1.sub2        <- 0
      rg_sub2.con_sub1.sub2           <- angle2rg(a.rad_sub2.con_sub1.sub2)
      rg_sub1.con_sub1.sub2           <- angle2rg(a.rad_sub1.con_sub1.sub2)
      a.deg_sub2.con_sub1.sub2        <- rad2deg(a.rad_sub2.con_sub1.sub2)
      a.deg_sub1.con_sub1.sub2        <- rad2deg(a.rad_sub1.con_sub1.sub2)

      # Step 5: Calculate allcases line length
      d_allcases.con                  <- d_sub1.con - d_allcases.sub1
      h2_allcases.con                 <- d2h2(d_allcases.con)

      # Step 6: Calculate angles allcases line with subtypes
      a.rad_allcases.con_sub2.con     <- 0
      a.rad_allcases.con_sub1.con     <- 0
      rg_allcases.con_sub2.con        <- angle2rg(a.rad_allcases.con_sub2.con)
      rg_allcases.con_sub1.con        <- angle2rg(a.rad_allcases.con_sub1.con)
      a.deg_allcases.con_sub2.con     <- rad2deg(a.rad_allcases.con_sub2.con)
      a.deg_allcases.con_sub1.con     <- rad2deg(a.rad_allcases.con_sub1.con)

      ### Calculate coordinates of subgroups
      x.con                           <- 0
      y.con                           <- -k * d_allcases.con
      z.con                           <- 0
      x.sub2                          <- 0
      y.sub2                          <- d_sub2.con + y.con
      z.sub2                          <- 0 + z.con
      x.sub1                          <- 0
      y.sub1                          <- d_sub1.con + y.con
      z.sub1                          <- 0 + z.con
      x.allcases                      <- 0
      y.allcases                      <- d_allcases.con + y.con
      z.allcases                      <- 0
      x.popmean                       <- 0
      y.popmean                       <- 0
      z.popmean                       <- 0

      ### Put coordinates in vectors
      coord_con                       <- c(x.con, y.con, z.con)
      coord_sub1                      <- c(x.sub1, y.sub1,z.sub1)
      coord_sub2                      <- c(x.sub2,y.sub2,z.sub2)
      coord_allcases                  <- c(x.allcases, y.allcases, z.allcases)
      coord_popmean                   <- c(x.popmean, y.popmean, z.popmean)
    }

    ### Internal 2D checks
    if( (f.length_from_coords(coord_sub1,coord_con)     -  d_sub1.con )      >  (0.02 * d_sub1.con) )      {stop("Error code: length check d_sub1.con not passed")}
    if( (f.length_from_coords(coord_sub2,coord_con)     -  d_sub2.con )      >  (0.02 * d_sub2.con) )      {stop("Error code: length check d_sub2.con not passed")}
    if( (f.length_from_coords(coord_sub1,coord_sub2)    -  d_sub1.sub2 )     >  (0.02 * d_sub1.sub2) )     {stop("Error code: length check d_sub1.sub2 not passed")}
    if( (f.length_from_coords(coord_allcases,coord_con) -  d_allcases.con )  >  (0.02 * d_allcases.con) )  {stop("Error code: length check d_allcases.con not passed")}
    # Change small value to 0 so that internal check is correctly passed
    a.rad_sub1.con_sub2.con_from_coords <- f.angle_from_coords(coord_con,coord_sub1,coord_con,coord_sub2)
    if(a.rad_sub1.con_sub2.con_from_coords < 0.000001) {a.rad_sub1.con_sub2.con_from_coords <- 0}
    if( (a.rad_sub1.con_sub2.con_from_coords  - a.rad_sub1.con_sub2.con ) > (0.02 * a.rad_sub1.con_sub2.con) )   {stop("Error code: angle check a.rad_sub1.con_sub2.con not passed")}



#### Calculate 3D parameters -----------------------------------------------------------------------------------------------



      if (plot_3D == TRUE) {

      # Make 3D filename
      filename3D  <- paste0(filename, ".", name_ext)

      ## Step 1: Convert input parameters to line lengths and angles
      a_sub1.con_ext         <- rg2angle(rg_sub1.con_ext)
      a_sub2.con_ext         <- rg2angle(rg_sub2.con_ext)
      a_sub1.sub2_ext        <- rg2angle(optional_LDSC_rg_sub1.sub2_ext)
      a.deg_sub1.con_ext     <- rad2deg(a_sub1.con_ext)
      a.deg_sub2.con_ext     <- rad2deg(a_sub2.con_ext)
      a.deg_sub1.sub2_ext    <- rad2deg(a_sub1.sub2_ext)
      k.ext                  <- pop.prev_ext
      d_ext.total            <- h22d(h2_ext)
      d_extcontrol.popmean   <- k.ext * d_ext.total
      d_ext.popmean          <- (1-k.ext) * d_ext.total

      ## Step 2: Calculate helper line lengths from ext to both subtypes
      d_sub1_ext             <- f.length_from_angle_and_lines(d_ext.popmean,d_sub1.con,a_sub1.con_ext)
      d_sub2_ext             <- f.length_from_angle_and_lines(d_ext.popmean,d_sub2.con,a_sub2.con_ext)

      ## Extra step: Check if the triangle is possible with a helper triangle
      # Get helper triangle
      line1      <- f.length_from_angle_and_lines(1,1,a_sub1.con_ext)
      line2      <- f.length_from_angle_and_lines(1,1,a_sub2.con_ext)
      line3      <- f.length_from_angle_and_lines(1,1,a.rad_sub1.con_sub2.con)

      # Update save_data with the function result
      save_data <- f.triangle_possible(line1, line2, line3)
      message_triangle_not_possible <- "GDVIS determined that the triangle is not possible, possibly due to inconsistency in input parameters. Please double-check your input. If the issue persists, it may be due to inconsistency in LDSC estimates, for instance because the rg between sub1.con and sub2.con is close to 1"
      if(f.triangle_possible(line1, line2, line3) == FALSE & webversion == FALSE) { cli::cli_alert_danger(cli::col_red(message_triangle_not_possible))   }
      if(f.triangle_possible(line1, line2, line3) == FALSE & webversion == TRUE) { log_fun(message_triangle_not_possible, type = "info")   }

      # Only check the coordinates if triangle is possible
      if (save_data) {

      ## Print start computing coordinates
      if (webversion == FALSE) { cli::cli_alert_info(" Calculating 3D coordinates, this may take some time ... "  ) }

      ## Step 3: Find ext coordinates using brute force

      # Function to find best fitting external trait coordinates in the matrix
      test_matrix_fit <- function(xyz_matrix) {

        # Function to find best fitting external trait coordinates
        test_row_fit <- function(x.ext,y.ext,z.ext){
          fit_sub1   <- abs(  {d_sub1_ext^2}  - { {(x.ext-x.sub1)^2} + {(y.ext-y.sub1)^2} + {(z.ext-z.sub1)^2}  }   )
          fit_sub2   <- abs(  {d_sub2_ext^2}  - { {(x.ext-x.sub2)^2} + {(y.ext-y.sub2)^2} + {(z.ext-z.sub2)^2}  }   )
          fit_con    <- abs(  {d_ext.popmean^2}       - { {(x.ext-x.con)^2}  + {(y.ext-y.con)^2}  + {(z.ext-z.con)^2}   }   )
          return(fit_sub1+fit_sub2+fit_con)}

        # Find the differences for each combination of x,y,z and add it to the matrix
        xyz_fit <- tibble::as_tibble(xyz_matrix) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(fit = test_row_fit(x, y, z))

        # Get the lowest fit
        min_fit <- min(xyz_fit$fit)
        best_fit <- xyz_fit %>%
          dplyr::filter(fit == min_fit)

        return(best_fit)      }

      # Function to generate matrix and find best fit
      generate_best_fit <- function(xs, ys, zs) {
        xyz_matrix <- as.matrix(tidyr::crossing(x = xs, y = ys, z = zs))
        best_fit <- test_matrix_fit(xyz_matrix) #%>% unlist()
        return(best_fit)}

      # Initial rough matrix
      if(rg_sub1.con_sub2.con == 1) {
        coords_rough <- generate_best_fit(0,
                                          seq(from=-1, to=1, by=0.01),
                                          seq(from=0, to=1, by=0.01))
      } else{

      coords_rough <- generate_best_fit(seq(from=-1, to=1, by=0.01),
                                        seq(from=-1, to=1, by=0.01),
                                        seq(from=0, to=1, by=0.01))      }

      # Finer matrix based on rough matrix (was 0.01)
      coords_finer <- generate_best_fit(seq(from = coords_rough['x'] %>% dplyr::pull() - 0.04, to = coords_rough['x'] %>% dplyr::pull() + 0.04, by = 0.001),
                                        seq(from = coords_rough['y'] %>% dplyr::pull() - 0.04, to = coords_rough['y'] %>% dplyr::pull() + 0.04, by = 0.001),
                                        seq(from = max(0, coords_rough['z'] %>% dplyr::pull() - 0.04), to = coords_rough['z'] %>% dplyr::pull() + 0.04, by = 0.001) )

      # Coords finest if there are multiple peaks of best fit
       coords_finest <- tibble::tibble()

       for (i in 1:nrow(coords_finer)) {

         # Get one row
         line <- coords_finer[i,]

         # Test the row
         out <- generate_best_fit(seq(from=line['x'] %>% dplyr::pull() - 0.004, to=line['x'] %>% dplyr::pull() + 0.004, by=0.0001),
                                  seq(from=line['y'] %>% dplyr::pull() - 0.004, to=line['y'] %>% dplyr::pull() + 0.004, by=0.0001),
                                  seq(from=line['z'] %>% dplyr::pull() - 0.004, to=line['z'] %>% dplyr::pull() + 0.004, by=0.0001))

         # Add fit to results
         coords_finest <- rbind(coords_finest, out)
       } # end for loop

        # Get the lowest fit
        min_fit <- min(coords_finest$fit)

        # Filter on lowest fit
        coords_finest <- coords_finest %>%
          dplyr::ungroup() %>%  # Remove rowwise structure
          dplyr::filter(fit == min_fit)

        # Lowest fit will sometimes still include two options, then just get the first one
        coords_finest <- coords_finest[1,]

    #    } # End if multiple peaks


      # Finest matrix based on finer matrix (was 0.002)
      # coords_finest <- generate_best_fit(seq(from=coords_finer['x'] %>% dplyr::pull() - 0.004, to=coords_finer['x'] %>% dplyr::pull() + 0.004, by=0.0001),
      #                                    seq(from=coords_finer['y'] %>% dplyr::pull() - 0.004, to=coords_finer['y'] %>% dplyr::pull() + 0.004, by=0.0001),
      #                                    seq(from=coords_finer['z'] %>% dplyr::pull() - 0.004, to=coords_finer['z'] %>% dplyr::pull() + 0.004, by=0.0001))

      # Finest matrix based on finer matrix (was 0.002)
      coords_finestest <- generate_best_fit(seq(from=coords_finest['x'] %>% dplyr::pull() - 0.0004, to=coords_finest['x'] %>% dplyr::pull() + 0.0004, by=0.00001),
                                         seq(from=coords_finest['y'] %>% dplyr::pull() - 0.0004, to=coords_finest['y'] %>% dplyr::pull() + 0.0004, by=0.00001),
                                         seq(from=coords_finest['z'] %>% dplyr::pull() - 0.0004, to=coords_finest['z'] %>% dplyr::pull() + 0.0004, by=0.00001))

      # Coordinates
      #coords_matrix <- rbind(coords_rough, coords_finer, coords_finest, coords_finestest)

      # Extract coordinates
      x.ext <- coords_finestest %>% dplyr::select(x) %>% dplyr::pull()
      y.ext <- coords_finestest %>% dplyr::select(y) %>% dplyr::pull()
      z.ext <- coords_finestest %>% dplyr::select(z) %>% dplyr::pull()

      ## Step 4: Move the line to popmean (x and z coordinates stay the same)
      y.ext  <- y.ext + abs(y.con)
      coord_ext <- c(x.ext, y.ext, z.ext)

      ## Step 5: Get coordinates of ext line to ext controls

      # Get directional vector for ext
      dirvect.ext <- c(0 - x.ext, 0 - y.ext, 0 - z.ext)

      # Normalize to get unit vector
      unitvect.ext <- dirvect.ext / sqrt(sum(dirvect.ext^2))

      # Scale to length
      coord_ext.controls <- unitvect.ext * d_extcontrol.popmean

      # Extract coordinates
      x.extcontrol <- coord_ext.controls[1]
      y.extcontrol <- coord_ext.controls[2]
      z.extcontrol <- coord_ext.controls[3]
      coord_ext.con <- c(x.extcontrol, y.extcontrol,z.extcontrol)

      } # End if save_data == TRUE

    }


#### Calculate checks ---------------------------------------------------------------------------------------------------



    ## Convert rg to angles
    optional_LDSC_a.rad_sub1.con_allcases.con                         <- acos(optional_LDSC_rg_allcases.con_sub1.con)
    optional_LDSC_a.rad_sub2.con_allcases.con                         <- acos(optional_LDSC_rg_allcases.con_sub2.con)

    ## Convert from rad to degrees
    optional_LDSC_a.deg_sub1.con_allcases.con                         <- rad2deg(optional_LDSC_a.rad_sub1.con_allcases.con)
    optional_LDSC_a.deg_sub2.con_allcases.con                         <- rad2deg(optional_LDSC_a.rad_sub2.con_allcases.con)
    a.deg_sub1.con_sub2.con                                           <- rad2deg(a.rad_sub1.con_sub2.con)

    # Function to check whether a value falls within the SE
    fallswithin <- function(new = new, old = old, se = se) {
      if(is.na(old) | is.na(se)) {
        outcome <- NA
      return(outcome)  }
      lower = old - 1.96 * se
      upper = old + 1.96 * se
      if (new > lower & new < upper) {outcome <- TRUE}
      else {outcome <- FALSE}
      return(outcome)  }

    # Check 1: subsub
    check1_difference_h2_sub1.sub2_and_optional_LDSC_h2_sub1.sub2                            <- h2_sub1.sub2 - optional_LDSC_h2_sub1.sub2
    check1_passed_h2_sub1.sub2_and_optional_LDSC_h2_sub1.sub2                                <- fallswithin(h2_sub1.sub2,optional_LDSC_h2_sub1.sub2,optional_LDSC_h2_se_sub1.sub2)
    check1.propSE                                                                            <- check1_difference_h2_sub1.sub2_and_optional_LDSC_h2_sub1.sub2 / optional_LDSC_h2_se_sub1.sub2

    # Check 2: allcases
    check2_difference_h2_allcases.con_and_optional_LDSC_h2_allcases.con                      <- h2_allcases.con - optional_LDSC_h2_allcases.con
    check2_passed_h2_allcases.con_and_optional_LDSC_h2_allcases.con                          <- fallswithin(h2_allcases.con, optional_LDSC_h2_allcases.con, optional_LDSC_h2_se_allcases.con)
    check2.propSE                                                                            <- check2_difference_h2_allcases.con_and_optional_LDSC_h2_allcases.con / optional_LDSC_h2_se_allcases.con

    # Check 3: top corners
    check3a_difference_rg_sub1.con_sub1.sub2_and_optional_LDSC_rg_sub1.con_sub1.sub2         <- rg_sub1.con_sub1.sub2 - optional_LDSC_rg_sub1.con_sub1.sub2
    check3b_difference_rg_sub2.con_sub1.sub2_and_optional_LDSC_rg_sub2.con_sub1.sub2         <- rg_sub2.con_sub1.sub2 - optional_LDSC_rg_sub2.con_sub1.sub2
    check3a_passed_rg_sub1.con_sub1.sub2_and_optional_LDSC_rg_sub1.con_sub1.sub2             <- fallswithin(rg_sub1.con_sub1.sub2,optional_LDSC_rg_sub1.con_sub1.sub2,optional_LDSC_rg_se_sub1.con_sub1.sub2)
    check3b_passed_rg_sub2.con_sub1.sub2_and_optional_LDSC_rg_sub2.con_sub1.sub2             <- fallswithin(rg_sub2.con_sub1.sub2,optional_LDSC_rg_sub2.con_sub1.sub2,optional_LDSC_rg_se_sub2.con_sub1.sub2)
    check3a.propSE                                                                           <- check3a_difference_rg_sub1.con_sub1.sub2_and_optional_LDSC_rg_sub1.con_sub1.sub2 / optional_LDSC_rg_se_sub1.con_sub1.sub2
    check3b.propSE                                                                           <- check3b_difference_rg_sub2.con_sub1.sub2_and_optional_LDSC_rg_sub2.con_sub1.sub2 / optional_LDSC_rg_se_sub2.con_sub1.sub2

    # Check 4: allcases with subgroups
    check4a_difference_rg_sub1.con_allcases.con_and_optional_LDSC_rg_sub1.con_allcases.con   <- rg_allcases.con_sub1.con - optional_LDSC_rg_allcases.con_sub1.con
    check4b_difference_rg_sub2.con_allcases.con_and_optional_LDSC_rg_sub2.con_allcases.con   <- rg_allcases.con_sub2.con - optional_LDSC_rg_allcases.con_sub2.con
    check4a_passed_rg_sub1.con_allcases.con_and_optional_LDSC_rg_sub1.con_allcases.con       <- fallswithin(rg_allcases.con_sub1.con,optional_LDSC_rg_allcases.con_sub1.con,optional_LDSC_rg_se_allcases.con_sub1.con)
    check4b_passed_rg_sub2.con_allcases.con_and_optional_LDSC_rg_sub2.con_allcases.con       <- fallswithin(rg_allcases.con_sub2.con,optional_LDSC_rg_allcases.con_sub2.con,optional_LDSC_rg_se_allcases.con_sub2.con)
    check4a.propSE                                                                           <- check4a_difference_rg_sub1.con_allcases.con_and_optional_LDSC_rg_sub1.con_allcases.con / optional_LDSC_rg_se_allcases.con_sub1.con
    check4b.propSE                                                                           <- check4b_difference_rg_sub2.con_allcases.con_and_optional_LDSC_rg_sub2.con_allcases.con / optional_LDSC_rg_se_allcases.con_sub2.con

    # 3D checks
    if (plot_3D == TRUE & save_data == TRUE) {
    a.rad_sub1.sub2_ext                                                                      <- f.angle_from_coords(coord.start.line.a = coord_popmean, coord.end.line.a = coord_ext, coord.start.line.b = coord_sub2, coord.end.line.b = coord_sub1)
    a.deg_sub1.sub2_ext                                                                      <- rad2deg(a.rad_sub1.sub2_ext)
    rg_sub1.sub2_ext                                                                         <- angle2rg(a.rad_sub1.sub2_ext)
    check5_difference_rg_sub1.sub2_ext_and_optional_LDSC_rg_sub1.sub2_ext                    <- rg_sub1.sub2_ext -  optional_LDSC_rg_sub1.sub2_ext
    check5_passed_rg_sub1.sub2_ext_and_optional_LDSC_rg_sub1.sub2_ext                        <- fallswithin(rg_sub1.sub2_ext,optional_LDSC_rg_sub1.sub2_ext,optional_LDSC_rg_se_sub1.sub2_ext)
    check5.propSE                                                                            <- check5_difference_rg_sub1.sub2_ext_and_optional_LDSC_rg_sub1.sub2_ext / optional_LDSC_rg_se_sub1.sub2_ext

    a.rad_allcases.con_ext                                                                   <- f.angle_from_coords(coord.start.line.a = coord_popmean, coord.end.line.a = coord_ext, coord.start.line.b = coord_con, coord.end.line.b = coord_allcases)
    a.deg_allcases.con_ext                                                                   <- rad2deg(a.rad_allcases.con_ext)
    rg_allcases.con_ext                                                                      <- angle2rg(a.rad_allcases.con_ext)
    check6_difference_rg_allcases.con_ext_and_optional_LDSC_rg_allcases.con_ext              <- rg_allcases.con_ext -  optional_LDSC_rg_allcases.con_ext
    check6_passed_rg_allcases.con_ext_and_optional_LDSC_rg_allcases.con_ext                  <- fallswithin(rg_allcases.con_ext,optional_LDSC_rg_allcases.con_ext,optional_LDSC_rg_se_allcases.con_ext)
    check6.propSE                                                                            <- check6_difference_rg_allcases.con_ext_and_optional_LDSC_rg_allcases.con_ext / optional_LDSC_rg_se_allcases.con_ext

    # Internal 3D checks
    internal_check_a_sub1.con_ext                                                            <- f.angle_from_coords(coord.start.line.a = coord_popmean, coord.end.line.a = coord_ext, coord.start.line.b = coord_con, coord.end.line.b = coord_sub1)
    internal_check_rg_sub1.con_ext                                                           <- angle2rg(internal_check_a_sub1.con_ext)
    check5a.internal.propSE                                                                  <- (internal_check_rg_sub1.con_ext - rg_sub1.con_ext) / rg_se_sub1.con_ext
    internal_check_a_sub2.con_ext                                                            <- f.angle_from_coords(coord.start.line.a = coord_popmean, coord.end.line.a = coord_ext, coord.start.line.b = coord_con, coord.end.line.b = coord_sub2)
    internal_check_rg_sub2.con_ext                                                           <- angle2rg(internal_check_a_sub2.con_ext)
    check5b.internal.propSE                                                                  <- (internal_check_rg_sub2.con_ext - rg_sub2.con_ext) / rg_se_sub2.con_ext

    if(save_data) {

      # Calculate checks
      if( abs((f.length_from_coords(coord_ext,coord_popmean)     -  d_ext.popmean ))          >  ( 0.02 * d_ext.popmean ))  {
        if (webversion == FALSE) {cli::cli_alert_danger(cli::col_red("GDVIS was not able to accurately estimate 3D coordinates, possibly due to inconsistency in input parameters. Please double-check your input. If the issue persists, it may be due to inconsistency in LDSC estimates, for instance because the rg between sub1.con and sub2.con is close to 1. Specific error: calculated d_ext.popmean differs more than 2% from LDSC value")) }
        if (webversion == TRUE) { log_fun("GDVIS was not able to accurately estimate 3D coordinates, possibly due to inconsistency in input parameters. Please double-check your input. If the issue persists, it may be due to inconsistency in LDSC estimates, for instance because the rg between sub1.con and sub2.con is close to 1. Specific error: calculated d_ext.popmean differs more than 2% from LDSC value", type = "info") }
        save_data <- FALSE }
      if( abs((f.length_from_coords(coord_ext.con,coord_popmean) -  d_extcontrol.popmean ))   >  ( 0.02 * d_extcontrol.popmean ))  {
        if (webversion == FALSE) { cli::cli_alert_danger(cli::col_red("GDVIS was not able to accurately estimate 3D coordinates, possibly due to inconsistency in input parameters. Please double-check your input. If the issue persists, it may be due to inconsistency in LDSC estimates, for instance because the rg between sub1.con and sub2.con is close to 1. Specific error: calculated d_extcontrol.popmean differs more than 2% from LDSC value")) }
        if (webversion == TRUE) {log_fun("GDVIS was not able to accurately estimate 3D coordinates, possibly due to inconsistency in input parameters. Please double-check your input. If the issue persists, it may be due to inconsistency in LDSC estimates, for instance because the rg between sub1.con and sub2.con is close to 1. Specific error: calculated d_extcontrol.popmean differs more than 2% from LDSC value", type = "info") }
        save_data <- FALSE }
      if( abs((internal_check_a_sub1.con_ext  -  a_sub1.con_ext)) > 0.05) {
        if (webversion == FALSE) { cli::cli_alert_danger(cli::col_red("GDVIS was not able to accurately estimate 3D coordinates, possibly due to inconsistency in input parameters. Please double-check your input. If the issue persists, it may be due to inconsistency in LDSC estimates, for instance because the rg between sub1.con and sub2.con is close to 1. Specific error: calculated angle a_sub1.con_ext differs more than 0.05 from LDSC value")) }
        if (webversion == TRUE) {log_fun("GDVIS was not able to accurately estimate 3D coordinates, possibly due to inconsistency in input parameters. Please double-check your input. If the issue persists, it may be due to inconsistency in LDSC estimates, for instance because the rg between sub1.con and sub2.con is close to 1. Specific error: calculated angle a_sub1.con_ext differs more than 0.05 from LDSC value", type = "info") }
        save_data <- FALSE }
      if( abs((internal_check_a_sub2.con_ext  -  a_sub2.con_ext)) > 0.05) {
        if (webversion == FALSE) { cli::cli_alert_danger(cli::col_red("GDVIS was not able to accurately estimate 3D coordinates, possibly due to inconsistency in input parameters. Please double-check your input. If the issue persists, it may be due to inconsistency in LDSC estimates, for instance because the rg between sub1.con and sub2.con is close to 1. Specific error: calculated angle a_sub2.con_ext differs more than 0.05 from LDSC value")) }
        if (webversion == TRUE) {log_fun("GDVIS was not able to accurately estimate 3D coordinates, possibly due to inconsistency in input parameters. Please double-check your input. If the issue persists, it may be due to inconsistency in LDSC estimates, for instance because the rg between sub1.con and sub2.con is close to 1. Specific error: calculated angle a_sub2.con_ext differs more than 0.05 from LDSC value", type = "info") }
        save_data <- FALSE }
    }
    }



#### Save parameters in output list -------------------------------------------------------------------------------------------------



    triangle.output.list <- list(
      h2_sub1.con                                                          = h2_sub1.con,
      h2_se_sub1.con                                                       = h2_se_sub1.con,
      h2_sub2.con                                                          = h2_sub2.con,
      h2_se_sub2.con                                                       = h2_se_sub2.con,
      h2_sub1.sub2                                                         = h2_sub1.sub2,
      h2_allcases.con                                                      = h2_allcases.con,
      optional_LDSC_h2_sub1.sub2                                           = optional_LDSC_h2_sub1.sub2,
      optional_LDSC_h2_se_sub1.sub2                                        = optional_LDSC_h2_se_sub1.sub2,
      optional_LDSC_h2_allcases.con                                        = optional_LDSC_h2_allcases.con,
      optional_LDSC_h2_se_allcases.con                                     = optional_LDSC_h2_se_allcases.con,
      d_sub1.con                                                           = d_sub1.con,
      d_sub2.con                                                           = d_sub2.con,

      rg_sub1.con_sub2.con                                                 = rg_sub1.con_sub2.con,
      rg_se_sub1.con_sub2.con                                              = rg_se_sub1.con_sub2.con,
      rg_sub1.con_sub1.sub2                                                = rg_sub1.con_sub1.sub2,
      optional_LDSC_rg_sub1.con_sub1.sub2                                  = optional_LDSC_rg_sub1.con_sub1.sub2,
      optional_LDSC_rg_se_sub1.con_sub1.sub2                               = optional_LDSC_rg_se_sub1.con_sub1.sub2,
      rg_sub2.con_sub1.sub2                                                = rg_sub2.con_sub1.sub2,
      optional_LDSC_rg_sub2.con_sub1.sub2                                  = optional_LDSC_rg_sub2.con_sub1.sub2,
      optional_LDSC_rg_se_sub2.con_sub1.sub2                               = optional_LDSC_rg_se_sub2.con_sub1.sub2,
      rg_allcases.con_sub1.con                                             = rg_allcases.con_sub1.con,
      optional_LDSC_rg_allcases.con_sub1.con                               = optional_LDSC_rg_allcases.con_sub1.con,
      optional_LDSC_rg_se_allcases.con_sub1.con                            = optional_LDSC_rg_se_allcases.con_sub1.con,
      rg_allcases.con_sub2.con                                             = rg_allcases.con_sub2.con,
      optional_LDSC_rg_allcases.con_sub2.con                               = optional_LDSC_rg_allcases.con_sub2.con,
      optional_LDSC_rg_se_allcases.con_sub2.con                            = optional_LDSC_rg_se_allcases.con_sub2.con,

      a.rad_sub1.con_sub2.con                                              = a.rad_sub1.con_sub2.con,
      a.rad_sub1.con_sub1.sub2                                             = a.rad_sub1.con_sub1.sub2,
      a.rad_sub2.con_sub1.sub2                                             = a.rad_sub2.con_sub1.sub2,
      optional_LDSC_a.rad_sub1.con_allcases.con                            = optional_LDSC_a.rad_sub1.con_allcases.con,
      optional_LDSC_a.rad_sub2.con_allcases.con                            = optional_LDSC_a.rad_sub2.con_allcases.con,

      a.deg_sub1.con_sub2.con                                              = a.deg_sub1.con_sub2.con,
      a.deg_sub1.con_sub1.sub2                                             = a.deg_sub1.con_sub1.sub2,
      a.deg_sub2.con_sub1.sub2                                             = a.deg_sub2.con_sub1.sub2,
      optional_LDSC_a.deg_sub1.con_allcases.con                            = optional_LDSC_a.deg_sub1.con_allcases.con,
      optional_LDSC_a.deg_sub2.con_allcases.con                            = optional_LDSC_a.deg_sub2.con_allcases.con,


      check1_difference_h2_sub1.sub2_and_optional_LDSC_h2_sub1.sub2                             = check1_difference_h2_sub1.sub2_and_optional_LDSC_h2_sub1.sub2,
      check1.propSE                                                                             = check1.propSE,
      check2_difference_h2_allcases.con_and_optional_LDSC_h2_allcases.con                       = check2_difference_h2_allcases.con_and_optional_LDSC_h2_allcases.con,
      check2.propSE                                                                             = check2.propSE,
      check3a_difference_rg_sub1.con_sub1.sub2_and_optional_LDSC_rg_sub1.con_sub1.sub2          = check3a_difference_rg_sub1.con_sub1.sub2_and_optional_LDSC_rg_sub1.con_sub1.sub2,
      check3b_difference_rg_sub2.con_sub1.sub2_and_optional_LDSC_rg_sub2.con_sub1.sub2          = check3b_difference_rg_sub2.con_sub1.sub2_and_optional_LDSC_rg_sub2.con_sub1.sub2,
      check3a.propSE                                                                            = check3a.propSE,
      check3b.propSE                                                                            = check3b.propSE,
      check4a_difference_rg_sub1.con_allcases.con_and_optional_LDSC_rg_sub1.con_allcases.con    = check4a_difference_rg_sub1.con_allcases.con_and_optional_LDSC_rg_sub1.con_allcases.con,
      check4a.propSE                                                                            = check4a.propSE,
      check4b_difference_rg_sub2.con_allcases.con_and_optional_LDSC_rg_sub2.con_allcases.con    = check4b_difference_rg_sub2.con_allcases.con_and_optional_LDSC_rg_sub2.con_allcases.con,
      check4b.propSE                                                                            = check4b.propSE,
      check5_difference_rg_sub1.sub2_ext_and_optional_LDSC_rg_sub1.sub2_ext                     = check5_difference_rg_sub1.sub2_ext_and_optional_LDSC_rg_sub1.sub2_ext,
      check5.propSE                                                                             = check5.propSE,
      check5a.internal.propSE                                                                   = check5a.internal.propSE,
      check5b.internal.propSE                                                                   = check5b.internal.propSE,
      check6_difference_rg_allcases.con_ext_and_optional_LDSC_rg_allcases.con_ext               = check6_difference_rg_allcases.con_ext_and_optional_LDSC_rg_allcases.con_ext,
      check6.propSE                                                                             = check6.propSE,

      optional_LDSC_rg_sub1.sub2_ext                                                            = optional_LDSC_rg_sub1.sub2_ext,
      optional_LDSC_rg_se_sub1.sub2_ext                                                         = optional_LDSC_rg_se_sub1.sub2_ext,
      optional_LDSC_rg_allcases.con_ext                                                         = optional_LDSC_rg_allcases.con_ext,
      optional_LDSC_rg_se_allcases.con_ext                                                      = optional_LDSC_rg_se_allcases.con_ext,
      rg_sub2.con_ext                                                                           = rg_sub2.con_ext,
      rg_sub1.con_ext                                                                           = rg_sub1.con_ext,
      rg_se_sub1.con_ext                                                                        = rg_se_sub1.con_ext,
      rg_se_sub2.con_ext                                                                        = rg_se_sub2.con_ext,

      a.deg_sub1.con_ext                                                                        = a.deg_sub1.con_ext,
      a.deg_sub2.con_ext                                                                        = a.deg_sub2.con_ext,
      a.deg_sub1.sub2_ext                                                                       = a.deg_sub1.sub2_ext,

      filter1                                                                                   = filter1,
      filter2                                                                                   = filter2,
      filter3                                                                                   = filter3,

      x.con                                                                                     = x.con,
      y.con                                                                                     = y.con,
      z.con                                                                                     = z.con,
      x.sub2                                                                                    = x.sub2,
      y.sub2                                                                                    = y.sub2,
      z.sub2                                                                                    = z.sub2,
      x.sub1                                                                                    = x.sub1,
      y.sub1                                                                                    = y.sub1,
      z.sub1                                                                                    = z.sub1,
      x.allcases                                                                                = x.allcases,
      y.allcases                                                                                = y.allcases,
      z.allcases                                                                                = z.allcases,
      x.popmean                                                                                 = x.popmean,
      y.popmean                                                                                 = y.popmean,
      z.popmean                                                                                 = z.popmean,
      x.ext                                                                                     = x.ext,
      y.ext                                                                                     = y.ext,
      z.ext                                                                                     = z.ext,
      x.extcontrol                                                                              = x.extcontrol,
      y.extcontrol                                                                              = y.extcontrol,
      z.extcontrol                                                                              = z.extcontrol,

      name_con                                                                                  = name_con,
      name_allcases                                                                             = name_allcases,
      name_sub1                                                                                 = name_sub1,
      name_sub2                                                                                 = name_sub2,
      name_ext                                                                                  = name_ext,
      N_sub1                                                                                    = N_sub1,
      N_sub2                                                                                    = N_sub2,
      filename                                                                                  = filename,
      filename3D                                                                                = filename3D,
      folder_location                                                                           = folder_location,
      plot_title                                                                                = plot_title,
      pop.prev_ext                                                                              = pop.prev_ext,

      # 3D
      h2_ext                                                                                    = h2_ext,
      h2_se_ext                                                                                 = h2_se_ext,
      rg_sub1.sub2_ext                                                                          = rg_sub1.sub2_ext,
      rg_allcases.con_ext                                                                       = rg_allcases.con_ext,
      internal_check_rg_sub1.con_ext                                                            = internal_check_rg_sub1.con_ext,
      internal_check_rg_sub2.con_ext                                                            = internal_check_rg_sub2.con_ext,
      internal_check_a_sub1.con_ext                                                             = internal_check_a_sub1.con_ext,
      internal_check_a_sub2.con_ext                                                             = internal_check_a_sub2.con_ext,
      a.deg_allcases.con_ext                                                                    = a.deg_allcases.con_ext,
      plot_3D                                                                                   = plot_3D,

      # Potential updated variables
      rg_sub1.con_sub2.con_original                                                             = rg_sub1.con_sub2.con_original,
      optional_LDSC_rg_allcases.con_sub1.con_original                                           = optional_LDSC_rg_allcases.con_sub1.con_original,
      optional_LDSC_rg_allcases.con_sub2.con_original                                           = optional_LDSC_rg_allcases.con_sub2.con_original,
      optional_LDSC_rg_sub1.con_sub1.sub2_original                                              = optional_LDSC_rg_sub1.con_sub1.sub2_original,
      optional_LDSC_rg_sub2.con_sub1.sub2_original                                              = optional_LDSC_rg_sub2.con_sub1.sub2_original,
      optional_LDSC_rg_sub1.sub2_ext_original                                                   = optional_LDSC_rg_sub1.sub2_ext_original,
      optional_LDSC_rg_allcases.con_ext_original                                                = optional_LDSC_rg_allcases.con_ext_original,
      rg_sub1.sub2_ext_original                                                                 = rg_sub1.sub2_ext_original,
      rg_allcases_con_ext_original                                                              = rg_allcases_con_ext_original,
      rg_sub1.con_ext_original                                                                  = rg_sub1.con_ext_original,
      rg_sub2.con_ext_original                                                                  = rg_sub2.con_ext_original,
      rg_allcases.con_ext_original                                                              = rg_allcases.con_ext_original)






#### Save, print to console and return --------------------------------------------------------------------------------------------------------



    # Save
    if(save_data == TRUE ) {

      if (plot_3D == TRUE & webversion == FALSE) {
        save_path <- file.path(folder_location, paste0(filename, "_", name_ext, ".3D.triangle_parameters.RData")) %>% gsub("^/", "", .)
        save(triangle.output.list, file = save_path)
        cli::cli_alert_info(paste0(" Data saved as ", folder_location,"/",filename,"_", name_ext,".3D.triangle_parameters.RData")) }
      if (plot_3D == TRUE & webversion == TRUE) {
        save_path <- file.path(folder_location, paste0(filename, "_", name_ext,".3D.triangle_parameters.RData")) %>% gsub("^/", "", .)
        save(triangle.output.list, file = save_path)
        log_fun("For more info see the logfile", type = "info") }
      if (plot_3D == FALSE & webversion == FALSE) {
        save_path <- file.path(folder_location, paste0(filename, ".2D.triangle_parameters.RData")) %>% gsub("^/", "", .)
        save(triangle.output.list, file = save_path)
        cli::cli_alert_info(paste0(" Data saved as ", folder_location,"/",filename,".2D.triangle_parameters.RData"))    }
      if (plot_3D == FALSE & webversion == TRUE) {
        save_path <- file.path(folder_location, paste0(filename, ".2D.triangle_parameters.RData")) %>% gsub("^/", "", .)
        save(triangle.output.list, file = save_path)
        log_fun("For more info see the logfile", type = "info") }

    } # End save_data = TRUE


#### Write log file  ---------------------------------------------------------------------------------------------------------


    # Get info that needs to be written

    # Change coordinates
    "coordinates_con" <- paste0("(", round(x.con, 3), ",", round(y.con, 3), ",",  round(z.con, 3), ")")
    "coordinates_sub1" <- paste0("(", round(x.sub1, 3), ",", round(y.sub1, 3), ",",  round(z.sub1, 3), ")")
    "coordinates_sub2" <- paste0("(", round(x.sub2, 3), ",", round(y.sub2, 3), ",",  round(z.sub2, 3), ")")
    "coordinates_allcases" <- paste0("(", round(x.allcases, 3), ",", round(y.allcases, 3), ",",  round(z.allcases, 3), ")")
    "coordinates_popmean" <- paste0("(", round(x.popmean, 3), ",", round(y.popmean, 3), ",",  round(z.popmean, 3), ")")
    "coordinates_ext" <- paste0("(", round(x.ext, 3), ",", round(y.ext, 3), ",",  round(z.ext, 3), ")")
    "coordinates_extcontrol" <- paste0("(", round(x.extcontrol, 3), ",", round(y.extcontrol, 3), ",",  round(z.extcontrol, 3), ")")

    # Required parameters for 2D
    if (webversion == FALSE) {
      to_logfile_required_input_pars_2D <- c(
        "h2_sub1.con", "h2_se_sub1.con",  "h2_sub2.con", "h2_se_sub2.con","rg_sub1.con_sub2.con", "rg_se_sub1.con_sub2.con",
        "name_sub1", "N_sub1", "name_sub2", "N_sub2","name_allcases", "name_con",  "plot_title","folder_location", "filename", "pop.prev_case") }
    if (webversion == TRUE) { # don't show folder_location
      to_logfile_required_input_pars_2D <- c(
        "h2_sub1.con", "h2_se_sub1.con",  "h2_sub2.con", "h2_se_sub2.con","rg_sub1.con_sub2.con", "rg_se_sub1.con_sub2.con",
        "name_sub1", "N_sub1", "name_sub2", "N_sub2","name_allcases", "name_con",  "plot_title", "filename", "pop.prev_case")  }

    # Optional parameters for 2D
    to_logfile_optional_input_pars_2D <- c(
      "optional_LDSC_h2_sub1.sub2", "optional_LDSC_h2_se_sub1.sub2","optional_LDSC_h2_allcases.con","optional_LDSC_h2_se_allcases.con",
      "optional_LDSC_rg_sub1.con_sub1.sub2", "optional_LDSC_rg_se_sub1.con_sub1.sub2","optional_LDSC_rg_sub2.con_sub1.sub2", "optional_LDSC_rg_se_sub2.con_sub1.sub2",
      "optional_LDSC_rg_allcases.con_sub1.con", "optional_LDSC_rg_se_allcases.con_sub1.con","optional_LDSC_rg_allcases.con_sub2.con","optional_LDSC_rg_se_allcases.con_sub2.con")

    # Output parameters for 2D
    to_logfile_output_pars_2D <- c(
      "h2_sub1.sub2","h2_allcases.con",
      "rg_sub1.con_sub1.sub2" ,"rg_sub2.con_sub1.sub2","rg_allcases.con_sub1.con","rg_allcases.con_sub2.con",
      "d_sub1.con","d_sub2.con","d_sub1.sub2" ,"d_allcases.con" ,"d_allcases.sub1","d_allcases.sub2",
      "a.rad_sub1.con_sub2.con", "a.deg_sub1.con_sub2.con",
      "a.rad_sub1.con_sub1.sub2","a.deg_sub1.con_sub1.sub2",
      "a.rad_sub2.con_sub1.sub2","a.deg_sub2.con_sub1.sub2",
      "a.rad_allcases.con_sub1.con", "a.deg_allcases.con_sub1.con",
      "a.rad_allcases.con_sub2.con", "a.deg_allcases.con_sub2.con",
      "coordinates_con", "coordinates_sub1", "coordinates_sub2", "coordinates_allcases", "coordinates_popmean")

    # Required parameters for 3D
    if (webversion == FALSE) {
      to_logfile_required_input_pars_3D <- c(
        "h2_sub1.con", "h2_se_sub1.con",  "h2_sub2.con", "h2_se_sub2.con","rg_sub1.con_sub2.con", "rg_se_sub1.con_sub2.con", # Same as 2D
        "name_sub1", "N_sub1", "name_sub2", "N_sub2","name_allcases", "name_con",  "plot_title","folder_location", "filename", "pop.prev_case", # Same as 2D
        "h2_ext","h2_se_ext", "rg_sub1.con_ext","rg_se_sub1.con_ext","rg_sub2.con_ext","rg_se_sub2.con_ext","name_ext","pop.prev_ext") } # Required 3D
    if (webversion == TRUE) { # Don't show folder location
      to_logfile_required_input_pars_3D <- c(
        "h2_sub1.con", "h2_se_sub1.con",  "h2_sub2.con", "h2_se_sub2.con","rg_sub1.con_sub2.con", "rg_se_sub1.con_sub2.con", # Same as 2D
        "name_sub1", "N_sub1", "name_sub2", "N_sub2","name_allcases", "name_con",  "plot_title","filename", "pop.prev_case", # Same as 2D
        "h2_ext","h2_se_ext", "rg_sub1.con_ext","rg_se_sub1.con_ext","rg_sub2.con_ext","rg_se_sub2.con_ext","name_ext","pop.prev_ext") } # Required 3D

    # Optional parameters for 3D
    to_logfile_optional_input_pars_3D <- c(
      "optional_LDSC_rg_sub1.sub2_ext", "optional_LDSC_rg_se_sub1.sub2_ext","optional_LDSC_rg_allcases.con_ext", "optional_LDSC_rg_se_allcases.con_ext")

    # Output parameters for 3D
    to_logfile_output_pars_3D <- c(
      "h2_sub1.sub2","h2_allcases.con", # Same as 2D
      "rg_sub1.con_sub1.sub2" ,"rg_sub2.con_sub1.sub2","rg_allcases.con_sub1.con","rg_allcases.con_sub2.con", # Same as 2D
      "rg_sub1.sub2_ext", "rg_allcases.con_ext", # 3D
      "d_sub1.con","d_sub2.con","d_sub1.sub2" ,"d_allcases.con" ,"d_allcases.sub1","d_allcases.sub2", # Same as 2D
      "a.rad_sub1.con_sub2.con", "a.deg_sub1.con_sub2.con", # Same as 2D
      "a.rad_sub1.con_sub1.sub2","a.deg_sub1.con_sub1.sub2", # Same as 2D
      "a.rad_sub2.con_sub1.sub2","a.deg_sub2.con_sub1.sub2", # Same as 2D
      "a.rad_allcases.con_sub1.con", "a.deg_allcases.con_sub1.con", # Same as 2D
      "a.rad_allcases.con_sub2.con", "a.deg_allcases.con_sub2.con", # Same as 2D
      "a.rad_sub1.sub2_ext", "a.deg_sub1.sub2_ext",# 3D
      "a.rad_allcases.con_ext","a.deg_allcases.con_ext",# 3D
      "coordinates_con", "coordinates_sub1", "coordinates_sub2", "coordinates_allcases", "coordinates_popmean", # Same as 2D
      "coordinates_ext", "coordinates_extcontrol") # 3D
    #"internal_check_rg_sub1.con_ext","internal_check_rg_sub2.con_ext","internal_check_a_sub1.con_ext", "internal_check_a_sub2.con_ext") # 3D

    # Get the time
    rounded_time <- as.POSIXct(round(as.numeric(Sys.time()) / 2) * 2, origin = "1970-01-01")

    # Mode in which GDVIS was run
    if (plot_2D.2D) {
      mode <- "2D   "
      to_logfile_required_input_pars <- to_logfile_required_input_pars_2D
      to_logfile_optional_input_pars <- to_logfile_optional_input_pars_2D
      to_logfile_output_pars <- to_logfile_output_pars_2D
    } else if (plot_3D) {
      mode <- "3D   "
      to_logfile_required_input_pars <- to_logfile_required_input_pars_3D
      to_logfile_optional_input_pars <- to_logfile_optional_input_pars_3D
      to_logfile_output_pars <- to_logfile_output_pars_3D
    } else {
      mode <- "2D   "
      to_logfile_required_input_pars <- to_logfile_required_input_pars_2D
      to_logfile_optional_input_pars <- to_logfile_optional_input_pars_2D
      to_logfile_output_pars <- to_logfile_output_pars_2D
    }

    # Open connection
    if (webversion == FALSE) {
      if (plot_3D == TRUE) {
        path_to_log <- file.path(folder_location,paste0(filename,"_", name_ext,".3D.log.txt")) %>% gsub("^/", "", .)
        con_log <- file(path_to_log, open = "w")
        cli::cli_alert_info(paste0(" For more info see ",folder_location,"/",filename,"_", name_ext,".3D.log.txt"))
      } else {
        path_to_log <- file.path(folder_location,paste0(filename,".2D.log.txt")) %>% gsub("^/", "", .)
        con_log <- file(path_to_log, open = "w")
        cli::cli_alert_info(paste0(" For more info see ",folder_location,"/",filename,".2D.log.txt"))
      }}

    if (webversion == TRUE) {
      if (plot_3D == TRUE) {
        con_path <- file.path(folder_location,paste0(filename,"_", name_ext,".3D.log.txt")) %>% gsub("^/", "", .)
        con_log <- file(con_path, open = "w")
        cli::cli_alert_info(paste0(" For more info see ",con_path))
      } else {
        con_path <- file.path(folder_location,paste0(filename,".2D.log.txt")) %>% gsub("^/", "", .)
        con_log <- file(con_path, open = "w")
        cli::cli_alert_info(paste0(" For more info see ",con_path))
      } }

    # Write general info
    writeLines("............................................................", con_log)
    writeLines(".   GDVIS software created by: A.B. Thijssen               .", con_log)
    writeLines(paste(sprintf(".   Log file created on: %s", rounded_time),"              ."), con_log)
    if (webversion == FALSE) {
      writeLines(paste(".   GDVIS mode used:", mode, "                                ."), con_log) }
    if (webversion == TRUE) {
      writeLines(paste(".   GDVIS mode used: online webversion                     ."), con_log)  }
    writeLines("............................................................", con_log)
    writeLines("    ", con_log)
    writeLines("    ", con_log)
    writeLines("    ", con_log)

    # Write message triangle not possible
    if (save_data == FALSE) {
      writeLines("Error:", con_log)
      writeLines(message_triangle_not_possible, con_log)
      writeLines("    ", con_log)
      writeLines("    ", con_log)
    }

    # Write required input data
    writeLines("Required input data", con_log)
    writeLines(strrep("-", 60), con_log)

    # Iterate over all input parameters and write to the logfile
    for (i in to_logfile_required_input_pars) {
      name <- i
      value <-  get(i)
      writeLines(sprintf("%-45s %s", name, value), con_log)  }

    # Write updated variables if exceeded [-1,1]
    if (length(log_messages) > 1) {
      writeLines("    ", con_log)
      writeLines(log_messages, con_log)
      writeLines("    ", con_log)
      writeLines("    ", con_log)
    } else {
      writeLines("    ", con_log)
      writeLines("    ", con_log)
    }

    # Write optional input data
    writeLines("Optional input data", con_log)
    writeLines(strrep("-", 60), con_log)

    # Iterate over all input parameters and write to the logfile
    for (i in to_logfile_optional_input_pars) {
      name <- i
      value <-  get(i)
      writeLines(sprintf("%-45s %s", name, value), con_log)  }
    writeLines("    ", con_log)
    writeLines("    ", con_log)

    # Write output data
    writeLines("Calculated output data", con_log)
    writeLines(strrep("-", 60), con_log)

    for (i in to_logfile_output_pars) {
      name <- i
      value <-  get(i)
      if (grepl("coordinates", name)) {
        rounded_value <- value    } else {
          rounded_value <- round(as.numeric(value), 3)      }
      writeLines(sprintf("%-45s %s", name, rounded_value), con_log)  }
    writeLines("    ", con_log)
    writeLines("    ", con_log)

    # Write filters
    writeLines("Check whether subgroups are genetically different", con_log)
    writeLines(strrep("-", 60), con_log)
    writeLines("When there is no suggested evidence of significant difference between subgroups, the triangle will not be plotted. Specifically, when rg is not significantly different from 1 (filter 1a) AND h2.sub1 is not clearly different from h2.sub2 (filter 1b and 1c), the triangle will not be plotted. In some occasions, a triangle will not be plotted while it would have been plotted with larger GWAS sample sizes (and therefore smaller standard errors). ", con_log)
    writeLines("    ", con_log)

    writeLines("Filter 1a:", con_log)
    if (outcome_filter_1a > 1) {
      writeLines(paste0("The genetic correlation between ", name_sub1, "_", name_con, "' and '", name_sub2, "_", name_con, "' is not significantly different from 1"), con_log) }
    if (outcome_filter_1a < 1) {
      writeLines(paste0("The genetic correlation between ", name_sub1, "_", name_con, "' and '", name_sub2, "_", name_con, "' is significantly different from 1"), con_log) }
    writeLines(paste0("The upper bound is calculated as ", rg_sub1.con_sub2.con, " (the genetic correlation) + 1.96 * ", rg_se_sub1.con_sub2.con, " (se) = ", outcome_filter_1a), con_log)
    writeLines("    ", con_log)
    writeLines("Filter 1b:", con_log)
    if (outcome_filter_1b == FALSE) {
      writeLines(paste0("The heritability of '", name_sub1, "_", name_con, "' does not lie within the confidence intervals of the heritability of '",  name_sub2, "_", name_con, "'"), con_log)    }
    if (outcome_filter_1b == TRUE) {
      writeLines(paste0("The heritability of '", name_sub1, "_", name_con, "' lies within the confidence intervals of the heritability of '",  name_sub2, "_", name_con, "'"), con_log)   }
    writeLines(paste0("The heritability of ", name_sub1,"_", name_con," is ", h2_sub1.con, ", the interval of ", name_sub2, "_", name_con, " is ", interval_filter_1b), con_log)
    writeLines("    ", con_log)
    writeLines("Filter 1c:", con_log)
    if (outcome_filter_1c == FALSE) {
      writeLines(paste0("The heritability of '", name_sub2, "_", name_con, "' does not lie within the confidence intervals of the heritability of '",  name_sub1, "_", name_con, "'"), con_log)    }
    if (outcome_filter_1c == TRUE) {
      writeLines(paste0("The heritability of '", name_sub2, "_", name_con, "' lies within the confidence intervals of the heritability of '",  name_sub1, "_", name_con, "'"), con_log)   }
    writeLines(paste0("The heritability of '", name_sub2, "_",name_con, "' is ", h2_sub2.con, ", the interval is ", interval_filter_1c), con_log)
    writeLines("    ", con_log)
    if (filter1 == 1) {
      writeLines("--> Filter 1 is not passed: The triangle will not be plotted.", con_log)   }
    if (filter1 == 0) {
      writeLines("--> Filter 1 is passed: subgroups are genetically different from each other", con_log) }
    writeLines("    ", con_log)
    writeLines("Filter 2:", con_log)
    if(filter2 == 1) {
      writeLines(paste0("The heritability of '", name_sub1, "_", name_con, "' is not significantly different from 0"), con_log)
      writeLines(paste0("The heritability is ",h2_sub1.con, " with a lower bound of the confidence interval of ",round(outcome_filter_2,3)), con_log)
      writeLines("--> Filter 2 is not passed: The triangle will not be plotted.", con_log) }
    if(filter2 == 0) {
      writeLines(paste0("The heritability of '", name_sub1, "_", name_con, "' is significantly different from 0"), con_log)
      writeLines(paste0("The heritability is ",h2_sub1.con, " with a lower bound of the confidence interval of ",round(outcome_filter_2,3)), con_log)
      writeLines("--> Filter 2 is passed", con_log) }
    writeLines("    ", con_log)
    writeLines("Filter 3:", con_log)
    if(filter3 == 1) {
      writeLines(paste0("The heritability of '", name_sub2, "_", name_con, "' is not significantly different from 0"), con_log)
      writeLines(paste0("The heritabiltiy is ", h2_sub2.con, " where the lower bound of the confidence interval is ",round(outcome_filter_3,3)), con_log)
      writeLInes("--> Filter 3 is not passed: The triangle will not be plotted.", con_log) }
    if(filter3 == 0) {
      writeLines(paste0("The heritability of '", name_sub2, "_", name_con, "' is significantly different from 0"), con_log)
      writeLines(paste0("The heritabiltiy is ", h2_sub2.con, " where the lower bound of the confidence interval is ",round(outcome_filter_3,3)), con_log)
      writeLines("--> Filter 3 is passed", con_log) }
    writeLines("    ", con_log)
    writeLines("    ", con_log)


    # Write double checks
    writeLines("Optional double checks", con_log)
    writeLines(strrep("-", 60), con_log)


    writeLines("Double check 1:", con_log)
    writeLines(paste0("GDVIS calculated the heritability of ", name_sub1, " versus ", name_sub2, " as ", round(h2_sub1.sub2, 4)), con_log)
    # Don't write this stuff if no optional input was provided
    if (is.na(check1_passed_h2_sub1.sub2_and_optional_LDSC_h2_sub1.sub2)) {writeLines("No optional input for this check was provided. ", con_log)}  else {
      # What to write as dubbelcheck
      writeLines(paste0("The LDSC estimated value provided as input is ", optional_LDSC_h2_sub1.sub2, " (se = ", optional_LDSC_h2_se_sub1.sub2, ")"), con_log)
      if (check1_passed_h2_sub1.sub2_and_optional_LDSC_h2_sub1.sub2 == TRUE) { writeLines(paste0("The GDVIS value lies within the 95% CI ",optional_LDSC_h2_sub1.sub2, " +- 1.96 * ", optional_LDSC_h2_se_sub1.sub2,"."), con_log)}
      if (check1_passed_h2_sub1.sub2_and_optional_LDSC_h2_sub1.sub2 == FALSE) {  writeLines(paste0("The GDVIS value does not lie within the 95% CI +- 1.96 * ", optional_LDSC_h2_se_sub1.sub2,"."), con_log)}
      writeLines(paste0("The difference between ",round(h2_sub1.sub2, 4), " and ", round(optional_LDSC_h2_sub1.sub2,4), " is ",  round(check1_difference_h2_sub1.sub2_and_optional_LDSC_h2_sub1.sub2, 4), ", this is ", round(check1.propSE, 4), " * se.") , con_log) }
    writeLines("    ", con_log)


    writeLines("Double check 2:", con_log)
    writeLines(paste0("GDVIS calculated the heritability of allcases versus ",name_con," as ", round(h2_allcases.con,4) ), con_log)
    # Don't write this stuff if no optional input was provided
    if (is.na(check2_passed_h2_allcases.con_and_optional_LDSC_h2_allcases.con)) {writeLines("No optional input for this check was provided. ", con_log)}  else {
      # What to write as dubbelcheck
      writeLines(paste0("The LDSC estimated value provided as input is  ", optional_LDSC_h2_allcases.con, " (se = ", optional_LDSC_h2_se_allcases.con, ")"), con_log)
      if (check2_passed_h2_allcases.con_and_optional_LDSC_h2_allcases.con == TRUE)  { writeLines(paste0("The GDVIS value lies within the 95%  CI +- 1.96 * ", optional_LDSC_h2_se_allcases.con,"."), con_log)}
      if (check2_passed_h2_allcases.con_and_optional_LDSC_h2_allcases.con == FALSE) {  writeLines(paste0("The GDVIS value does not lie within the 95% CI +- 1.96 * ", optional_LDSC_h2_se_allcases.con,"."), con_log)}
      writeLines(paste0("The difference between ",round(h2_allcases.con, 4), " and ", round(optional_LDSC_h2_allcases.con, 4), " is ",  round(check2_difference_h2_allcases.con_and_optional_LDSC_h2_allcases.con,4), ", this is ", round(check2.propSE, 4), " * se.") , con_log) }
    writeLines("    ", con_log)

    writeLines("Double check 3a:", con_log)
    writeLines(paste0("GDVIS calculated the genetic correlation between ", name_sub1, ".vs." , name_con, " and ", name_sub1, ".vs.", name_sub2, " as ", round(rg_sub1.con_sub1.sub2,4)), con_log)
    # Don't write this stuff if no optional input was provided
    if (is.na(check3a_passed_rg_sub1.con_sub1.sub2_and_optional_LDSC_rg_sub1.con_sub1.sub2)) {writeLines("No optional input for this check was provided. ", con_log)}  else {
      # What to write as dubbelcheck
      writeLines(paste0("The LDSC estimated value provided as input is ", optional_LDSC_rg_sub1.con_sub1.sub2, " (se = ", optional_LDSC_rg_se_sub1.con_sub1.sub2, ")"), con_log)
      if (check3a_passed_rg_sub1.con_sub1.sub2_and_optional_LDSC_rg_sub1.con_sub1.sub2 == TRUE)  { writeLines(paste0("The GDVIS value lies within the 95% CI +- 1.96 * ", round(optional_LDSC_rg_se_sub1.con_sub1.sub2, 4),"."), con_log)}
      if (check3a_passed_rg_sub1.con_sub1.sub2_and_optional_LDSC_rg_sub1.con_sub1.sub2 == FALSE) { writeLines(paste0("The GDVIS value does not lie within the 95% CI +- 1.96 * ", round(optional_LDSC_rg_se_sub1.con_sub1.sub2, 4),"."), con_log)}
      writeLines(paste0("The difference between ",round(rg_sub1.con_sub1.sub2, 4), " and ", round(optional_LDSC_rg_sub1.con_sub1.sub2, 4), " is ",  round(check3a_difference_rg_sub1.con_sub1.sub2_and_optional_LDSC_rg_sub1.con_sub1.sub2, 4), ", this is ", round(check3a.propSE, 4), " * se.") , con_log) }
    writeLines("    ", con_log)

    writeLines("Double check 3b:", con_log)
    writeLines(paste0("GDVIS calculated the genetic correlation between ", name_sub2, ".vs." , name_con, " and ", name_sub1, ".vs.", name_sub2, " as ", round(rg_sub2.con_sub1.sub2, 4)), con_log)
    # Don't write this stuff if no optional input was provided
    if (is.na(check3b_passed_rg_sub2.con_sub1.sub2_and_optional_LDSC_rg_sub2.con_sub1.sub2)) {writeLines("No optional input for this check was provided. ", con_log)}  else {
      # What to write as dubbelcheck
      writeLines(paste0("The LDSC estimated value provided as input is ", optional_LDSC_rg_sub2.con_sub1.sub2, " (se = ", optional_LDSC_rg_se_sub2.con_sub1.sub2, ")"), con_log)
      if (check3b_passed_rg_sub2.con_sub1.sub2_and_optional_LDSC_rg_sub2.con_sub1.sub2 == TRUE)  { writeLines(paste0("The GDVIS value lies within the 95%  CI +- 1.96 * ", round(optional_LDSC_rg_se_sub1.con_sub1.sub2, 4),"."), con_log)}
      if (check3b_passed_rg_sub2.con_sub1.sub2_and_optional_LDSC_rg_sub2.con_sub1.sub2 == FALSE) { writeLines(paste0("The GDVIS value does not lie within the 95% CI +- 1.96 * ", round(optional_LDSC_rg_se_sub2.con_sub1.sub2, 4),"."), con_log)}
      writeLines(paste0("The difference between ",round(rg_sub2.con_sub1.sub2, 4), " and ", round(optional_LDSC_rg_sub2.con_sub1.sub2, 4), " is ",  round(check3b_difference_rg_sub2.con_sub1.sub2_and_optional_LDSC_rg_sub2.con_sub1.sub2, 4), ", this is ", round(check3b.propSE, 4), " * se.") , con_log)  }
    writeLines("    ", con_log)


    writeLines("Double check 4a:", con_log)
    writeLines(paste0("GDVIS calculated the genetic correlation between allcases.vs." , name_con, " and ", name_sub1, ".vs.", name_con, " as ", round(rg_allcases.con_sub1.con,4)), con_log)
    # Don't write this stuff if no optional input was provided
    if (is.na(check4a_passed_rg_sub1.con_allcases.con_and_optional_LDSC_rg_sub1.con_allcases.con)) {writeLines("No optional input for this check was provided. ", con_log)}  else {
      # What to write as dubbelcheck
      writeLines(paste0("The LDSC estimated value provided as input is ", optional_LDSC_rg_allcases.con_sub1.con, " (se = ", optional_LDSC_rg_se_allcases.con_sub1.con, ")"), con_log)
      if (check4a_passed_rg_sub1.con_allcases.con_and_optional_LDSC_rg_sub1.con_allcases.con == TRUE)  { writeLines(paste0("The GDVIS value lies within the 95% CI +- 1.96 * ", round(optional_LDSC_rg_se_allcases.con_sub1.con, 4),"."), con_log)}
      if (check4a_passed_rg_sub1.con_allcases.con_and_optional_LDSC_rg_sub1.con_allcases.con == FALSE) { writeLines(paste0("The GDVIS value does not lie within the 95% CI +- 1.96 * ", round(optional_LDSC_rg_se_allcases.con_sub1.con, 4),"."), con_log)}
      writeLines(paste0("The difference between ",round(rg_allcases.con_sub1.con, 4), " and ", round(optional_LDSC_rg_allcases.con_sub1.con, 4), " is ",  round(check4a_difference_rg_sub1.con_allcases.con_and_optional_LDSC_rg_sub1.con_allcases.con, 4), ", this is ", round(check1.propSE, 4), " * se.") , con_log) }
    writeLines("    ", con_log)

    writeLines("Double check 4b:", con_log)
    writeLines(paste0("GDVIS calculated the genetic correlation between allcases.vs." , name_con, " and ", name_sub2, ".vs.", name_con, " is calculated to be ", round(rg_allcases.con_sub2.con, 4)), con_log)
    # Don't write this stuff if no optional input was provided
    if (is.na(check4b_passed_rg_sub2.con_allcases.con_and_optional_LDSC_rg_sub2.con_allcases.con)) {writeLines("No optional input for this check was provided. ", con_log)}  else {
      # What to write as dubbelcheck
      writeLines(paste0("The LDSC estimated value provided as input is ", optional_LDSC_rg_allcases.con_sub2.con, " (se = ", optional_LDSC_rg_se_allcases.con_sub2.con, ")"), con_log)
      if (check4b_passed_rg_sub2.con_allcases.con_and_optional_LDSC_rg_sub2.con_allcases.con == TRUE)  { writeLines(paste0("The GDVIS value lies within the CI +- 1.96 * ", round(optional_LDSC_rg_se_allcases.con_sub2.con, 4),"."), con_log)}
      if (check4b_passed_rg_sub2.con_allcases.con_and_optional_LDSC_rg_sub2.con_allcases.con == FALSE) { writeLines(paste0("The GDVIS value does not lie within the CI +- 1.96 * ", round(optional_LDSC_rg_se_allcases.con_sub2.con, 4),"."), con_log)}
      writeLines(paste0("The difference between ",round(rg_allcases.con_sub2.con, 4), " and ", round(optional_LDSC_rg_allcases.con_sub2.con, 4), " is ",  round(check4b_difference_rg_sub2.con_allcases.con_and_optional_LDSC_rg_sub2.con_allcases.con, 4), ", this is ", round(check1.propSE, 4), " * se.") , con_log)  }
    writeLines("    ", con_log)

    # 3D checks
    if (plot_3D == TRUE & save_data == TRUE) {

      writeLines("Internal 3D double checks", con_log)
      writeLines(strrep("-", 60), con_log)
      writeLines("Double check 5:", con_log)
      writeLines(paste0("GDVIS calculated the genetic correlation between ", name_sub1,".vs.",name_sub2, " and ", name_ext, " as ", round(rg_sub1.sub2_ext, 4)), con_log)

      if (is.na(optional_LDSC_rg_sub1.sub2_ext)) {writeLines("No optional input for this check was provided. ", con_log)}  else {
        # What to write as dubbelcheck
        writeLines(paste0("The LDSC estimated value provided as input is ", optional_LDSC_rg_sub1.sub2_ext, " (se = ", optional_LDSC_rg_se_sub1.sub2_ext, ")"), con_log)
        if (check5_passed_rg_sub1.sub2_ext_and_optional_LDSC_rg_sub1.sub2_ext == TRUE)  { writeLines(paste0("The GDVIS value lies within the CI +- 1.96 * ", round(optional_LDSC_rg_se_sub1.sub2_ext, 4),"."), con_log)}
        if (check5_passed_rg_sub1.sub2_ext_and_optional_LDSC_rg_sub1.sub2_ext == FALSE) { writeLines(paste0("The GDVIS value does not lie within the CI +- 1.96 * ", round(optional_LDSC_rg_se_sub1.sub2_ext, 4),"."), con_log)}
        writeLines(paste0("The difference between ", round(rg_sub1.sub2_ext, 4), " and ", optional_LDSC_rg_sub1.sub2_ext, " is ",  round(check5_difference_rg_sub1.sub2_ext_and_optional_LDSC_rg_sub1.sub2_ext, 4), ", this is ", round(check5.propSE, 4), " * se.") , con_log)  }
      writeLines("    ", con_log)

      writeLines("Double check 6:", con_log)
      writeLines(paste0("GDVIS calculated the genetic correlation between ", name_allcases,".vs.", name_con, " and ", name_ext, " as ", round(a.rad_allcases.con_ext,4) ), con_log)

      if (is.na(optional_LDSC_rg_allcases.con_ext)) {writeLines("No optional input for this check was provided. ", con_log)}  else {
        # What to write as dubbelcheck
        writeLines(paste0("The LDSC estimated value provided as input is ", optional_LDSC_rg_allcases.con_ext, " (se = ", optional_LDSC_rg_se_allcases.con_ext, ")"), con_log)
        if (check6_passed_rg_allcases.con_ext_and_optional_LDSC_rg_allcases.con_ext == TRUE)  { writeLines(paste0("The GDVIS value lies within the CI +- 1.96 * ", round(optional_LDSC_rg_se_allcases.con_ext, 4),"."), con_log)}
        if (check6_passed_rg_allcases.con_ext_and_optional_LDSC_rg_allcases.con_ext == FALSE) { writeLines(paste0("The GDVIS value does not lie within the CI +- 1.96 * ", round(optional_LDSC_rg_se_allcases.con_ext, 4),"."), con_log)}
        writeLines(paste0("The difference between ", round(rg_allcases.con_ext, 4), " and ", optional_LDSC_rg_allcases.con_ext, " is ",  round(check6_difference_rg_allcases.con_ext_and_optional_LDSC_rg_allcases.con_ext, 4), ", this is ", round(check6.propSE, 4), " * se.") , con_log)  }
      writeLines("    ", con_log)

      writeLines("Internal check rg_sub1.con_ext", con_log)
      writeLines(paste0("After finding 3D coordinates, GDVIS re-estimated the genetic correlation between ", name_sub1,".vs.", name_con, " and ", name_ext, " as ", round(internal_check_rg_sub1.con_ext,4) ), con_log)
      writeLines(paste0("The original LDSC estimated value is ", rg_sub1.con_ext, " (se = ", rg_se_sub1.con_ext, ")"), con_log)
      writeLines(paste0("The difference between ", round(internal_check_rg_sub1.con_ext, 4), " and ", rg_sub1.con_ext, " is ",  round(abs((internal_check_rg_sub1.con_ext  -  rg_sub1.con_ext)), 4)) , con_log)
      writeLines("    ", con_log)

      writeLines("Internal check rg_sub2.con_ext", con_log)
      writeLines(paste0("After finding 3D coordinates, GDVIS re-estimated the genetic correlation between ", name_sub2,".vs.", name_con, " and ", name_ext, " as ", round(internal_check_rg_sub2.con_ext,4) ), con_log)
      writeLines(paste0("The original LDSC estimated value is ", rg_sub2.con_ext, " (se = ", rg_se_sub2.con_ext, ")"), con_log)
      writeLines(paste0("The difference between ", round(internal_check_rg_sub2.con_ext, 4), " and ", rg_sub2.con_ext, " is ",  round(abs((internal_check_rg_sub2.con_ext  -  rg_sub2.con_ext)), 4)) , con_log)
      writeLines("    ", con_log)

    }

    if (plot_3D == TRUE & save_data == FALSE) {
      writeLines("Internal 3D double checks", con_log)
      writeLines(strrep("-", 60), con_log)
      writeLines("No internal 3D checks done", con_log)
    }

    # Close the connection
    close(con_log)

    # Message if 2D.2D
    if (plot_2D.2D == T) {
      cli::cli_alert_success("GDVIS calc succesfully finished 2D calculations")
    }

    return(mget(ls()))



    } # END CODE CHUNK


#### Run code chunk ----------------------------------------------------------------------------------------------------------


    # Run this for 2D.2D
    if (plot_2D.2D == T) {
      for (triangle in c("triangle1", "triangle2")) {   returned_variables <- run_code_chunk(triangle)  }}

    # Run this for normal 2D and for 3D
    if (plot_2D.2D == F & plot_CD == F) {
      if (webversion == TRUE) {
        returned_variables <- run_code_chunk()
        list2env(returned_variables, envir = temp_triangle_env) # Return the file path of the saved RData object
        if(save_data == TRUE) {  log_fun("GDVIS calc successfully finished!", type = "success") }
        return(list(
          save_path = save_path,
          log_path = con_path    ))      }

      if (webversion == FALSE) {
        returned_variables <- run_code_chunk() }
      list2env(returned_variables, envir = temp_triangle_env) # Return the file path of the saved RData object
      if(save_data == TRUE) {
        cli::cli_alert_success("GDVIS calc succesfully finished!")
        return(save_path) }   }

   # list2env(returned_variables, envir = temp_triangle_env)


#### Calc 2D-2D --------------------------------------------------------------------------------------------------------------



    if (plot_2D.2D == T) {

    cli::cli_inform(paste0(" Starting 2D.2D calculations" ))

    # Load calculated parameters
    if (triangle1.folder_location == "") {    parameters_triangle1 <- list.files(pattern = paste0(triangle1.plot_title, ".2D.triangle_parameters"), full.names = T) }
    if (triangle2.folder_location == "") {    parameters_triangle2 <- list.files(pattern = paste0(triangle2.plot_title, ".2D.triangle_parameters"), full.names = T) }
    if (triangle1.folder_location != "") {    parameters_triangle1 <- list.files(path = triangle1.folder, pattern = paste0(triangle1.plot_title, ".2D.triangle_parameters"), full.names = T) }
    if (triangle2.folder_location != "") {    parameters_triangle2 <- list.files(path = triangle2.folder, pattern = paste0(triangle2.plot_title, ".2D.triangle_parameters"), full.names = T) }

    # Load the output list and immediately rename it
    load(parameters_triangle1, envir = temp_triangle_env)
    triangle1_list <- triangle.output.list
    load(parameters_triangle2, envir = temp_triangle_env)
    triangle2_list <- triangle.output.list

    # Rename the elements with a prefix
    names(triangle1_list) <- paste0("triangle1.", names(triangle1_list))
    names(triangle2_list) <- paste0("triangle2.", names(triangle2_list))

    # Make lists available
    list2env(triangle1_list, envir = temp_triangle_env)
    list2env(triangle2_list, envir = temp_triangle_env)

    # Check that the allcases h2's are not too different
    if( abs(as.numeric(triangle1.h2_allcases.con) - as.numeric(triangle2.h2_allcases.con)) > 0.01) {
      different_h2 <- TRUE
      cli::cli_alert_danger("Triangles have different allcases h2's and do thus not compare similar cases")
      invokeRestart("abort")  # Stops execution without an error message
    }

    # Put coordinates in vectors
    triangle1.coord_sub1 <- c(triangle1.x.sub1,triangle1.y.sub1,triangle1.z.sub1)
    triangle1.coord_sub2 <- c(triangle1.x.sub2,triangle1.y.sub2,triangle1.z.sub2)
    triangle2.coord_sub1 <- c(triangle2.x.sub1,triangle2.y.sub1,triangle2.z.sub1)
    triangle2.coord_sub2 <- c(triangle2.x.sub2,triangle2.y.sub2,triangle2.z.sub2)
    coord_allcases <- c(triangle2.x.allcases,triangle2.y.allcases,triangle2.z.allcases)
    coord_con <- c(triangle2.x.con,triangle2.y.con, triangle2.z.con)

    # Save the lengths before rotation
    length_subsub_2r <- f.length_from_coords(triangle2.coord_sub1, triangle2.coord_sub2)
    length_sub1_con_2r <- f.length_from_coords(triangle2.coord_sub1, coord_con)
    length_sub2_con_2r <- f.length_from_coords(triangle2.coord_sub2, coord_con)


    # Step 2: Check if values exceed [-1,1] and restrain if necessary

    # Get the variables that need to be checked
    rg_variables <- c("rg_triangle1.sub1.con_triangle2.sub1.con", "rg_triangle1.sub2.con_triangle2.sub1.con", "rg_triangle1.sub1.con_triangle2.sub2.con", "rg_triangle1.sub2.con_triangle2.sub2.con" )

    # Initialize log variable
    log_messages <- c()

    # Loop through the variables and clip to [-1,1] if needed
    for (var in rg_variables) {
      original <- get(var)  # Retrieve the original values
      if (is.na(original)) {next } # Skip to the next variable in the loop
      updated <- pmin(pmax(original, -1), 1)  # Perform the update if needed

      # Compare and log a message if values are updated
      if (!(original == updated)) {
        log_messages <- c(log_messages, sprintf(paste(" Values in %s were updated from", get(var), "to", updated, "to perform calculations"), var))
        assign(var, as.numeric(updated))  # Save the updated values back to the variable
        assign(paste0(var, "_original"), original) # Keep track of the original variable
      }}

    # Print log messages to the console
    if (length(log_messages)>0) { cli::cli_alert_info(paste(log_messages, collapse = "\n")) }

    # Step 3: Test if triangles are possible
    f.triangle_possible_from_rg <- function(rg1, rg2, rg3) {
      length1 <- f.length_from_angle_and_lines(1, 1, rg2angle(rg1))
      length2 <- f.length_from_angle_and_lines(1, 1, rg2angle(rg2))
      length3 <- f.length_from_angle_and_lines(1, 1, rg2angle(rg3))
      test <- f.triangle_possible(length1, length2, length3)
      return(test)    }

    test1 <- f.triangle_possible_from_rg(rg_triangle1.sub2.con_triangle2.sub1.con, rg_triangle1.sub1.con_triangle2.sub1.con, triangle1.rg_sub1.con_sub2.con)
    test2 <- f.triangle_possible_from_rg(rg_triangle1.sub1.con_triangle2.sub2.con, rg_triangle1.sub2.con_triangle2.sub2.con, triangle1.rg_sub1.con_sub2.con)
    test3 <- f.triangle_possible_from_rg(rg_triangle1.sub1.con_triangle2.sub1.con, rg_triangle1.sub1.con_triangle2.sub2.con, triangle2.rg_sub1.con_sub2.con)
    test4 <- f.triangle_possible_from_rg(rg_triangle1.sub2.con_triangle2.sub1.con, rg_triangle1.sub2.con_triangle2.sub2.con, triangle2.rg_sub1.con_sub2.con)


    # Update save_data with the function result
    #save_data <- f.triangle_possible(line1, line2, line3)
    if (any(!c(test1, test2, test3, test4))) {
      message_triangle_not_possible_2D.2D <- "GDVIS determined that the triangle is not possible, possibly due to inconsistency in input parameters. Please double-check your input. If the issue persists, it may be due to inconsistency in LDSC estimates, for instance because the rg between sub1.con and sub2.con is close to 1"
      cli::cli_alert_danger(cli::col_red(message_triangle_not_possible_2D.2D))
      save_data <- FALSE
      return(save_data)
    }
    if (any(!c(test1, test2, test3, test4)) & webversion == TRUE) {
      log_fun(message_triangle_not_possible_2D.2D, type = "info")
    }


    # Step 4: Rotate the second triangle
    if (save_data) {

    # List the wcontrols se objects
    se_list <- list(
      rg_se_triangle1.sub1.con_triangle2.sub1.con = rg_se_triangle1.sub1.con_triangle2.sub1.con,
      rg_se_triangle1.sub2.con_triangle2.sub1.con = rg_se_triangle1.sub2.con_triangle2.sub1.con,
      rg_se_triangle1.sub1.con_triangle2.sub2.con = rg_se_triangle1.sub1.con_triangle2.sub2.con,
      rg_se_triangle1.sub2.con_triangle2.sub2.con = rg_se_triangle1.sub2.con_triangle2.sub2.con  )

    # Find the SE label with the smallest value and the corresponding object
    smallest_se_name <- names(se_list)[which.min(unlist(se_list))]
    smallest_se_value <- se_list[[smallest_se_name]]
    rg_name <- sub("_se", "", smallest_se_name)
    rg_smallest_se_value <- get(rg_name)

    # Get the corresponding subtypes
    flat_subtype <- stringr::str_extract(smallest_se_name, "(?<=triangle1\\.)[^.]+")
    rotated_subtype <- stringr::str_extract(smallest_se_name, "(?<=triangle2\\.)[^.]+")

    # Make vectors for wcontrols lines
    coord_2r <- paste0("triangle2.coord_", rotated_subtype)
    vector_2r <- get(coord_2r) - coord_con

    # Vector to align to
    coord_flat_subtype <- paste0("triangle1.coord_", flat_subtype)
    vector_flat <- get(coord_flat_subtype) - coord_con

    # Determine how far it should rotate
    rotation_angle <- rg2angle(rg_smallest_se_value)

    # Function to rotate the triangle
    adjust_angle_to_target <- function(vector1, vector2, target_angle, max_iterations = 10000000, tolerance = 1e-8) {
      current_angle <- acos(sum(vector1 * vector2) / (sqrt(sum(vector1^2)) * sqrt(sum(vector2^2))))
      iteration <- 0

      while (abs(current_angle - target_angle) > tolerance && iteration < max_iterations) {
        delta_theta <- target_angle - current_angle

        # Rotation matrix around the Y-axis
        rotation_matrix_Y <- matrix(c(
          cos(delta_theta), 0, sin(delta_theta),  0, 1, 0,
          -sin(delta_theta), 0, cos(delta_theta)), nrow = 3, byrow = TRUE)

        # Rotate the vector and recalculate the angle
        vector2 <- as.vector(rotation_matrix_Y %*% vector2)
        current_angle <- acos(sum(vector1 * vector2) / (sqrt(sum(vector1^2)) * sqrt(sum(vector2^2))))

        iteration <- iteration + 1  }

      return(list(rotated_vector = vector2, final_angle = current_angle))    }


    # Adjust the angle iteratively
    results <- adjust_angle_to_target(vector_flat, vector_2r, rotation_angle)

    # Extract the final rotated vector and angle
    rotated_vector <- results$rotated_vector
    final_angle <- results$final_angle
    angle_difference <- rotation_angle - final_angle

    # Add rotated vector to the controls point
    coord_rotated <- coord_con + rotated_vector

    # Get dir vector from sub_2r to allcases
    dir.vect_rsub_allcases <- coord_allcases - coord_rotated

    # Normalize the direction vector
    dir.vect_rsub_allcases <- dir.vect_rsub_allcases / sqrt(sum(dir.vect_rsub_allcases^2))

    # Update the coordinates for triangle 2
    if (rotated_subtype == "sub1") {
      triangle2.x.sub1 <- coord_rotated[1]
      triangle2.y.sub1 <- coord_rotated[2]
      triangle2.z.sub1 <- coord_rotated[3]

      # Update coord
      triangle2.coord_sub1 <- coord_rotated

      # Calculate the new end coordinates for other subtype
      triangle2.coord_sub2 <- coord_rotated + (dir.vect_rsub_allcases * length_subsub_2r)

      triangle2.x.sub2 <- triangle2.coord_sub2[1]
      triangle2.y.sub2 <- triangle2.coord_sub2[2]
      triangle2.z.sub2 <- triangle2.coord_sub2[3]

    } else if (rotated_subtype == "sub2") {
      triangle2.x.sub2 <- coord_rotated[1]
      triangle2.y.sub2 <- coord_rotated[2]
      triangle2.z.sub2 <- coord_rotated[3]

      # Update coord
      triangle2.coord_sub2 <- coord_rotated

      # Calculate the new end coordinates for other subtype
      triangle2.coord_sub1 <- coord_rotated + (dir.vect_rsub_allcases * length_subsub_2r)

      triangle2.x.sub1 <- triangle2.coord_sub1[1]
      triangle2.y.sub1 <- triangle2.coord_sub1[2]
      triangle2.z.sub1 <- triangle2.coord_sub1[3]
    }

    # New rg's
    rg_rsub1.rsub2_sub1.sub2 <- angle2rg(f.angle_from_coords(triangle1.coord_sub2, triangle1.coord_sub1, triangle2.coord_sub2, triangle2.coord_sub1))
    rg_rsub1.con_sub1.con    <- angle2rg(f.angle_from_coords(coord_con, triangle1.coord_sub1, coord_con, triangle2.coord_sub1))
    rg_rsub1.con_sub2.con    <- angle2rg(f.angle_from_coords(coord_con, triangle1.coord_sub2, coord_con, triangle2.coord_sub1))
    rg_rsub2.con_sub1.con    <- angle2rg(f.angle_from_coords(coord_con, triangle1.coord_sub1, coord_con, triangle2.coord_sub2))
    rg_rsub2.con_sub2.con    <- angle2rg(f.angle_from_coords(coord_con, triangle1.coord_sub2, coord_con, triangle2.coord_sub2))

    # Length checks
    if( abs(f.length_from_coords(triangle2.coord_sub1, triangle2.coord_sub2) - length_subsub_2r) > 0.01) { stop(cli::cli_alert_danger("Length check rotated sub1_sub2 failed"))} # Check sub1_sub2 length still the same
    if( abs(f.length_from_coords(coord_con, triangle2.coord_sub1) - length_sub1_con_2r)  > 0.01) { stop(cli::cli_alert_danger("Length check rotated sub1_con failed"))} # Check sub1_con length still the same
    if( abs(f.length_from_coords(coord_con, triangle2.coord_sub2) - length_sub2_con_2r) > 0.01) { stop(cli::cli_alert_danger("Length check rotated sub2_con failed"))} # Check sub2_con length still the same

    # Angle checks
    delta.rsub1.rsub2_sub1.sub2 <- rg_rsub1.rsub2_sub1.sub2 - rg_triangle1.sub1.sub2_triangle2.sub1.sub2
    delta.rsub1.con_sub1.con    <- rg_rsub1.con_sub1.con    - rg_triangle1.sub1.con_triangle2.sub1.con
    delta.rsub1.con_sub2.con    <- rg_rsub1.con_sub2.con    - rg_triangle1.sub2.con_triangle2.sub1.con
    delta.rsub2.con_sub1.con    <- rg_rsub2.con_sub1.con    - rg_triangle1.sub1.con_triangle2.sub2.con
    delta.rsub2.con_sub2.con    <- rg_rsub2.con_sub2.con    - rg_triangle1.sub2.con_triangle2.sub2.con

    propSE.rsub1.rsub2_sub1.sub2 <- (rg_rsub1.rsub2_sub1.sub2 - rg_triangle1.sub1.sub2_triangle2.sub1.sub2)  / rg_se_triangle1.sub1.sub2_triangle2.sub1.sub2
    propSE.rsub1.con_sub1.con    <- (rg_rsub1.con_sub1.con    - rg_triangle1.sub1.con_triangle2.sub1.con )   / rg_se_triangle1.sub1.con_triangle2.sub1.con
    propSE.rsub1.con_sub2.con    <- (rg_rsub1.con_sub2.con    - rg_triangle1.sub2.con_triangle2.sub1.con )   / rg_se_triangle1.sub2.con_triangle2.sub1.con
    propSE.rsub2.con_sub1.con    <- (rg_rsub2.con_sub1.con    - rg_triangle1.sub1.con_triangle2.sub2.con )   / rg_se_triangle1.sub1.con_triangle2.sub2.con
    propSE.rsub2.con_sub2.con    <- (rg_rsub2.con_sub2.con    - rg_triangle1.sub2.con_triangle2.sub2.con )   / rg_se_triangle1.sub2.con_triangle2.sub2.con

    # to degs
    a.deg_rsub1.rsub2_sub1.sub2  <- rad2deg(rg2angle(rg_rsub1.rsub2_sub1.sub2))
    a.deg_rsub1.con_sub1.con     <- rad2deg(rg2angle(rg_rsub1.con_sub1.con))
    a.deg_rsub1.con_sub2.con     <- rad2deg(rg2angle(rg_rsub1.con_sub2.con))
    a.deg_rsub2.con_sub2.con     <- rad2deg(rg2angle(rg_rsub2.con_sub2.con))
    a.deg_rsub2.con_sub1.con     <- rad2deg(rg2angle(rg_rsub2.con_sub1.con))

    # Make the list
    double.triangle.output.list <- list(

      # from first part
      triangle1.filename                                                   = triangle1.filename,
      triangle1.folder_location                                            = triangle1.folder_location,
      triangle1.name_sub1                                                  = triangle1.name_sub1,
      triangle1.name_sub2                                                  = triangle1.name_sub2,
      triangle1.name_allcases                                              = triangle1.name_allcases,
      triangle1.name_con                                                   = triangle1.name_con,
      triangle2.filename                                                   = triangle2.filename,
      triangle2.folder_location                                            = triangle2.folder_location,
      triangle2.name_sub1                                                  = triangle2.name_sub1,
      triangle2.name_sub2                                                  = triangle2.name_sub2,
      triangle2.name_allcases                                              = triangle2.name_allcases,
      triangle2.name_con                                                   = triangle2.name_con,

      triangle1.h2_sub1.con                                                          = triangle1.h2_sub1.con,
      triangle1.h2_se_sub1.con                                                       = triangle1.h2_se_sub1.con,
      triangle1.h2_sub2.con                                                          = triangle1.h2_sub2.con,
      triangle1.h2_se_sub2.con                                                       = triangle1.h2_se_sub2.con,
      triangle1.h2_sub1.sub2                                                         = triangle1.h2_sub1.sub2,
     # triangle1.h2_allcases.con                                                      = triangle1.h2_allcases.con,

      triangle1.rg_sub1.con_sub2.con                                                 = triangle1.rg_sub1.con_sub2.con,
      triangle1.rg_sub1.con_sub1.sub2                                                = triangle1.rg_sub1.con_sub1.sub2,
      triangle1.rg_sub2.con_sub1.sub2                                                = triangle1.rg_sub2.con_sub1.sub2,

      triangle1.a.deg_sub1.con_sub2.con                                              = triangle1.a.deg_sub1.con_sub2.con,
      triangle1.a.deg_sub1.con_sub1.sub2                                             = triangle1.a.deg_sub1.con_sub1.sub2,
      triangle1.a.deg_sub2.con_sub1.sub2                                             = triangle1.a.deg_sub2.con_sub1.sub2,

      triangle2.h2_sub1.con                                                          = triangle2.h2_sub1.con,
      triangle2.h2_se_sub1.con                                                       = triangle2.h2_se_sub1.con,
      triangle2.h2_sub2.con                                                          = triangle2.h2_sub2.con,
      triangle2.h2_se_sub2.con                                                       = triangle2.h2_se_sub2.con,
      triangle2.h2_sub1.sub2                                                         = triangle2.h2_sub1.sub2,
      triangle2.h2_allcases.con                                                      = triangle2.h2_allcases.con,

      triangle2.rg_sub1.con_sub2.con                                                 = triangle2.rg_sub1.con_sub2.con,
      triangle2.rg_sub1.con_sub1.sub2                                                = triangle2.rg_sub1.con_sub1.sub2,
      triangle2.rg_sub2.con_sub1.sub2                                                = triangle2.rg_sub2.con_sub1.sub2,

      triangle2.a.deg_sub1.con_sub2.con                                              = triangle2.a.deg_sub1.con_sub2.con,
      triangle2.a.deg_sub1.con_sub1.sub2                                             = triangle2.a.deg_sub1.con_sub1.sub2,
      triangle2.a.deg_sub2.con_sub1.sub2                                             = triangle2.a.deg_sub2.con_sub1.sub2,

      # From 2D-2D part
      plot_2D.2D                                                           = plot_2D.2D,
      triangle1.x.sub1                                                     = triangle1.x.sub1,
      triangle1.y.sub1                                                     = triangle1.y.sub1,
      triangle1.z.sub1                                                     = triangle1.z.sub1,
      triangle1.x.sub2                                                     = triangle1.x.sub2,
      triangle1.y.sub2                                                     = triangle1.y.sub2,
      triangle1.z.sub2                                                     = triangle1.z.sub2,
      triangle2.x.sub1                                                     = triangle2.x.sub1,
      triangle2.y.sub1                                                     = triangle2.y.sub1,
      triangle2.z.sub1                                                     = triangle2.z.sub1,
      triangle2.x.sub2                                                     = triangle2.x.sub2,
      triangle2.y.sub2                                                     = triangle2.y.sub2,
      triangle2.z.sub2                                                     = triangle2.z.sub2,
      triangle1.x.allcases                                                 = triangle1.x.allcases,
      triangle1.y.allcases                                                 = triangle1.y.allcases,
      triangle1.z.allcases                                                 = triangle1.z.allcases,
      triangle1.x.con                                                      = triangle1.x.con,
      triangle1.y.con                                                      = triangle1.y.con,
      triangle1.z.con                                                      = triangle1.z.con,
      triangle2.x.con                                                      = triangle2.x.con,
      triangle2.y.con                                                      = triangle2.y.con,
      triangle2.z.con                                                      = triangle2.z.con,
      triangle1.x.popmean                                                  = triangle1.x.popmean,
      triangle1.y.popmean                                                  = triangle1.y.popmean,
      triangle1.z.popmean                                                  = triangle1.z.popmean,
      triangle2.x.popmean                                                  = triangle2.x.popmean,
      triangle2.y.popmean                                                  = triangle2.y.popmean,
      triangle2.z.popmean                                                  = triangle2.z.popmean,
      rotated                                                              = paste0(flat_subtype, ".con_", rotated_subtype, ".con"),
      angle_difference                                                     = angle_difference,

      rg_triangle1.sub1.sub2_triangle2.sub1.sub2                           = rg_triangle1.sub1.sub2_triangle2.sub1.sub2,
      rg_se_triangle1.sub1.sub2_triangle2.sub1.sub2                        = rg_se_triangle1.sub1.sub2_triangle2.sub1.sub2,
      rg_rsub1.rsub2_sub1.sub2                                             = rg_rsub1.rsub2_sub1.sub2,
      a.deg_rsub1.rsub2_sub1.sub2                                          = a.deg_rsub1.rsub2_sub1.sub2,
      delta.rsub1.rsub2_sub1.sub2                                          = delta.rsub1.rsub2_sub1.sub2,
      propSE.rsub1.rsub2_sub1.sub2                                         = propSE.rsub1.rsub2_sub1.sub2,

      rg_triangle1.sub1.con_triangle2.sub1.con                             = rg_triangle1.sub1.con_triangle2.sub1.con,
      rg_se_triangle1.sub1.con_triangle2.sub1.con                          = rg_se_triangle1.sub1.con_triangle2.sub1.con,
      rg_rsub1.con_sub1.con                                                = rg_rsub1.con_sub1.con,
      a.deg_rsub1.con_sub1.con                                             = a.deg_rsub1.con_sub1.con,
      delta.rsub1.con_sub1.con                                             = delta.rsub1.con_sub1.con,
      propSE.rsub1.con_sub1.con                                            = propSE.rsub1.con_sub1.con,

      rg_triangle1.sub2.con_triangle2.sub1.con                             = rg_triangle1.sub2.con_triangle2.sub1.con,
      rg_se_triangle1.sub2.con_triangle2.sub1.con                          = rg_se_triangle1.sub2.con_triangle2.sub1.con,
      rg_rsub1.con_sub2.con                                                = rg_rsub1.con_sub2.con,
      a.deg_rsub1.con_sub2.con                                             = a.deg_rsub1.con_sub2.con,
      delta.rsub1.con_sub2.con                                             = delta.rsub1.con_sub2.con,
      propSE.rsub1.con_sub2.con                                            = propSE.rsub1.con_sub2.con,

      rg_triangle1.sub1.con_triangle2.sub2.con                             = rg_triangle1.sub1.con_triangle2.sub2.con,
      rg_se_triangle1.sub1.con_triangle2.sub2.con                          = rg_se_triangle1.sub1.con_triangle2.sub2.con,
      rg_rsub2.con_sub1.con                                                = rg_rsub2.con_sub1.con,
      a.deg_rsub2.con_sub1.con                                             = a.deg_rsub2.con_sub1.con,
      delta.rsub2.con_sub1.con                                             = delta.rsub2.con_sub1.con,
      propSE.rsub2.con_sub1.con                                            = propSE.rsub2.con_sub1.con,

      rg_triangle1.sub2.con_triangle2.sub2.con                             = rg_triangle1.sub2.con_triangle2.sub2.con,
      rg_se_triangle1.sub2.con_triangle2.sub2.con                          = rg_se_triangle1.sub2.con_triangle2.sub2.con,
      rg_rsub2.con_sub2.con                                                = rg_rsub2.con_sub2.con,
      a.deg_rsub2.con_sub2.con                                             = a.deg_rsub2.con_sub2.con,
      delta.rsub2.con_sub2.con                                             = delta.rsub2.con_sub2.con,
      propSE.rsub2.con_sub2.con                                            = propSE.rsub2.con_sub2.con)



    } # end if save_data == TRUE



    # Write log file --------------------------------------------------------------

    # Get the time
    rounded_time <- as.POSIXct(round(as.numeric(Sys.time()) / 2) * 2, origin = "1970-01-01")

    # Input
    required_pars <- c(
      "triangle1.h2_sub1.con","triangle1.h2_se_sub1.con", "triangle1.name_sub1","triangle1.N_sub1","triangle1.h2_sub2.con","triangle1.h2_se_sub2.con","triangle1.name_sub2",
      "triangle1.N_sub2", "triangle1.name_allcases","triangle1.name_con", "triangle1.rg_sub1.con_sub2.con", "triangle1.rg_se_sub1.con_sub2.con", #"triangle1.h2_allcases.con","triangle1.h2_se_allcases.con" ,
      "triangle1.plot_title", "triangle1.folder_location","triangle1.filename","triangle1.pop.prev_case",
      "triangle2.h2_sub1.con","triangle2.h2_se_sub1.con", "triangle2.name_sub1","triangle2.N_sub1","triangle2.h2_sub2.con","triangle2.h2_se_sub2.con","triangle2.name_sub2",
      "triangle2.N_sub2", "triangle2.name_allcases","triangle2.name_con", "triangle2.rg_sub1.con_sub2.con", "triangle2.rg_se_sub1.con_sub2.con",# "triangle2.h2_allcases.con","triangle2.h2_se_allcases.con" ,
      "triangle2.plot_title", "triangle2.folder_location","triangle2.filename","triangle2.pop.prev_case",
      "rg_triangle1.sub1.sub2_triangle2.sub1.sub2","rg_se_triangle1.sub1.sub2_triangle2.sub1.sub2","rg_triangle1.sub1.con_triangle2.sub1.con","rg_se_triangle1.sub1.con_triangle2.sub1.con","rg_triangle1.sub2.con_triangle2.sub1.con" ,
      "rg_se_triangle1.sub2.con_triangle2.sub1.con","rg_triangle1.sub1.con_triangle2.sub2.con","rg_se_triangle1.sub1.con_triangle2.sub2.con","rg_triangle1.sub2.con_triangle2.sub2.con","rg_se_triangle1.sub2.con_triangle2.sub2.con"  )

    if (webversion == TRUE) {
      required_pars <- c(
        "triangle1.h2_sub1.con","triangle1.h2_se_sub1.con", "triangle1.name_sub1","triangle1.N_sub1","triangle1.h2_sub2.con","triangle1.h2_se_sub2.con","triangle1.name_sub2",
        "triangle1.N_sub2", "triangle1.name_allcases","triangle1.name_con", "triangle1.rg_sub1.con_sub2.con", "triangle1.rg_se_sub1.con_sub2.con", #"triangle1.h2_allcases.con","triangle1.h2_se_allcases.con" ,
        "triangle1.plot_title", "triangle1.filename","triangle1.pop.prev_case",
        "triangle2.h2_sub1.con","triangle2.h2_se_sub1.con", "triangle2.name_sub1","triangle2.N_sub1","triangle2.h2_sub2.con","triangle2.h2_se_sub2.con","triangle2.name_sub2",
        "triangle2.N_sub2", "triangle2.name_allcases","triangle2.name_con", "triangle2.rg_sub1.con_sub2.con", "triangle2.rg_se_sub1.con_sub2.con",# "triangle2.h2_allcases.con","triangle2.h2_se_allcases.con" ,
        "triangle2.plot_title", "triangle2.filename","triangle2.pop.prev_case",
        "rg_triangle1.sub1.sub2_triangle2.sub1.sub2","rg_se_triangle1.sub1.sub2_triangle2.sub1.sub2","rg_triangle1.sub1.con_triangle2.sub1.con","rg_se_triangle1.sub1.con_triangle2.sub1.con","rg_triangle1.sub2.con_triangle2.sub1.con" ,
        "rg_se_triangle1.sub2.con_triangle2.sub1.con","rg_triangle1.sub1.con_triangle2.sub2.con","rg_se_triangle1.sub1.con_triangle2.sub2.con","rg_triangle1.sub2.con_triangle2.sub2.con","rg_se_triangle1.sub2.con_triangle2.sub2.con"  )
    }

    # Put in coordinates
    "triangle1.coordinates_sub1" <- paste0("(", round(triangle1.x.sub1, 3), ",", round(triangle1.y.sub1, 3), ",",  round(triangle1.z.sub1, 3), ")")
    "triangle1.coordinates_sub2" <- paste0("(", round(triangle1.x.sub2, 3), ",", round(triangle1.y.sub2, 3), ",",  round(triangle1.z.sub2, 3), ")")
    "triangle2.coordinates_sub1" <- paste0("(", round(triangle2.x.sub1, 3), ",", round(triangle2.y.sub1, 3), ",",  round(triangle2.z.sub1, 3), ")")
    "triangle2.coordinates_sub2" <- paste0("(", round(triangle2.x.sub2, 3), ",", round(triangle2.y.sub2, 3), ",",  round(triangle2.z.sub2, 3), ")")
    "triangle1.coordinates_allcases" <- paste0("(", round(triangle1.x.allcases, 3), ",", round(triangle1.y.allcases, 3), ",",  round(triangle1.z.allcases, 3), ")")
    "triangle2.coordinates_allcases" <- paste0("(", round(triangle2.x.allcases, 3), ",", round(triangle2.y.allcases, 3), ",",  round(triangle2.z.allcases, 3), ")")
    "triangle1.coordinates_con" <- paste0("(", round(triangle1.x.con, 3), ",", round(triangle1.y.con, 3), ",",  round(triangle1.z.con, 3), ")")
    "triangle2.coordinates_con" <- paste0("(", round(triangle2.x.con, 3), ",", round(triangle2.y.con, 3), ",",  round(triangle2.z.con, 3), ")")
    "triangle1.coordinates_popmean" <- paste0("(", round(triangle1.x.popmean, 3), ",", round(triangle1.y.popmean, 3), ",",  round(triangle1.z.popmean, 3), ")")
    "triangle2.coordinates_popmean" <- paste0("(", round(triangle2.x.popmean, 3), ",", round(triangle2.y.popmean, 3), ",",  round(triangle2.z.popmean, 3), ")")

    # Round
    a.deg_rsub1.rsub2_sub1.sub2 <- round(a.deg_rsub1.rsub2_sub1.sub2, 0)
    a.deg_rsub1.con_sub1.con    <- round(a.deg_rsub1.con_sub1.con, 0)
    a.deg_rsub1.con_sub2.con    <- round(a.deg_rsub1.con_sub2.con, 0)
    a.deg_rsub2.con_sub1.con    <- round(a.deg_rsub2.con_sub1.con, 0)
    a.deg_rsub2.con_sub2.con    <- round(a.deg_rsub2.con_sub2.con, 0)

    to_logfile_required_input_pars <- required_pars

    to_logfile_output_triangle2D.2D <- c("triangle1.coordinates_sub1", "triangle1.coordinates_sub2", "triangle1.coordinates_allcases", "triangle1.coordinates_con", "triangle1.coordinates_popmean",
                                         "triangle2.coordinates_sub1", "triangle2.coordinates_sub2", "triangle2.coordinates_allcases", "triangle2.coordinates_con", "triangle2.coordinates_popmean",
                                         "a.deg_rsub1.rsub2_sub1.sub2", "a.deg_rsub1.con_sub1.con", "a.deg_rsub1.con_sub2.con", "a.deg_rsub2.con_sub1.con", "a.deg_rsub2.con_sub2.con")

    # Open connection
    if (webversion == FALSE) {
      con_path <- file.path(triangle1.folder_location,paste0(triangle1.filename, ".with.",triangle2.filename,".2D.2D.log.txt")) %>% gsub("^/", "", .)
      con_log <- file(con_path, open = "w")
      }

    if (webversion == TRUE) {
        con_path <- file.path(triangle1.folder_location,paste0(triangle1.filename, ".with.",triangle2.filename,".2D.2D.log.txt")) %>% gsub("^/", "", .)
        con_log <- file(con_path, open = "w")
      }



    # Write general info
    writeLines("............................................................", con_log)
    writeLines(".   GDVIS software created by: A.B. Thijssen               .", con_log)
    writeLines(paste(sprintf(".   Log file created on: %s", rounded_time),"              ."), con_log)
    if(webversion == FALSE) {
     writeLines(paste(".   GDVIS mode used: 2D.2D                                 ."), con_log) }
    if(webversion == TRUE) {
      writeLines(paste(".   GDVIS mode used: 2D.2D online webversion                 ."), con_log)     }
    writeLines("............................................................", con_log)
    writeLines("    ", con_log)
    writeLines("    ", con_log)

    # Write required input data
    writeLines("Note that the 2D.2D mode will also output two individual 2D logfiles. This logfile is only for the 2D.2D relations", con_log)
    writeLines("    ", con_log)
    writeLines("    ", con_log)

    # Write message triangle not possible
    if (save_data == FALSE) {
      writeLines("Error:", con_log)
      writeLines("GDVIS determined that the triangle is not possible, possibly due to inconsistency in input parameters. Please double-check your input. If the issue persists, it may be due to inconsistency in LDSC estimates, for instance because the rg between sub1.con and sub2.con is close to 1", con_log)
      writeLines("    ", con_log)
      writeLines("    ", con_log)
    }

    # Write required input data
    writeLines("Required input data", con_log)
    writeLines(strrep("-", 60), con_log)

    # Iterate over all input parameters and write to the logfile
    for (i in required_pars) {
      name <- i
      value <-  get(i)
      writeLines(sprintf("%-45s %s", name, value), con_log)  }

    # Write updated variables if exceeded [-1,1]
    if (length(log_messages) > 1) {
      writeLines("    ", con_log)
      writeLines(log_messages, con_log)
      writeLines("    ", con_log)
      writeLines("    ", con_log)
    } else {
      writeLines("    ", con_log)
      writeLines("    ", con_log)
    }

    # Write output data
    writeLines("Calculated output data", con_log)
    writeLines(strrep("-", 60), con_log)

    for (i in to_logfile_output_triangle2D.2D) {
      name <- i
      value <-  get(i)
      writeLines(sprintf("%-45s %s", name, value), con_log)  }
    writeLines("    ", con_log)
    writeLines("    ", con_log)


    # Write double checks
    writeLines("Double checks", con_log)
    writeLines(strrep("-", 60), con_log)


    writeLines("Double check 1:", con_log)
    writeLines(paste0("GDVIS rotated ", triangle1.filename, " based on the most precise rg estimate with ", triangle2.filename, " which was ", rg_name, " (", rg_smallest_se_value, " (se ", smallest_se_value, ")"), con_log)
    writeLines(paste0("The second triangle needed to rotate ", round(rotation_angle, 3), " and rotated ", round(final_angle, 3), " which gives a rotation difference of ", round(angle_difference, 3)), con_log)
    writeLines("    ", con_log)

    writeLines("Double check 2:", con_log)
    writeLines(paste0("The provided rg of rg_triangle1.sub1.sub2_triangle2.sub1.sub2 is ", rg_triangle1.sub1.sub2_triangle2.sub1.sub2), con_log)
    writeLines(paste0("After rotation the rg is ", round(rg_rsub1.rsub2_sub1.sub2, 4)), con_log)
    writeLines(paste0("The difference is ", round(delta.rsub1.rsub2_sub1.sub2, 4), " which is ", round(propSE.rsub1.rsub2_sub1.sub2, 4), " * se"), con_log)
    writeLines("    ", con_log)

    writeLines("Double check 3:", con_log)
    writeLines(paste0("The provided rg of rg_triangle1.sub1.con_triangle2.sub1.con is ", rg_triangle1.sub1.con_triangle2.sub1.con), con_log)
    writeLines(paste0("After rotation the rg is ", round(rg_rsub1.con_sub1.con, 4)), con_log)
    writeLines(paste0("The difference is ", round(delta.rsub1.con_sub1.con, 4), " which is ", round(propSE.rsub1.con_sub1.con, 4), " * se"), con_log)
    writeLines("    ", con_log)

    writeLines("Double check 4:", con_log)
    writeLines(paste0("The provided rg of rg_triangle1.sub2.con_triangle2.sub1.con is ", rg_triangle1.sub2.con_triangle2.sub1.con), con_log)
    writeLines(paste0("After rotation the rg is ", round(rg_rsub1.con_sub2.con, 4)), con_log)
    writeLines(paste0("The difference is ", round(delta.rsub1.con_sub2.con, 4), " which is ", round(propSE.rsub1.con_sub2.con, 4), " * se"), con_log)
    writeLines("    ", con_log)

    writeLines("Double check 5:", con_log)
    writeLines(paste0("The provided rg of rg_triangle1.sub1.con_triangle2.sub2.con is ", rg_triangle1.sub1.con_triangle2.sub2.con), con_log)
    writeLines(paste0("After rotation the rg is ", round(rg_rsub2.con_sub1.con, 4)), con_log)
    writeLines(paste0("The difference is ", round(delta.rsub2.con_sub1.con, 4), " which is ", round(propSE.rsub2.con_sub1.con, 4), " * se"), con_log)
    writeLines("    ", con_log)

    writeLines("Double check 6:", con_log)
    writeLines(paste0("The provided rg of rg_triangle1.sub2.con_triangle2.sub2.con is ", rg_triangle1.sub2.con_triangle2.sub2.con), con_log)
    writeLines(paste0("After rotation the rg is ", round(rg_rsub2.con_sub2.con, 4)), con_log)
    writeLines(paste0("The difference is ", round(delta.rsub2.con_sub2.con, 4), " which is ", round(propSE.rsub2.con_sub2.con, 4), " * se"), con_log)
    writeLines("    ", con_log)

    # Close the connection
    close(con_log)


    cli::cli_alert_info(paste0(" For more info see ",triangle1.folder_location,"/",triangle1.filename, ".with.",triangle2.filename,".2D.2D.log.txt"))




  #  return(mget(ls()))



    ### Save ----------------


    if (save_data) {
      save_path <-  file.path(triangle1.folder_location, paste0(triangle1.filename, ".with.",triangle2.filename,".2D_2D.triangle_parameters.RData")) %>% gsub("^/", "", .)
          save(double.triangle.output.list, file = save_path)
          cli::cli_alert_info(file.path(" Data saved as ", triangle1.folder_location, paste0(triangle1.filename, ".with.",triangle2.filename,".2D_2D.triangle_parameters.RData")))
    }

     # Save for webversion
    if (webversion == TRUE & save_data == TRUE) {
      save_path <- file.path(triangle1.folder_location, paste0(triangle1.filename, ".with.",triangle2.filename,".2D_2D.triangle_parameters.RData"))
      save(double.triangle.output.list, file = save_path)
      log_fun("For more info see the logfile", type = "info")
      return(list(
        save_path = save_path,
        log_path = con_path    ))
    }

    cli::cli_alert_success("GDVIS calc 2D.2D succesfully finished!")
    if(webversion == FALSE & plot_2D.2D == TRUE){
      return(save_path)}


    } # End if 2D.2D = T


    #### Calc CD --------------------------------------------------------------------------------------------------------------


    if (plot_CD == TRUE) {

      ### Calculate angles and line lengths of triangle

      # Step 1: Convert input parameters to line lengths and angles
      a.rad_trait1_trait2             <- rg2angle(rg_trait1_trait2)
      a.deg_trait1_trait2             <- rad2deg(a.rad_trait1_trait2)
      a.rad_trait1_trait3             <- rg2angle(rg_trait1_trait3)
      a.deg_trait1_trait3             <- rad2deg(a.rad_trait1_trait3)
      a.rad_trait2_trait3             <- rg2angle(rg_trait2_trait3)
      a.deg_trait2_trait3             <- rad2deg(a.rad_trait2_trait3)
      d_trait1.total                  <- h22d(h2_trait1)
      d_trait2.total                  <- h22d(h2_trait2)
      d_trait3.total                  <- h22d(h2_trait3)
      d_trait1_popmean                <- (1 - pop.prev_trait1) * d_trait1.total
      d_trait2_popmean                <- (1 - pop.prev_trait2) * d_trait2.total
      d_trait3_popmean                <- (1 - pop.prev_trait3) * d_trait3.total
      d_trait1.con_popmean            <- pop.prev_trait1 * d_trait1.total
      d_trait2.con_popmean            <- pop.prev_trait2 * d_trait2.total
      d_trait3.con_popmean            <- pop.prev_trait3 * d_trait3.total

      # Step 2: Calculate trait.trait line lengths
      d_trait1.trait2                 <- f.length_from_angle_and_lines(d_trait1_popmean,d_trait2_popmean,a.rad_trait1_trait2)
      h2_trait1.trait2                <- d2h2(d_trait1.trait2)
      d_trait1.trait3                 <- f.length_from_angle_and_lines(d_trait1_popmean,d_trait3_popmean,a.rad_trait1_trait3)
      h2_trait1.trait3                <- d2h2(d_trait1.trait3)
      d_trait2.trait3                 <- f.length_from_angle_and_lines(d_trait2_popmean,d_trait3_popmean,a.rad_trait2_trait3)
      h2_trait2.trait3                <- d2h2(d_trait2.trait3)

      ### Calculate coordinates of trait1 and trait2
      x.popmean                       <- 0
      y.popmean                       <- 0
      z.popmean                       <- 0
      coord_popmean                   <- c(x.popmean, y.popmean, z.popmean)
      x.trait1                        <- 0
      y.trait1                        <- d_trait1_popmean
      z.trait1                        <- 0
      coord_trait1                    <- c(x.trait1, y.trait1, z.trait1)
      x.trait2                        <- sin(a.rad_trait1_trait2) * d_trait2_popmean
      y.trait2                        <- cos(a.rad_trait1_trait2) * d_trait2_popmean
      z.trait2                        <- 0
      coord_trait2                    <- c(x.trait2, y.trait2, z.trait2)

      # controls trait 1
      x.trait1.con                    <- 0
      y.trait1.con                    <- - d_trait1.con_popmean
      z.trait1.con                    <- 0
      coord_trait1.con                <- c(x.trait1.con,y.trait1.con,z.trait1.con)

      # Controls trait 2
      coord_trait2.con                <- f.extend_line(x.trait2, y.trait2, z.trait2, d_trait2.con_popmean)
      x.trait2.con                    <- coord_trait2.con[1]
      y.trait2.con                    <- coord_trait2.con[2]
      z.trait2.con                    <- coord_trait2.con[3]



      # Function to find best fitting external trait coordinates in the matrix
      test_matrix_fit <- function(xyz_matrix) {

        # Function to find best fitting external trait coordinates
        test_row_fit <- function(x.ext,y.ext,z.ext){
          fit_trait1   <- abs(  {d_trait1.trait3^2}  - { {(x.ext-x.trait1)^2} + {(y.ext-y.trait1)^2} + {(z.ext-z.trait1)^2}  }   )
          fit_trait2   <- abs(  {d_trait2.trait3^2}  - { {(x.ext-x.trait2)^2} + {(y.ext-y.trait2)^2} + {(z.ext-z.trait2)^2}  }   )
          fit_popmean  <- abs(  {d_trait3_popmean^2} - { {(x.ext-x.popmean)^2}  + {(y.ext-y.popmean)^2}  + {(z.ext-z.popmean)^2}   }   )
          return(fit_trait1+fit_trait2+fit_popmean)}

        # Find the differences for each combination of x,y,z and add it to the matrix
        xyz_fit <- tibble::as_tibble(xyz_matrix) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(fit = test_row_fit(x, y, z))

        # Get the lowest fit
        min_fit <- min(xyz_fit$fit)
        best_fit <- xyz_fit %>%
          dplyr::filter(fit == min_fit)

        return(best_fit)      }

      # Function to generate matrix and find best fit
      generate_best_fit <- function(xs, ys, zs) {
        xyz_matrix <- as.matrix(tidyr::crossing(x = xs, y = ys, z = zs))
        best_fit <- test_matrix_fit(xyz_matrix) #%>% unlist()
        return(best_fit)}

      # Initial rough matrix
      coords_rough <- generate_best_fit(seq(from=-1, to=1, by=0.01),
                                        seq(from=-1, to=1, by=0.01),
                                        seq(from=0, to=1, by=0.01))

      # Finer matrix based on rough matrix (was 0.01)
      coords_finer <- generate_best_fit(seq(from = coords_rough['x'] %>% dplyr::pull() - 0.04, to = coords_rough['x'] %>% dplyr::pull() + 0.04, by = 0.001),
                                        seq(from = coords_rough['y'] %>% dplyr::pull() - 0.04, to = coords_rough['y'] %>% dplyr::pull() + 0.04, by = 0.001),
                                        seq(from = max(0, coords_rough['z'] %>% dplyr::pull() - 0.04), to = coords_rough['z'] %>% dplyr::pull() + 0.04, by = 0.001) )

      # Coords finest
      coords_finest <- tibble::tibble()

      for (i in 1:nrow(coords_finer)) {

        # Get one row
        line <- coords_finer[i,]

        # Test the row
        out <- generate_best_fit(seq(from=line['x'] %>% dplyr::pull() - 0.004, to=line['x'] %>% dplyr::pull() + 0.004, by=0.0001),
                                 seq(from=line['y'] %>% dplyr::pull() - 0.004, to=line['y'] %>% dplyr::pull() + 0.004, by=0.0001),
                                 seq(from=line['z'] %>% dplyr::pull() - 0.004, to=line['z'] %>% dplyr::pull() + 0.004, by=0.0001))

        # Add fit to results
        coords_finest <- rbind(coords_finest, out)
      } # end for loop

      # Get the lowest fit
      min_fit <- min(coords_finest$fit)

      # Filter on lowest fit
      coords_finest <- coords_finest %>%
        dplyr::ungroup() %>%  # Remove rowwise structure
        dplyr::filter(fit == min_fit)

      # Lowest fit will sometimes still include two options, then just get the first one
      coords_finest <- coords_finest[1,]

      # Finest matrix based on finer matrix (was 0.002)
      coords_finestest <- generate_best_fit(seq(from=coords_finest['x'] %>% dplyr::pull() - 0.0004, to=coords_finest['x'] %>% dplyr::pull() + 0.0004, by=0.00001),
                                            seq(from=coords_finest['y'] %>% dplyr::pull() - 0.0004, to=coords_finest['y'] %>% dplyr::pull() + 0.0004, by=0.00001),
                                            seq(from=coords_finest['z'] %>% dplyr::pull() - 0.0004, to=coords_finest['z'] %>% dplyr::pull() + 0.0004, by=0.00001))

      # Extract coordinates
      x.trait3 <- coords_finestest %>% dplyr::select(x) %>% dplyr::pull()
      y.trait3 <- coords_finestest %>% dplyr::select(y) %>% dplyr::pull()
      z.trait3 <- coords_finestest %>% dplyr::select(z) %>% dplyr::pull()
      coord_trait3 <- c(x.trait3, y.trait3, z.trait3)

      ## Get coordinates of trait cases to trait controls
      coord_trait3.con                <- f.extend_line(x.trait3, y.trait3, z.trait3, d_trait3.con_popmean)
      x.trait3.con                    <- coord_trait3.con[1]
      y.trait3.con                    <- coord_trait3.con[2]
      z.trait3.con                    <- coord_trait3.con[3]


      ### Internal checks
      if( (f.length_from_coords(coord_popmean,coord_trait1)     -  d_trait1_popmean )      >  (0.02 * d_trait1_popmean) )      {stop("Error code: length check d_trait1_popmean not passed")}
      if( (f.length_from_coords(coord_popmean,coord_trait2)     -  d_trait2_popmean )      >  (0.02 * d_trait2_popmean) )      {stop("Error code: length check d_trait2_popmean not passed")}
      if( (f.length_from_coords(coord_popmean,coord_trait2.con)     -  d_trait2.con_popmean )      >  (0.02 * d_trait2.con_popmean) )      {stop("Error code: length check d_trait2.con_popmean not passed")}
      if( (f.length_from_coords(coord_trait2,coord_trait2.con)     -  d_trait2.total )      >  (0.02 * d_trait2.total) )      {stop("Error code: length check d_trait2.total not passed")}
      if( (f.length_from_coords(coord_popmean,coord_trait1.con)     -  d_trait1.con_popmean )      >  (0.02 * d_trait1.con_popmean) )      {stop("Error code: length check d_trait1.con_popmean not passed")}
      if( (f.length_from_coords(coord_trait1,coord_trait1.con)     -  d_trait1.total )      >  (0.02 * d_trait1.total) )      {stop("Error code: length check d_trait1.total not passed")}
      if( (f.length_from_coords(coord_popmean,coord_trait3.con)     -  d_trait3.con_popmean )      >  (0.02 * d_trait3.con_popmean) )      {stop("Error code: length check d_trait3.con_popmean not passed")}
      if( (f.length_from_coords(coord_trait3,coord_trait3.con)     -  d_trait3.total )      >  (0.02 * d_trait3.total) )      {stop("Error code: length check d_trait3.total not passed")}
      # Internal 3D checks
      internal_check_rg_trait1.con_trait3.con   <- angle2rg(f.angle_from_coords(coord.start.line.a = coord_trait1.con, coord.end.line.a = coord_trait1, coord.start.line.b = coord_trait3.con, coord.end.line.b = coord_trait3))
      internal_check_rg_trait1_trait3           <- angle2rg(f.angle_from_coords(coord.start.line.a = coord_popmean, coord.end.line.a = coord_trait1, coord.start.line.b = coord_popmean, coord.end.line.b = coord_trait3))
      if( internal_check_rg_trait1.con_trait3.con - internal_check_rg_trait1_trait3 > 0.01) {stop("Error code: extension check 1 not passed")}
      if( internal_check_rg_trait1_trait3 - rg_trait1_trait3 > 0.01 ) {stop("Error code: internal check rg_trait1_trait3 not passed")}
      internal_check_rg_trait2.con_trait3.con   <- angle2rg(f.angle_from_coords(coord.start.line.a = coord_trait2.con, coord.end.line.a = coord_trait2, coord.start.line.b = coord_trait3.con, coord.end.line.b = coord_trait3))
      internal_check_rg_trait2_trait3           <- angle2rg(f.angle_from_coords(coord.start.line.a = coord_popmean, coord.end.line.a = coord_trait2, coord.start.line.b = coord_popmean, coord.end.line.b = coord_trait3))
      if( internal_check_rg_trait2.con_trait3.con - internal_check_rg_trait2_trait3 > 0.01) {stop("Error code: extension check 2 not passed")}
      if( internal_check_rg_trait2_trait3 - rg_trait2_trait3 > 0.01 ) {stop("Error code: internal check rg_trait2_trait3 not passed")}
      internal_check_rg_trait1.con_trait2.con   <- angle2rg(f.angle_from_coords(coord.start.line.a = coord_trait1.con, coord.end.line.a = coord_trait1, coord.start.line.b = coord_trait2.con, coord.end.line.b = coord_trait2))
      internal_check_rg_trait1_trait2           <- angle2rg(f.angle_from_coords(coord.start.line.a = coord_popmean, coord.end.line.a = coord_trait1, coord.start.line.b = coord_popmean, coord.end.line.b = coord_trait2))
      if( internal_check_rg_trait1.con_trait2.con - internal_check_rg_trait1_trait2 > 0.01) {stop("Error code: extension check 3 not passed")}
      if( internal_check_rg_trait1_trait2 - rg_trait1_trait2 > 0.01 ) {stop("Error code: internal check rg_trait1_trait2 not passed")}

      # Other h2's
      h2_trait1.cases_trait2.cases <- d2h2(f.length_from_coords(coord_trait1, coord_trait2))
      h2_trait1.cases_trait3.cases <- d2h2(f.length_from_coords(coord_trait1, coord_trait3))
      h2_trait2.cases_trait3.cases <- d2h2(f.length_from_coords(coord_trait2, coord_trait3))

      # All the other rg's
      a.deg_trait1.cases_trait2.cases_vs_trait1.cases_trait3.cases <- rad2deg(f.angle_from_coords(coord_trait1, coord_trait2, coord_trait1, coord_trait3))
      a.deg_trait1.cases_trait2.cases_vs_trait2.cases_trait3.cases <- rad2deg(f.angle_from_coords(coord_trait1, coord_trait2, coord_trait2, coord_trait3))
      a.deg_trait1.cases_trait3.cases_vs_trait2.cases_trait3.cases <- rad2deg(f.angle_from_coords(coord_trait1, coord_trait3, coord_trait2, coord_trait3))
      rg_trait1.cases_trait2.cases_vs_trait1.cases_trait3.cases <- angle2rg(f.angle_from_coords(coord_trait1, coord_trait2, coord_trait1, coord_trait3))
      rg_trait1.cases_trait2.cases_vs_trait2.cases_trait3.cases <- angle2rg(f.angle_from_coords(coord_trait1, coord_trait2, coord_trait2, coord_trait3))
      rg_trait1.cases_trait3.cases_vs_trait2.cases_trait3.cases <- angle2rg(f.angle_from_coords(coord_trait1, coord_trait3, coord_trait2, coord_trait3))



      # Make output list --------------------
      CD.triangle.output.list <- list(

        # input
        name_trait1                                            = name_trait1,
        name_trait2                                            = name_trait2,
        name_trait3                                            = name_trait3,
        name_trait1.con                                        = name_trait1.con,
        name_trait2.con                                        = name_trait2.con,
        name_trait3.con                                        = name_trait3.con,
        rg_trait1_trait2                                       = rg_trait1_trait2,
        rg_trait1_trait3                                       = rg_trait1_trait3,
        rg_trait2_trait3                                       = rg_trait2_trait3,
        h2_trait1                                              = h2_trait1,
        h2_trait2                                              = h2_trait2,
        h2_trait3                                              = h2_trait3,
        pop.prev_trait1                                        = pop.prev_trait1,
        pop.prev_trait2                                        = pop.prev_trait2,
        pop.prev_trait3                                        = pop.prev_trait3,

        # Calculated
        a.rad_trait1_trait2                                     = a.rad_trait1_trait2,
        a.deg_trait1_trait2                                     = a.deg_trait1_trait2,
        a.rad_trait1_trait3                                     = a.rad_trait1_trait3 ,
        a.deg_trait1_trait3                                     = a.deg_trait1_trait3 ,
        a.rad_trait2_trait3                                     = a.rad_trait2_trait3 ,
        a.deg_trait2_trait3                                     = a.deg_trait2_trait3,
        d_trait1.total                                          = d_trait1.total,
        d_trait2.total                                          = d_trait2.total,
        d_trait3.total                                          = d_trait3.total,
        d_trait1_popmean                                        = d_trait1_popmean,
        d_trait2_popmean                                        = d_trait2_popmean,
        d_trait3_popmean                                        = d_trait3_popmean,
        d_trait1.con_popmean                                    = d_trait1.con_popmean,
        d_trait2.con_popmean                                    = d_trait2.con_popmean,
        d_trait3.con_popmean                                    = d_trait3.con_popmean ,
        d_trait1.trait2                                         = d_trait1.trait2,
        h2_trait1.trait2                                        = h2_trait1.trait2,
        d_trait1.trait3                                         = d_trait1.trait3,
        h2_trait1.trait3                                        = h2_trait1.trait3,
        d_trait2.trait3                                         = d_trait2.trait3,
        h2_trait2.trait3                                        = h2_trait2.trait3,
        x.popmean                                               = x.popmean,
        y.popmean                                               = y.popmean,
        z.popmean                                               = z.popmean,
        x.trait1                                                = x.trait1,
        y.trait1                                                = y.trait1,
        z.trait1                                                = z.trait1,
        x.trait2                                                = x.trait2,
        y.trait2                                                = y.trait2,
        z.trait2                                                = z.trait2,
        x.trait1.con                                            = x.trait1.con,
        y.trait1.con                                            = y.trait1.con,
        z.trait1.con                                            = z.trait1.con,
        x.trait2.con                                            = x.trait2.con,
        y.trait2.con                                            = y.trait2.con,
        z.trait2.con                                            = z.trait2.con ,
        x.trait3                                                = x.trait3,
        y.trait3                                                = y.trait3,
        z.trait3                                                = z.trait3,
        x.trait3.con                                            = x.trait3.con,
        y.trait3.con                                            = y.trait3.con,
        z.trait3.con                                            = z.trait3.con,
        h2_trait1.cases_trait2.cases                            = h2_trait1.cases_trait2.cases,
        h2_trait1.cases_trait3.cases                            = h2_trait1.cases_trait3.cases,
        h2_trait2.cases_trait3.cases                            = h2_trait2.cases_trait3.cases,
        rg_trait1.cases_trait2.cases_vs_trait1.cases_trait3.cases = rg_trait1.cases_trait2.cases_vs_trait1.cases_trait3.cases,
        rg_trait1.cases_trait2.cases_vs_trait2.cases_trait3.cases = rg_trait1.cases_trait2.cases_vs_trait2.cases_trait3.cases,
        rg_trait1.cases_trait3.cases_vs_trait2.cases_trait3.cases = rg_trait1.cases_trait3.cases_vs_trait2.cases_trait3.cases,
        a.deg_trait1.cases_trait2.cases_vs_trait1.cases_trait3.cases = a.deg_trait1.cases_trait2.cases_vs_trait1.cases_trait3.cases,
        a.deg_trait1.cases_trait2.cases_vs_trait2.cases_trait3.cases = a.deg_trait1.cases_trait2.cases_vs_trait2.cases_trait3.cases,
        a.deg_trait1.cases_trait3.cases_vs_trait2.cases_trait3.cases = a.deg_trait1.cases_trait3.cases_vs_trait2.cases_trait3.cases,
        folder_location                                         = folder_location,
        plot_CD                                                 = plot_CD

      )


      # Save the list
      save_path <- file.path(folder_location, paste0(name_trait1, ".with.",name_trait1,".with.",name_trait3,".CD.triangle_parameters.RData")) %>% gsub("^/", "", .)
      save(CD.triangle.output.list, file = save_path)
      cli::cli_alert_info(paste0(" Data saved as ", folder_location,"/",name_trait1, ".with.",name_trait1,".with.",name_trait3,".CD.triangle_parameters.RData"  ))


      ## Make logfile ----------------------------------------------------------------

      # Get info that needs to be written

      # Change coordinates
      "coordinates_popmean" <- paste0("(", round(x.popmean, 3), ",", round(y.popmean, 3), ",",  round(z.popmean, 3), ")")
      "coordinates_trait1" <- paste0("(", round(x.trait1, 3), ",", round(y.trait1, 3), ",",  round(z.trait1, 3), ")")
      "coordinates_trait2" <- paste0("(", round(x.trait2, 3), ",", round(y.trait2, 3), ",",  round(z.trait2, 3), ")")
      "coordinates_trait3" <- paste0("(", round(x.trait3, 3), ",", round(y.trait3, 3), ",",  round(z.trait3, 3), ")")
      "coordinates_trait1.con" <- paste0("(", round(x.trait1.con, 3), ",", round(y.trait1.con, 3), ",",  round(z.trait1.con, 3), ")")
      "coordinates_trait2.con" <- paste0("(", round(x.trait2.con, 3), ",", round(y.trait2.con, 3), ",",  round(z.trait2.con, 3), ")")
      "coordinates_trait3.con" <- paste0("(", round(x.trait3.con, 3), ",", round(y.trait3.con, 3), ",",  round(z.trait3.con, 3), ")")

      # Input parameters for CD
      to_logfile_input_pars_CD <- c(
        "name_trait1",  "name_trait2", "name_trait3",  "name_trait1.con", "name_trait2.con", "name_trait3.con", "rg_trait1_trait2", "rg_trait1_trait3",  "rg_trait2_trait3",
        "h2_trait1","h2_trait2","h2_trait3","pop.prev_trait1","pop.prev_trait2","pop.prev_trait3", "folder_location")
      if (webversion == TRUE) {
        to_logfile_input_pars_CD <- c(
          "name_trait1",  "name_trait2", "name_trait3",  "name_trait1.con", "name_trait2.con", "name_trait3.con", "rg_trait1_trait2", "rg_trait1_trait3",  "rg_trait2_trait3",
          "h2_trait1","h2_trait2","h2_trait3","pop.prev_trait1","pop.prev_trait2","pop.prev_trait3")}

      # Output parameters for CD
      to_logfile_output_pars_CD <- c(
        "a.rad_trait1_trait2",   "a.deg_trait1_trait2", "a.rad_trait1_trait3","a.deg_trait1_trait3","a.rad_trait2_trait3","a.deg_trait2_trait3","d_trait1.total", "d_trait2.total",
        "d_trait3.total",   "d_trait1_popmean","d_trait2_popmean", "d_trait3_popmean","d_trait1.con_popmean", "d_trait2.con_popmean","d_trait3.con_popmean",  "d_trait1.trait2", "h2_trait1.trait2","d_trait1.trait3",
        "h2_trait1.trait3", "d_trait2.trait3", "h2_trait2.trait3", "h2_trait1.cases_trait2.cases", "h2_trait1.cases_trait3.cases", "h2_trait2.cases_trait3.cases","rg_trait1.cases_trait2.cases_vs_trait1.cases_trait3.cases",
        "rg_trait1.cases_trait2.cases_vs_trait2.cases_trait3.cases","rg_trait1.cases_trait3.cases_vs_trait2.cases_trait3.cases"        )

      # Get the time
      rounded_time <- as.POSIXct(round(as.numeric(Sys.time()) / 2) * 2, origin = "1970-01-01")

      # Open connection
      con_path <- file.path(folder_location, paste0(name_trait1,"_", name_trait2,"_", name_trait3,".CD.log.txt")) %>% gsub("^/", "", .)
      con_log <- file(con_path, open = "w")
      cli::cli_alert_info(paste0(" For more info see ",folder_location,"/", name_trait1,"_", name_trait2,"_", name_trait3,".CD.log.txt"))

      # Write general info
      writeLines("............................................................", con_log)
      writeLines(".   GDVIS software created by: A.B. Thijssen               .", con_log)
      writeLines(paste(sprintf(".   Log file created on: %s", rounded_time),"              ."), con_log)
      if (webversion == FALSE) {
        writeLines(paste(".   GDVIS mode used: cross-disorder                    ."), con_log) }
      if (webversion == TRUE) {
        writeLines(paste(".   GDVIS mode used: cross-disorder webversion             ."), con_log) }
      writeLines("............................................................", con_log)
      writeLines("    ", con_log)
      writeLines("    ", con_log)
      writeLines("    ", con_log)

      writeLines("Note that no significance checks are done, make sure to only visualize sensible data.", con_log)
      writeLines("    ", con_log)
      writeLines("    ", con_log)

      # Write required input data
      writeLines("Required input data", con_log)
      writeLines(strrep("-", 60), con_log)

      # Iterate over all input parameters and write to the logfile
      for (i in to_logfile_input_pars_CD) {
        name <- i
        value <-  get(i)
        writeLines(sprintf("%-45s %s", name, value), con_log)  }
      writeLines("    ", con_log)

      # Write output data
      writeLines("Calculated output data", con_log)
      writeLines(strrep("-", 60), con_log)

      for (i in to_logfile_output_pars_CD) {
        name <- i
        value <-  get(i)
        if (grepl("coordinates", name)) {
          rounded_value <- value    } else {
            rounded_value <- round(as.numeric(value), 3)      }
        writeLines(sprintf("%-45s %s", name, rounded_value), con_log)  }
      writeLines("    ", con_log)
      writeLines("    ", con_log)


      # Close the connection
      close(con_log)

      if(webversion == FALSE) {
        cli::cli_alert_success("GDVIS calc succesfully finished CD calculations!")
        return(save_path) }

      if(webversion == TRUE) {
        log_fun("GDVIS calc successfully finished!", type = "success")
        return(list(
          save_path = save_path,
          log_path = con_path,
          log_path = "Hello"))      }



    } # End calc CD


  }) # End with env

}  # End function
