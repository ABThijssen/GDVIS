devtools::load_all()

## 2D -------------------------------------------------------------------

input.list <- list(
  h2_sub1.con                                   = 0.1606,
  h2_se_sub1.con                                = 0.0139,
  name_sub1                                     = "With CT",
  N_sub1                                        = 14146,
  h2_sub2.con                                   = 0.0839,
  h2_se_sub2.con                                = 0.0098,
  name_sub2                                     = "Without CT",
  N_sub2                                        = 16285,
  name_allcases                                 = "all MDD cases",
  name_con                                      = "controls",
  rg_sub1.con_sub2.con                          = 0.6011,
  rg_se_sub1.con_sub2.con                       = 0.0610,
  plot_title                                    = "Childhood Trauma",
  filename                                      = "Childhood Trauma",
  pop.prev_case                                 = 0.16,

  # All variables below can optionally be added to the input.list
  optional_LDSC_rg_allcases.con_sub1.con        = 0.9111,
  optional_LDSC_rg_se_allcases.con_sub1.con     = 0.0167,
  optional_LDSC_rg_allcases.con_sub2.con        = 0.8765,
  optional_LDSC_rg_se_allcases.con_sub2.con     = 0.0198,
  optional_LDSC_h2_sub1.sub2                    = 0.1036,
  optional_LDSC_h2_se_sub1.sub2                 = 0.0162,
  optional_LDSC_rg_sub1.con_sub1.sub2           = 0.6995,
  optional_LDSC_rg_se_sub1.con_sub1.sub2        = 0.0389,
  optional_LDSC_rg_sub2.con_sub1.sub2           = -0.1544,
  optional_LDSC_rg_se_sub2.con_sub1.sub2        = 0.0983,
  optional_LDSC_h2_allcases.con                 = 0.0915,
  optional_LDSC_h2_se_allcases.con              = 0.0079)

output <- GDVIS::GDVIS_calc(input.list)
GDVIS::GDVIS_plot(output, x_lower = -0.35, x_upper = 0.35, y_lower = -0.15, y_upper = 0.35)
GDVIS::GDVIS_plot_2D(output)


## 3D -------------------------------------------------------------------

input.list.3D <- list(

  # 2D required input
  h2_sub1.con                                   = 0.1606,
  h2_se_sub1.con                                = 0.0139,
  name_sub1                                     = "With CT",
  N_sub1                                        = 14146,
  h2_sub2.con                                   = 0.0839,
  h2_se_sub2.con                                = 0.0098,
  name_sub2                                     = "Without CT",
  N_sub2                                        = 16285,
  name_allcases                                 = "all MDD cases",
  name_con                                      = "controls",
  rg_sub1.con_sub2.con                          = 0.6011,
  rg_se_sub1.con_sub2.con                       = 0.0610,
  plot_title                                    = "Childhood Trauma",
  filename                                      = "Childhood Trauma",
  pop.prev_case                                 = 0.16,

  # 2D optional input
  optional_LDSC_rg_allcases.con_sub1.con        = 0.9111,
  optional_LDSC_rg_se_allcases.con_sub1.con     = 0.0167,
  optional_LDSC_rg_allcases.con_sub2.con        = 0.8765,
  optional_LDSC_rg_se_allcases.con_sub2.con     = 0.0198,
  optional_LDSC_h2_sub1.sub2                    = 0.1036,
  optional_LDSC_h2_se_sub1.sub2                 = 0.0162,
  optional_LDSC_rg_sub1.con_sub1.sub2           = 0.6995,
  optional_LDSC_rg_se_sub1.con_sub1.sub2        = 0.0389,
  optional_LDSC_rg_sub2.con_sub1.sub2           = -0.1544,
  optional_LDSC_rg_se_sub2.con_sub1.sub2        = 0.0983,
  optional_LDSC_h2_allcases.con                 = 0.0915,
  optional_LDSC_h2_se_allcases.con              = 0.0079,

  # 3D required input
  plot_3D                                       = TRUE,
  h2_ext                                        = 0.0738,
  h2_se_ext                                     = 0.0276,
  rg_sub1.con_ext                               = 0.5883,
  rg_se_sub1.con_ext                            = 0.1629,
  rg_sub2.con_ext                               = 0.6665,
  rg_se_sub2.con_ext                            = 0.1733,
  name_ext                                      = "Anxiety",
  pop.prev_ext                                  = 0.2,

  # 3D optional input
  optional_LDSC_rg_sub1.sub2_ext                = 0.1198,
  optional_LDSC_rg_se_sub1.sub2_ext             = 0.1797,
  optional_LDSC_rg_allcases.con_ext             = 0.6942,
  optional_LDSC_rg_se_allcases.con_ext          = 0.1733)


output.3D <- GDVIS_calc(input.list.3D)
GDVIS_plot(output.3D)



## 2D.2D -------------------------------------------------------------------

input.list.2D.2D <- list(

  # Subgroup info triangle 1
  triangle1.h2_sub1.con                                 = 0.1606,
  triangle1.h2_se_sub1.con                              = 0.0139,
  triangle1.name_sub1                                   = "With childhood trauma",
  triangle1.N_sub1                                      = 14146,
  triangle1.h2_sub2.con                                 = 0.0839,
  triangle1.h2_se_sub2.con                              = 0.0098,
  triangle1.name_sub2                                   = "without childhood trauma",
  triangle1.N_sub2                                      = 16285,
  triangle1.name_allcases                               = "all MDD cases",
  triangle1.name_con                                    = "MDD controls",
  triangle1.rg_sub1.con_sub2.con                        = 0.6011,
  triangle1.rg_se_sub1.con_sub2.con                     = 0.0610,

  # General info triangle 1
  triangle1.plot_title                                  = "Childhood Trauma",
  triangle1.filename                                    = "Childhood Trauma",
  triangle1.pop.prev_case                               = 0.16,

  # Optional input triangle 1
  triangle1.optional_LDSC_rg_allcases.con_sub1.con      = 0.9111,
  triangle1.optional_LDSC_rg_se_allcases.con_sub1.con   = 0.0167,
  triangle1.optional_LDSC_rg_allcases.con_sub2.con      = 0.8765,
  triangle1.optional_LDSC_rg_se_allcases.con_sub2.con   = 0.0198,
  triangle1.optional_LDSC_h2_sub1.sub2                  = 0.1036,
  triangle1.optional_LDSC_h2_se_sub1.sub2               = 0.0162,
  triangle1.optional_LDSC_rg_sub1.con_sub1.sub2         = 0.6995,
  triangle1.optional_LDSC_rg_se_sub1.con_sub1.sub2      = 0.0389,
  triangle1.optional_LDSC_rg_sub2.con_sub1.sub2         = -0.1544,
  triangle1.optional_LDSC_rg_se_sub2.con_sub1.sub2      = 0.0983,
  triangle1.optional_LDSC_h2_allcases.con               = 0.0915,
  triangle1.optional_LDSC_h2_se_allcases.con            = 0.0079,

  # Subgroup info triangle 2
  triangle2.h2_sub1.con                                 = 0.1394,
  triangle2.h2_se_sub1.con                              = 0.0147,
  triangle2.name_sub1                                   = "With anxiety",
  triangle2.N_sub1                                      = 11885,
  triangle2.h2_sub2.con                                 = 0.0798,
  triangle2.h2_se_sub2.con                              = 0.0098,
  triangle2.name_sub2                                   = "without anxiety",
  triangle2.N_sub2                                      = 18546,
  triangle2.name_allcases                               = "all MDD cases",
  triangle2.name_con                                    = "MDD controls",
  triangle2.rg_sub1.con_sub2.con                        = 0.8451,
  triangle2.rg_se_sub1.con_sub2.con                     = 0.0720,

  # General info triangle 2
  triangle2.plot_title                                  = "Comorbid anxiety",
  triangle2.filename                                    = "Comorbid anxiety",
  triangle2.pop.prev_case                               = 0.16,

  # Optional input triangle 2
  triangle2.optional_LDSC_rg_allcases.con_sub1.con      = 0.9504,
  triangle2.optional_LDSC_rg_se_allcases.con_sub1.con   = 0.0235,
  triangle2.optional_LDSC_rg_allcases.con_sub2.con      = 0.9686,
  triangle2.optional_LDSC_rg_se_allcases.con_sub2.con   = 0.0153,
  triangle2.optional_LDSC_h2_sub1.sub2                  = 0.0450,
  triangle2.optional_LDSC_h2_se_sub1.sub2               = 0.0166,
  triangle2.optional_LDSC_rg_sub1.con_sub1.sub2         = 0.6503,
  triangle2.optional_LDSC_rg_se_sub1.con_sub1.sub2      = 0.0752,
  triangle2.optional_LDSC_rg_sub2.con_sub1.sub2         = 0.1302,
  triangle2.optional_LDSC_rg_se_sub2.con_sub1.sub2      = 0.1870,
  triangle2.optional_LDSC_h2_allcases.con               = 0.0915,
  triangle2.optional_LDSC_h2_se_allcases.con            = 0.0079,

  # Required data on triangle relations
  rg_triangle1.sub1.sub2_triangle2.sub1.sub2            = 0.4084,
  rg_se_triangle1.sub1.sub2_triangle2.sub1.sub2         = 0.1793,
  rg_triangle1.sub1.con_triangle2.sub1.con              = 0.9139,
  rg_se_triangle1.sub1.con_triangle2.sub1.con           = 0.0326,
  rg_triangle1.sub2.con_triangle2.sub1.con              = 0.7831,
  rg_se_triangle1.sub2.con_triangle2.sub1.con           = 0.0431,
  rg_triangle1.sub1.con_triangle2.sub2.con              = 0.8452,
  rg_se_triangle1.sub1.con_triangle2.sub2.con           = 0.0316,
  rg_triangle1.sub2.con_triangle2.sub2.con              = 0.8882,
  rg_se_triangle1.sub2.con_triangle2.sub2.con           = 0.0328,

  # Type info
  plot_2D.2D                                            = TRUE)

output.2D.2D <- GDVIS_calc(input.list.2D.2D)
GDVIS_plot(output.2D.2D)




## CD ------------------------------------------------------------------------

input.list.CD <- list(
  plot_CD                                       = TRUE,
  name_trait1                                   = "SCZ_cases",
  name_trait2                                   = "BIP_cases",
  name_trait3                                   = "MDD_cases",
  name_trait1.con                               = "SCZ_con",
  name_trait2.con                               = "BIP_con",
  name_trait3.con                               = "MDD_con",
  rg_trait1_trait2                              = 0.6986,
  rg_trait1_trait3                              = 0.3598,
  rg_trait2_trait3                              = 0.4758,
  h2_trait1                                     = 0.4118,
  h2_trait2                                     = 0.2902,
  h2_trait3                                     = 0.0703,
  pop.prev_trait1                               = 0.004,
  pop.prev_trait2                               = 0.01,
  pop.prev_trait3                               = 0.16)

output.CD <- GDVIS_calc(input.list.CD)
GDVIS_plot(output.CD, show.names = F)
