# Plot GDVIS plot

This function plots a GDVIS visualization based on GDVIS parameters.

## Usage

``` r
GDVIS_plot(
  input_triangle_parameters,
  show.rendering = TRUE,
  show.names = TRUE,
  x_lower = NULL,
  x_upper = NULL,
  y_lower = NULL,
  y_upper = NULL
)
```

## Arguments

- input_triangle_parameters:

  the path to the RData GDVIS output, which ends with
  3D.triangle_parameters.RData

- show.rendering:

  option to show the rendering of the plot in an external window,
  default is TRUE

- show.names:

  option to show the names of the groups, default is TRUE

- x_lower:

  change axes to make the plots all the same size or make the legend
  appear outside of the plot (only for 2D plot)

- x_upper:

  change axes to make the plots all the same size or make the legend
  appear outside of the plot (only for 2D plot)

- y_lower:

  change axes to make the plots all the same size or make the legend
  appear outside of the plot (only for 2D plot)

- y_upper:

  change axes to make the plots all the same size or make the legend
  appear outside of the plot (only for 2D plot)
