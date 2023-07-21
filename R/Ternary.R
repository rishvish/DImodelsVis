#' @rdname Ternary
#' @title Create a 2-d simplex (Ternary) diagram showing change in
#'        response over the ternary surface
#'
#' @description
#' The data preparation function to create a gradient of proportions of three
#' entities representing a 2 dimensional simplex (the proportions sum to 1)
#' representing a ternary surface along with the mapping of each point within
#' the simplex on to the x-y plane. The final output can then used with a
#' relevant statistical model to make predictions for the response as the
#' proportions of the three entities change over the simplex.
#'
#' @param prop A character vector specifying the columns names of entities
#'             whose proportions to manipulate. Default is ".P1", ".P2",
#'             and ".P3".
#' @param x A character string specifying the name for the column containing
#'          the x component of the x-y projection of the simplex.
#' @param y A character string specifying the name for the column containing
#'          the y component of the x-y projection of the simplex.
#' @param exp_str A list specifying values for additional experimental
#'                structures in the model other than the proportions.
#'                This would be useful to compare the predictions across
#'                different values for a categorical variable.
#'                One plot will be generated for each unique combination
#'                of values specified here.
#' @param resolution A number between 1 and 5 describing the resolution of the
#'                   resultant graph.
#'                   A high value would result in a higher definition figure
#'                   but at the cost of being computationally expensive.
#' @param prediction A logical value indicating whether to pass the final data
#'                   to `\link{add_prediction}` and add predictions to the data.
#'                   Default value is \code{TRUE}, but often it would be
#'                   desirable to make additional changes to the data before
#'                   making any predictions, so the user can set this to
#'                   \code{FALSE} and manually call the `\link{add_prediction}`
#'                   function.
#' @inheritDotParams add_prediction -data
#'
#' @return A data-frame with the following columns and any additional columns
#'         specified in `exp_str` parameter
#'  \describe{
#'    \item{.x}{The x component of the x-y projection of the simplex point.}
#'    \item{.y}{The y component of the x-y projection of the simplex point.}
#'    \item{.P1}{The first entity whose proportion is varied across the simplex.}
#'    \item{.P2}{The second entity whose proportion is varied across the simplex.}
#'    \item{.P3}{The third entity whose proportion is varied across the simplex.}
#'    \item{.add_str_ID}{An identifier column for grouping the cartesian product
#'                       of all additional columns specified in `exp_str`
#'                       parameter (if `exp_str` is specified).}
#'    \item{.Pred}{The predicted repsonse for each observation
#'                (if `prediction` is \code{TRUE}).}
#'    \item{.Lower}{The lower limit of the prediction/confidence interval
#'                  for each observation.}
#'    \item{.Upper}{The upper limit of the prediction/confidence interval
#'                  for each observation.}
#'  }
#'
#' @export
#'
#' @examples
#' library(DImodels)
#' library(dplyr)
#'
#' ## Load data
#' data(sim0)
#'
#' ## Fit model
#' mod <- lm(response ~ 0 + (p1 + p2 + p3)^2, data = sim0)
#'
#' ## Create a contour map of predicted response over the ternary surface
#' ## Remember to specify prop with the same character values as the names
#' ## of the variables in the model containing the prop.
#' head(ternary_data(resolution = 1, model = mod,
#'                   prop = c("p1", "p2", "p3")))
#'
#' ## Can also add any additional experimental structures
#' head(ternary_data(prop = c("p1", "p2", "p3"),
#'                   exp_str = list("richness" = c("A", "B")),
#'                   resolution = 1,
#'                   prediction = FALSE))
#'
#' ## It could be desirable to take the output of this function and add
#' ## additional variables to the data before making predictions
#' ## Use `prediction = FALSE` to get data without any predictions
#' contour_data <- ternary_data(prop = c("p1", "p2", "p3"),
#'                   model = mod,
#'                   prediction = FALSE,
#'                   resolution = 1)
#' head(contour_data)
#'
#' ## Manually add the interaction terms
#' contour_data <- contour_data %>%
#'                   mutate(`p1:p2` = p1*p2,
#'                          `p2:p3` = p2*p3,
#'                          `p1:p3` = p1*p3)
#'
#' ## Add predictions using model coefficients
#' contour_data <- add_prediction(data = contour_data,
#'                                coefficient = mod$coefficient)
#' head(contour_data)
#' ## Note: Add predictions via coefficients would not give confidence intervals
ternary_data <- function(prop = c(".P1", ".P2", ".P3"),
                         x = ".x", y = ".y",
                         exp_str = list(),
                         resolution = 3, prediction = TRUE, ...){
  # Ensure inputs are proper
  sanity_checks(characters = list("prop" = prop, "x" = x, "y" = y),
                numerics = list("resolution" = resolution),
                booleans = list("prediction" = prediction),
                unit_lengths = list("resolution" = resolution, "x" = x,
                                    "prediction" = prediction, "y" = y))

  if(!isTRUE(all.equal(resolution, as.integer(resolution))) ||
                                     !between(resolution, 1, 10)){
    cli::cli_warn(c("{.var resolution} should be a whole number with values
                    between 1 and 10",
                    "i" = "The value specified for {.var resolution} was
                          {as.character(resolution)}.",
                    "i" = "Reverting back to the default value of 3."))
  }

  # Prepare data for creating conditional proportions
  base <- seq(0,1,l=100*2*resolution)
  high <- seq(0,sin(pi/3),l=86.6*2*resolution)
  triangle <- expand.grid(base = base, high = high)
  triangle <- subset(triangle, (((base*sin(pi/3)*2) > high) &
                                  (((1-base)*sin(pi/3)*2) > high)))

  # Extrapolate 2-d simplex coordinates to represent proportions
  # of the three species shown in the simplex
  triangle <- triangle %>%
    mutate(!! prop[1] := .data$high*2/sqrt(3),
           !! prop[3] := .data$base - .data$high/sqrt(3),
           !! prop[2] := 1 - .data$high*2/sqrt(3) -
                                    (.data$base - .data$high/sqrt(3))) %>%
    rename(!!x := base, !!y := high) %>%
    select(all_of(c(x, y, prop[1:3])))

  # Ensure experimental structure are specified correctly
  dotArgs <- rlang::dots_values(...)
  model <- if (!is.null(dotArgs$model)) dotArgs$model else NULL
  if(!is.null(model)){
    exp_str <- check_exp_str(model = model, exp_str = exp_str)
  }

  # Add any experimental structures
  if(length(exp_str) > 0){
    triangle <- add_exp_str(exp_str = exp_str, data = triangle)
  }

  if(prediction){
    triangle <- add_prediction(data = triangle, ...)
  }

  return(triangle)
}

#' @rdname Ternary
#' @title Create a 2-d simplex (Ternary) diagram showing change in response
#'        over the ternary surface
#'
#' @description
#' The plotting function to visualise response change over a 2-d simplex
#' surface (ternary diagram). Create the ternary surface using the
#' `\link{ternary_data}` function and add predictions of the
#' response using the `\link{add_prediction}`.
#' This function can then be used to visualise response change over the
#' ternary surface as a contour map.
#'
#' @param data A data-frame which is the output of the
#'             `\link{ternary_data}` function, consisting of the x-y
#'             plane projection of the 2-d simplex along with the predicted
#'             response.
#' @param x The column name or index containing the x-component of the x-y
#'          projection of the 2-d simplex.
#' @param y The column name or index containing the y-component of the x-y
#'          projection of the 2-d simplex.
#' @param pred The column name or index containing the predicted response.
#' @param nlevels The number of levels to show on the contour map.
#' @param colours A character vector or function specifying the colours for the
#'                contour map. The number of colours should be same as `nlevels`.
#'                The default colours scheme is the
#'                \code{\link[grDevices:terrain.colors]{terrain.colors()}}.
#' @param lower_lim A number to set a custom lower limit for the contour.
#'                  The default is minimum of the prediction.
#' @param upper_lim A number to set a custom upper limit for the contour.
#'                  The default is maximum of the prediction.
#' @param tern_labels A character vector containing the labels of the vertices
#'                    of the ternary. The default is "P1", "P2", and "P3".
#' @param contour_text A boolean value indicating whether to include labels on
#'                     the contour lines showing their values. The default is
#'                     \code{TRUE}.
#' @param show_axis_labels A boolean value indicating whether to show axis
#'                         labels along the edges of the ternary. The default
#'                         is \code{TRUE}.
#' @param show_axis_guides A boolean value indicating whether to show axis
#'                         guides within the interior of the ternary. The
#'                         default is \code{FALSE}.
#' @param axis_label_size A numeric value to adjust the size of the axis labels
#'                        in the ternary plot. The default size is 4.
#' @param vertex_label_size A numeric value to adjust the size of the vertex
#'                          labels in the ternary plot. The default size is 5.
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' library(DImodels)
#' library(dplyr)
#'
#' ## Load data
#' data(sim0)
#'
#' ## Fit model
#' mod <- lm(response ~ 0 + (p1 + p2 + p3)^2, data = sim0)
#'
#' ## Create a contour map of predicted response over the ternary surface
#' ## Remember to specify prop with the same character values as the names
#' ## of the variables in the model containing the prop.
#' plot_data <- ternary_data(resolution = 1, model = mod,
#'                           prop = c("p1", "p2", "p3"))
#'
#' ## Create a contour plot of response across the ternary space
#' ternary_plot(plot_data)
#'
#' ## Change colour scheme
#' cols <- hcl.colors(7) # because there are 7 contour levels by default
#' ternary_plot(plot_data, colours = cols)
#'
#' # Change number of contours and set custom upper and lower
#' # limits for the scale
#' ternary_plot(plot_data, nlevels = 10, colours = hcl.colors(10),
#'              lower_lim = 10, upper_lim = 35)
#'
#' # Add axis guides
#' ternary_plot(plot_data, show_axis_guides = TRUE)
#'
#' # Change ternary labels along with their font-size
#' ternary_plot(plot_data, tern_labels = c("Sp1", "Sp2", "Sp3"),
#'              vertex_label_size = 6, axis_label_size = 5)
ternary_plot <- function(data,
                         x = ".x",
                         y = ".y",
                         pred = ".Pred",
                         nlevels = 7,
                         colours = NULL,
                         lower_lim = NULL,
                         upper_lim = NULL,
                         tern_labels = c("P1", "P2", "P3"),
                         contour_text = TRUE,
                         show_axis_labels = TRUE,
                         show_axis_guides = FALSE,
                         axis_label_size = 4,
                         vertex_label_size = 5){
  # Check all inputs are appropriate. Print informative error messages if not
  sanity_checks(data = data,
                characters = list("tern_labels" = tern_labels),
                numerics = list("nlevels" = nlevels,
                                "axis_label_size" = axis_label_size,
                                "vertex_label_size" = vertex_label_size),
                booleans = list("contour_text" = contour_text,
                                "show_axis_guides" = show_axis_guides,
                                "show_axis_labels" = show_axis_labels),
                colours = colours)

  # Ensure x, y and predictions are present in data
  check_presence(data = data, col = pred,
                 message = c("The column name/index specified in {.var pred} is
                             not present in the data.",
                             "i" = "Specify the name/index of the column
                                   containing the predictions for the
                                   communities in the simplex in {.var pred}."))
  check_presence(data = data, col = x,
                 message = c("The column name/index specified in {.var x} is
                             not present in the data.",
                             "i" = "Specify the name/index of the column
                                   containing the x-coordinates of the 2-d
                                   projection of the communities in the simplex
                                   in {.var x}."))
  check_presence(data = data, col = y,
                 message = c("The column name/index specified in {.var y} is
                             not present in the data.",
                             "i" = "Specify the name/index of the column
                                   containing the y-coordinates of the 2-d
                                   projection of the communities in the simplex
                                   in {.var y}."))

  # If user didn't specify lower limit assume it to be min of predicted response
  if(is.null(lower_lim)){
    lower_lim <- round(min(data[, pred]), 2)
  }

  # If user didn't specify upper limit assume it to be max of predicted response
  if(is.null(upper_lim)){
    upper_lim <- round(max(data[, pred]), 2)
  }

  # Further checks for appropriateness of parameters
  sanity_checks(numerics = list("upper_lim" = upper_lim,
                                "lower_lim" = lower_lim),
                unit_lengths = list("upper_lim" = upper_lim,
                                    "lower_lim" = lower_lim,
                                    "x" = x, "y" = y, "pred" = pred,
                                    "nlevels" = nlevels,
                                    "axis_label_size" = axis_label_size,
                                    "vertex_label_size" = vertex_label_size,
                                    "contour_text" = contour_text,
                                    "show_axis_guides" = show_axis_guides,
                                    "show_axis_labels" = show_axis_labels))

  # Create colour-scale (legend) for plot
  # Create breaks between range of legend
  size <- nlevels + 1
  breaks <- round(seq(lower_lim, upper_lim, length.out= size), 2)

  # Choose colours for contour
  # If user didn't specify colours then use the default terrain colours
  if(is.null(colours)){
    colours <- terrain.colors(nlevels, rev = T)
  }

  # Labels for the ternary
  if(length(tern_labels) !=3){
    if(length(tern_labels) > 3){
      cli::cli_warn(c("More than three labels were specified for the ternary
                      diagram. The first three labels in {.var tern_labels}
                      will be chosen"))
    } else {
      cli::cli_abort(c("Three labels are needed for the ternary, only
                       {length(tern_labels)} were specified."))
    }
  }

  # Create plot
  pl <- ggplot(data, aes(x = .data[[x]], y = .data[[y]],
                         z = .data[[pred]]))+
    geom_raster(aes(fill = .data[[pred]]))+
    scale_fill_stepsn(colours = colours, breaks = breaks,
                      labels = function(val){
                        val
                      },
                      limits = c(lower_lim, upper_lim),
                      show.limits = T)+
    geom_contour(breaks = breaks, colour = 'black')+
    geom_text(data = tibble(x = c(0.5, 0, 1), y = c(sqrt(3)/2, 0,  0),
                            label = tern_labels,
                            !! pred := 0),
              aes(x= .data$x, y= .data$y, label = .data$label),
              size = vertex_label_size, fontface='plain',
              nudge_x = c(0, -0.05, 0.05),
              nudge_y = c(0.05, 0, 0))+
    geom_segment(data = tibble(x = c(0, 0, 1), y = c(0,0,0),
                               xend = c(1, 0.5, 0.5),
                               yend = c(0, sqrt(3)/2, sqrt(3)/2),
                               !! pred := 0),
                 aes(x=.data$x, y=.data$y, xend=.data$xend, yend=.data$yend),
                 linewidth = 1)+
    theme_void()+
    guides(fill = guide_colorbar(frame.colour = 'black',
                                 ticks.colour = 'black',
                                 title = 'Predicted\nResponse',
                                 show.limits = T))+
    #facet_wrap(~ Value, ncol = length(values))+
    coord_fixed()+
    theme(legend.key.size = unit(0.1, 'npc'),
          legend.key.height = unit(0.04, 'npc'),
          legend.title = element_text(size = 14, vjust = 0.75),
          plot.subtitle = element_text(hjust=0.5, size=14),
          strip.text = element_text(size =14, vjust = 0.5),
          legend.text = element_text(size = 12, angle = 45,
                                     vjust = 1.2, hjust = 1.2),
          legend.position = 'bottom')


  if(show_axis_labels){
    # Labels for the ternary axes
    axis_labels <- tibble(x1 = seq(0.2,0.8,0.2),
                          y1 = c(0,0,0,0),
                          x2 = .data$x1/2,
                          y2 = .data$x1*sqrt(3)/2,
                          x3 = (1-.data$x1)*0.5+.data$x1,
                          y3 = sqrt(3)/2-.data$x1*sqrt(3)/2,
                          label = .data$x1,
                          rev_label = rev(.data$label),
                          !! pred := 0)

    pl <- pl +
      geom_text(data = axis_labels,
                aes(x=.data$x1, y=.data$y1, label=.data$label),
                nudge_y=-0.055, size = axis_label_size)+
      geom_text(data = axis_labels,
                aes(x=.data$x2, y=.data$y2, label=.data$rev_label),
                nudge_x=-0.055, nudge_y=0.055, size = axis_label_size)+
      geom_text(data = axis_labels,
                aes(x=.data$x3, y=.data$y3, label=.data$rev_label),
                nudge_x=0.055, nudge_y=0.055, size = axis_label_size)
  }

  if(show_axis_guides){
    pl <- pl +
      geom_segment(data = axis_labels,
                   aes(x = .data$x1, y = .data$y1,
                       xend = .data$x2, yend = .data$y2), colour='grey',
                   linetype='dashed', linewidth=1, alpha = .75)+
      geom_segment(data = axis_labels,
                   aes(x = .data$x1, y = .data$y1,
                       xend = .data$x3, yend = .data$y3), colour='grey',
                   linetype='dashed', linewidth=1, alpha = .75)+
      geom_segment(data = axis_labels,
                   aes(x = .data$x2, y = .data$y2,
                       xend = rev(.data$x3), yend = rev(.data$y3)),
                   colour='grey',
                   linetype='dashed', linewidth=1, alpha = .75)
  }

  if(contour_text){
    pl <- pl +
      geom_text_contour(skip=0, breaks = breaks,
                        label.placer = metR::label_placer_fraction(0.15),
                        size=3.5, nudge_x = 0.015, nudge_y = 0.015)

  }
  return(pl)
}
