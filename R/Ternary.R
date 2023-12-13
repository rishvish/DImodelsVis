#' @title Ternary diagram with contours
#'
#' @description
#' The data preparation function to create a gradient of proportions of three
#' entities where their proportions sum to 1, thereby representing a 2 dimensional
#' simplex which can be shown on a ternary surface. The projections of the three
#' entities on the x-y plane is also calculated. The final output can then used with a
#' relevant statistical model to make predictions for the response and visualise it
#' as the proportions of the three entities change over the simplex (response surface
#' analysis). The output of this function can be passed to the
#' \code{\link{ternary_plot}} function to plot the results.
#'
#' @param prop A character vector specifying the columns names of entities
#'             whose proportions to manipulate. Default is ".P1", ".P2",
#'             and ".P3".
#' @param x A character string specifying the name for the column containing
#'          the x component of the x-y projection of the simplex. Default is ".x".
#' @param y A character string specifying the name for the column containing
#'          the y component of the x-y projection of the simplex. Default is ".y".
#' @param add_var A list specifying values for additional variables
#'                in the model other than the proportions (i.e. not part of the simplex design).
#'                This would be useful to compare the predictions across
#'                different values for a categorical variable.
#'                One plot will be generated for each unique combination
#'                of values specified here.
#' @param resolution A number between 1 and 5 describing the resolution of the
#'                   resultant graph.
#'                   A high value would result in a higher definition figure
#'                   but at the cost of being computationally expensive.
#' @param prediction A logical value indicating whether to pass the final data
#'                   to the `\link{add_prediction}` function and append the
#'                   predictions to the data. Default value is \code{TRUE}, but
#'                   often it would be desirable to make additional changes to
#'                   the data before making any predictions, so the user can set this to
#'                   \code{FALSE} and manually call the `\link{add_prediction}`
#'                   function.
#' @inheritDotParams add_prediction -data
#'
#' @return A data-frame with the following columns and any additional columns
#'         specified in `add_var` parameter
#'  \describe{
#'    \item{.x}{The x component of the x-y projection of the simplex point.}
#'    \item{.y}{The y component of the x-y projection of the simplex point.}
#'    \item{.P1}{The first entity whose proportion is varied across the simplex.}
#'    \item{.P2}{The second entity whose proportion is varied across the simplex.}
#'    \item{.P3}{The third entity whose proportion is varied across the simplex.}
#'    \item{.add_str_ID}{An identifier column for grouping the cartesian product
#'                       of all additional columns specified in `add_var`
#'                       parameter (if `add_var` is specified).}
#'    \item{.Pred}{The predicted response for each observation
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
#' ## Can also add any additional variables independent of the simplex
#' sim0$treatment <-  rep(c("A", "B", "C", "D"), each = 16)
#' new_mod <- update(mod, ~. + treatment, data = sim0)
#' head(ternary_data(prop = c("p1", "p2", "p3"),
#'                   add_var = list("treatment" = c("A", "B")),
#'                   resolution = 1, model = new_mod))
#'
#' ## It could be desirable to take the output of this function and add
#' ## additional variables to the data before making predictions
#' ## Use `prediction = FALSE` to get data without any predictions
#' contour_data <- ternary_data(prop = c("p1", "p2", "p3"),
#'                              model = mod,
#'                              prediction = FALSE,
#'                              resolution = 1)
#' head(contour_data)
#'
#' ## Manually add the treatment variable
#' contour_data$treatment <- "A"
#' ## Make predictions
#' head(add_prediction(data = contour_data, model = new_mod))#'
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
#'
#' ## Note: Add predictions via coefficients would not give confidence intervals
#' ## to get CIs using coefficients we need to specify the variance-covariance
#' ## matrix using `vcov`
#' contour_data <- add_prediction(data = contour_data,
#'                                coefficient = mod$coefficient,
#'                                vcov = vcov(mod),
#'                                interval = "confidence")
#' head(contour_data)
ternary_data <- function(prop = c(".P1", ".P2", ".P3"),
                         x = ".x", y = ".y",
                         add_var = list(),
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
  triangle <- subset(triangle, (((base*sin(pi/3)*2) >= high) &
                                  (((1-base)*sin(pi/3)*2) >= high)))

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
    add_var <- check_add_var(model = model, add_var = add_var)
  }

  # Add any experimental structures
  if(length(add_var) > 0){
    triangle <- add_add_var(add_var = add_var, data = triangle)
  }

  if(prediction){
    triangle <- add_prediction(data = triangle, ...)
  }

  # Order columns
  triangle <- triangle %>% select(all_of(c(prop, x, y)), everything())

  return(triangle)
}

#' @title Ternary diagram
#'
#' @description
#' The plotting function to visualise response change over a 2-d simplex
#' surface (ternary diagram). Create the ternary surface using the
#' `\link{ternary_data}` function and add predictions of the
#' response using the `\link{add_prediction}`. This function can then be used to visualise response change over the
#' ternary surface as a contour map. One could also visualise raw data-points
#' across the simplex by projecting them onto the ternary surface by using
#' `\link{prop_to_tern_proj}` function and using the `show = "points"`
#' parameter. See examples for more details.
#'
#' @param data A data-frame consisting of the x-y plane projection of the
#'             2-d simplex. This data could be the output of the
#'             `\link{ternary_data}` function, and contain the  predicted
#'             response at each point along the simplex to show the variation
#'             in response as a contour map.
#' @param x The column name or index containing the x-component of the x-y
#'          projection of the 2-d simplex.
#' @param y The column name or index containing the y-component of the x-y
#'          projection of the 2-d simplex.
#' @param tern_labels A character vector containing the labels of the vertices
#'                    of the ternary. The default is the column names of the
#'                    first three columns of the data, with the first column
#'                    corresponding to the top vertex, second column corresponding
#'                    to the left vertex and the third column corresponding to
#'                    the right vertex of the ternary.
#' @param pred The column name or index containing the predicted response
#'             (for showing contours).
#' @param show A character string indicating whether to show raw points or contours
#'             on the ternary. The default is to show "contours".
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
#' @param points_size If showing points, then a numeric value specifying the size
#'                    of the points.
#' @param nlevels The number of levels to show on the contour map.
#' @param colours A character vector or function specifying the colours for the
#'                contour map. The number of colours should be same as `nlevels`.
#'                The default colours scheme is the
#'                \code{\link[grDevices:terrain.colors]{terrain.colors()}}.
#' @param lower_lim A number to set a custom lower limit for the contour
#'                  (if `show = "contours"`). The default is minimum of the prediction.
#' @param upper_lim A number to set a custom upper limit for the contour
#'                  (if `show = "contours"`). The default is maximum of the prediction.
#' @param contour_text A boolean value indicating whether to include labels on
#'                     the contour lines showing their values
#'                     (if `show = "contours"`). The default is \code{TRUE}.
#' @param nrow Number of rows in which to arrange the final plot
#'             (when `add_var` is specified).
#' @param ncol Number of columns in which to arrange the final plot
#'             (when `add_var` is specified).
#'
#' @inherit prediction_contributions return
#'
#' @export
#'
#' @examples
#' library(DImodels)
#' library(dplyr)
#' library(ggplot2)
#'
#' ## Load data
#' data(sim0)
#'
#' ### Show raw data in ternary
#' ## Project the 3-d simplex data into x-y plane
#' raw_data <-  prop_to_tern_proj(data = sim0, prop = c("p1", "p2" , "p3"))
#'
#' ## ternary_plot shows contours by default, use `show = "points"` to show
#' ## points across the ternary
#' ternary_plot(data = raw_data, show = "points")
#'
#' ## Show axis guides and increase points size
#' ternary_plot(data = raw_data, show = "points",
#'              show_axis_guides = TRUE, points_size = 4)
#'
#' ## Overlay the output of the plot with geom_point to customise the points
#' ## further
#' ternary_plot(data = raw_data, show = "points",
#'              show_axis_guides = TRUE) +
#'              geom_point(aes(colour = response))
#'
#' ### Show contours of response
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
#' ## Change number of contours and set custom upper and lower
#' ## limits for the scale
#' ternary_plot(plot_data, nlevels = 10, colours = hcl.colors(10),
#'              lower_lim = 10, upper_lim = 35)
#'
#' ## Add axis guides
#' ternary_plot(plot_data, show_axis_guides = TRUE)
#'
#' ## Change ternary labels along with their font-size
#' ternary_plot(plot_data, tern_labels = c("Sp1", "Sp2", "Sp3"),
#'              vertex_label_size = 6, axis_label_size = 5)
#'
#' ## Add additional variables and create a separate plot for each
#' sim0$treatment <-  rep(c("A", "B", "C", "D"), each = 16)
#' new_mod <- update(mod, ~. + treatment, data = sim0)
#' tern_data <- ternary_data(resolution = 1, model = new_mod,
#'                           prop = c("p1", "p2", "p3"),
#'                           add_var = list("treatment" = c("A", "C")))
#' ## Arrange plot in 2 columns
#' ternary_plot(data = tern_data, ncol = 2)
ternary_plot <- function(data,
                         x = ".x",
                         y = ".y",
                         pred = ".Pred",
                         show = c("contours", "points"),
                         tern_labels = c("P1", "P2", "P3"),
                         show_axis_labels = TRUE,
                         show_axis_guides = FALSE,
                         axis_label_size = 4,
                         vertex_label_size = 5,
                         points_size = 2,
                         nlevels = 7,
                         colours = NULL,
                         lower_lim = NULL,
                         upper_lim = NULL,
                         contour_text = TRUE,
                         nrow = 0,
                         ncol = 0){
  if(missing(data)){
    cli::cli_abort(c("{.var data} cannot be empty.",
                     "i" = "Specify a data-frame or tibble
                     preferably the output of {.fn visualise_effects_data} or
                     a data-frame with a similar structure and column names."))
  }
  show <-  match.arg(show)
  # Multiple plots if additional variables were specified
  if(check_col_exists(data, ".add_str_ID")){
    ids <- unique(data$.add_str_ID)
    plots <- lapply(cli_progress_along(1:length(ids), name = "Creating plot",
                                       format = paste0(
                                         "{pb_spin} Creating plot ",
                                         "[{pb_current}/{pb_total}]   ETA:{pb_eta}"
                                       )),
                    function(i){
                      data_iter <- data %>% filter(.data$.add_str_ID == ids[i])
                      ternary_plot_internal(data = data_iter,
                                            x = x,
                                            y = y,
                                            pred = pred,
                                            show = show,
                                            points_size = points_size,
                                            nlevels = nlevels,
                                            colours = colours,
                                            lower_lim = lower_lim,
                                            upper_lim = upper_lim,
                                            tern_labels = tern_labels,
                                            contour_text = contour_text,
                                            show_axis_labels = show_axis_labels,
                                            show_axis_guides = show_axis_guides,
                                            axis_label_size = axis_label_size,
                                            vertex_label_size = vertex_label_size)+
                        labs(subtitle = ids[i])
                    })
    if(length(plots) > 1){
      plot <- new("ggmultiplot", plots = plots, nrow = nrow, ncol = ncol)
    } else {
      plot <- plots[[1]]
    }
    cli::cli_alert_success("Created all plots.")
  # Single plot otherwise
  } else {
    plot <- ternary_plot_internal(data = data,
                                  x = x,
                                  y = y,
                                  pred = pred,
                                  show = show,
                                  points_size = points_size,
                                  nlevels = nlevels,
                                  colours = colours,
                                  lower_lim = lower_lim,
                                  upper_lim = upper_lim,
                                  tern_labels = tern_labels,
                                  contour_text = contour_text,
                                  show_axis_labels = show_axis_labels,
                                  show_axis_guides = show_axis_guides,
                                  axis_label_size = axis_label_size,
                                  vertex_label_size = vertex_label_size)
    cli::cli_alert_success("Created plot.")
  }
  return(plot)
}

#' @keywords internal
#' Internal function for creating a ternary plot
#'
#' @usage NULL
NULL
ternary_plot_internal <- function(data,
                                  x = ".x",
                                  y = ".y",
                                  pred = ".Pred",
                                  show = c("contours", "points"),
                                  nlevels = 7,
                                  colours = NULL,
                                  lower_lim = NULL,
                                  upper_lim = NULL,
                                  tern_labels = c("P1", "P2", "P3"),
                                  show_contours = TRUE,
                                  contour_text = TRUE,
                                  show_axis_labels = TRUE,
                                  show_axis_guides = FALSE,
                                  points_size = 2,
                                  axis_label_size = 4,
                                  vertex_label_size = 5){
    # Check all inputs are appropriate. Print informative error messages if not
    sanity_checks(data = data,
                  characters = list("tern_labels" = tern_labels),
                  numerics = list("nlevels" = nlevels,
                                  "axis_label_size" = axis_label_size,
                                  "vertex_label_size" = vertex_label_size,
                                  "points_size" = points_size),
                  booleans = list("contour_text" = contour_text,
                                  "show_axis_guides" = show_axis_guides,
                                  "show_axis_labels" = show_axis_labels,
                                  "show_contours" = show_contours),
                  colours = colours)
    show <- match.arg(show)

    # Ensure x, y and predictions are present in data
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
    if(show == "contours"){
      check_presence(data = data, col = pred,
                     message = c("The column name/index specified in {.var pred} is
                               not present in the data.",
                                 "i" = "Specify the name/index of the column
                                     containing the predictions for the
                                     communities in the simplex in {.var pred}."))

      # If user didn't specify lower limit assume it to be min of predicted response
      if(is.null(lower_lim)){
        lower_lim <- round(min(data[, pred]), 2)
      }

      # If user didn't specify upper limit assume it to be max of predicted response
      if(is.null(upper_lim)){
        upper_lim <- round(max(data[, pred]), 2)
      }
    } else {
      lower_lim <- upper_lim <- 0
    }


    # Further checks for appropriateness of parameters
    sanity_checks(numerics = list("upper_lim" = upper_lim,
                                  "lower_lim" = lower_lim),
                  unit_lengths = list("upper_lim" = upper_lim,
                                      "lower_lim" = lower_lim,
                                      "x" = x, "y" = y, "pred" = pred,
                                      "nlevels" = nlevels,
                                      "points_size" = points_size,
                                      "axis_label_size" = axis_label_size,
                                      "vertex_label_size" = vertex_label_size,
                                      "contour_text" = contour_text,
                                      "show_contours" = show_contours,
                                      "show_axis_guides" = show_axis_guides,
                                      "show_axis_labels" = show_axis_labels))

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


    if(show == "contours"){
      # Create colour-scale (legend) for plot
      # Create breaks between range of legend
      size <- nlevels + 1
      breaks <- round(seq(lower_lim, upper_lim, length.out= size), 2)

      # Choose colours for contour
      # If user didn't specify colours then use the default terrain colours
      if(is.null(colours)){
        colours <- terrain.colors(nlevels, rev = T)
      }
      pl <- ggplot(data, aes(x = .data[[x]], y = .data[[y]],
                             z = .data[[pred]])) +
        geom_raster(aes(fill = .data[[pred]]))+
        scale_fill_stepsn(colours = colours, breaks = breaks,
                          labels = function(val){
                            val
                          },
                          limits = c(lower_lim, upper_lim),
                          show.limits = T)+
        geom_contour(breaks = breaks, colour = 'black')+
        guides(fill = guide_colorbar(frame.colour = 'black',
                                     ticks.colour = 'black',
                                     title = 'Predicted\nResponse',
                                     show.limits = T))+
        theme_void()+
        theme(legend.key.size = unit(0.1, 'npc'),
              legend.key.height = unit(0.04, 'npc'),
              legend.title = element_text(size = 14, vjust = 1.25),
              plot.subtitle = element_text(hjust=0.5, size=14),
              strip.text = element_text(size =14, vjust = 0.5),
              legend.text = element_text(size = 12, angle = 45,
                                         vjust = 1.2, hjust = 1.2),
              legend.position = 'bottom')
    } else {
      # Base of plot
      pl <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
        theme_void()
    }

    # Showing axis labels
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

    # Showing axis guides
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

    # Layering plot
    pl <- pl +
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
      #facet_wrap(~ Value, ncol = length(values))+
      coord_fixed()

    # Show points
    if(show == "points"){
      pl <- pl + geom_point(size = points_size)
    }

    if(show == "contours" && contour_text){
      pl <- pl +
        geom_text_contour(skip=0, breaks = breaks,
                          label.placer = metR::label_placer_fraction(0.15),
                          size=3.5, nudge_x = 0.015, nudge_y = 0.015)

    }
    return(pl)
  }
