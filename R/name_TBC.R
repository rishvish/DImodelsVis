#' @title
#' Creating data for visualising the change in a response variable between two
#' points in the simplex space
#'
#' @description
#' This is the helper function to prepare the underlying data for visualising
#' the change in a response variable between two points in a simplex space. The
#' two points specified by the `starts` and `ends` parameters are joined by a
#' straight line across the simplex space and the response is predicted for the
#' starting, ending and intermediate communities along this line. The associated
#' uncertainty along this prediction is also returned. The output of this function
#' can be passed to the \code{\link{suggest_name_plot}} function to visualise the
#' change in response.
#'
#' @param starts A data-frame specifying the starting proportions of the
#'               compositional variables.
#'               If a model object is specified then this data should contain all the
#'               variables present in the model object including any experimental structures.
#'               If a coefficient vector is specified then data should contain same number of
#'               columns as the number of elements in the coefficient vector and a one-to-one
#'               positional mapping would be assumed between the data columns and the
#'               elements of the coefficient vector.
#' @param ends A data-frame specifying the ending proportions of the
#'             compositional variables.
#'             If a model object is specified then this data should contain all the
#'             variables present in the model object including any experimental structures.
#'             If a coefficient vector is specified then data should contain same number of
#'             columns as the number of elements in the coefficient vector and a one-to-one
#'             positional mapping would be assumed between the data columns and the
#'             elements of the coefficient vector.
#' @param prop A vector of column names or indices identifying the columns containing the
#'             variable proportions (i.e., compositional columns) in the data.
#' @inheritParams ternary_data
#' @inheritDotParams add_prediction -data
#'
#' @return A data frame with the following columns appended at the end
#'  \describe{
#'    \item{.InterpConst}{The value of the interpolation constant for creating
#'                        the intermediate compositions between the start and end compositions.}
#'    \item{.Group}{An identifier column to discern between the different curves.}
#'    \item{.add_str_ID}{An identifier column for grouping the cartesian product
#'                       of all additional columns specified in `add_var`
#'                       parameter (if `add_var` is specified).}
#'    \item{.Pred}{The predicted response for each observation.}
#'    \item{.Lower}{The lower limit of the prediction/confidence interval for each observation.}
#'    \item{.Upper}{The upper limit of the prediction/confidence interval for each observation.}
#'  }
#'
#' @export
#'
#' @examples
#' library(DImodels)
#'
#' ## Load data
#' data(sim2)
#'
#' ## Fit model
#' mod <- glm(response ~ (p1 + p2 + p3 + p4)^2 + 0, data = sim2)
#'
#' ## Create data for visualising change in response as we move from
#' ## a species dominated by 70% of one species to a monoculture of
#' ## same species
#' head(suggest_name_data(starts = sim2[c(1, 5, 9, 13), 3:6],
#'                        ends = sim2[c(48, 52, 56, 60), 3:6],
#'                        prop = c("p1", "p2", "p3", "p4"),
#'                        model = mod))
#'
#' ## Create data for visualising change in response as we move from
#' ## the centroid mixture to each monoculture
#' ## If either of starts or ends have only row, then they'll be recycled
#' ## to match the number of rows in the other
#' ## Notice starts has only one row here, but will be recycled to have 4
#' ## since ends has 4 four rows
#' head(suggest_name_data(starts = sim2[c(18),3:6],
#'                        ends = sim2[c(48, 52, 56, 60),3:6],
#'                        prop = c("p1", "p2", "p3", "p4"),
#'                        model = mod))
#'
#' ## Changing the confidence level for the prediction interval
#' ## Use `conf.level` parameter
#' head(suggest_name_data(starts = sim2[c(18), 3:6],
#'                        ends = sim2[c(48, 52, 56, 60),3:6],
#'                        prop = c("p1", "p2", "p3", "p4"),
#'                        model = mod, conf.level = 0.99))
#'
#' ## Adding additional variables to the data using `add_var`
#' ## Notice the new .add_str_ID column in the output
#' sim2$block <- as.numeric(sim2$block)
#' new_mod <- update(mod, ~ . + block, data = sim2)
#' head(suggest_name_data(starts = sim2[c(18), 3:6],
#'                        ends = sim2[c(48, 52, 56, 60), 3:6],
#'                        prop = c("p1", "p2", "p3", "p4"),
#'                        model = new_mod, conf.level = 0.99,
#'                        add_var = list("block" = c(1, 2))))
#'
#' ## Use predict = FALSE to get raw data structure
#' out_data <- suggest_name_data(starts = sim2[c(18), 3:6],
#'                               ends = sim2[c(48, 52, 56, 60), 3:6],
#'                               prop = c("p1", "p2", "p3", "p4"),
#'                               model = new_mod,
#'                               prediction = FALSE)
#' head(out_data)
#' ## Manually add block
#' out_data$block = 3
#' ## Call `add_prediction` to get prediction
#' head(add_prediction(data = out_data, model = new_mod, interval = "conf"))
suggest_name_data <- function(starts, ends, prop,
                              add_var = list(),
                              prediction = TRUE, ...){
  if(missing(starts)){
    cli::cli_abort(c("{.var starts} cannot be empty.",
                     "i" = "Specify a data frame or tibble indicating the initial variable proportions in {.var data}."))
  }

  if(missing(ends)){
    cli::cli_abort(c("{.var ends} cannot be empty.",
                     "i" = "Specify a data frame or tibble indicating the final variable proportions in {.var data}."))
  }

  sanity_checks(data = starts, prop = prop,
                booleans = list("prediction" = prediction))
  sanity_checks(data = ends, prop = prop)

  # Ensure starts and ends have the same number of rows
  if(nrow(starts) != nrow(ends)){
    if(nrow(starts) != 1 & nrow(ends) != 1){
      cli::cli_abort(c("The number of rows of the data in {.var starts} and {.var ends}
                       should be the same or atleast one of {.var starts} or
                       {.var ends} should have only one row.",
                       "i" = "{.var starts} has {nrow(starts)} rows while
                       {.var ends} has {nrow(ends)} rows."))
    }
  }

  if(ncol(starts) != ncol(ends) | !isTRUE(all.equal(colnames(starts), colnames(ends)))){
    all_comm <- union(colnames(starts), colnames(ends))
    int_comm <- intersect(colnames(starts), colnames(ends))
    cli::cli_abort(c("The data in {.var starts} and {.var ends} should have identical columns.",
                     "i" = "The following column{?s} don't match between {.var starts}
                     and {.var ends} {.val {all_comm[which(!all_comm %in% int_comm)]}}"))
  }

  # If any of starts or ends has 1 row then expand
  if(nrow(starts) == 1){
    starts <- starts %>% slice(rep(1, times = nrow(ends)))
  }
  if(nrow(ends) == 1){
    ends <- ends %>% slice(rep(1, times = nrow(starts)))
  }

  # Ensure there's a one-to-one mapping between columns of starts and ends
  ends <- ends %>% select(all_of(colnames(starts)))

  # Get names of columns containing species proportions
  species_names <- colnames(starts[, prop])

  pvals <- seq(0, 100, 1)/100
  plot_data <- lapply(cli_progress_along(seq_len(nrow(starts)), name = "Preparing data"), function(i){
    # Interpolate data between the starting and ending communities
    interpolated_data <- interpolate_communities(starts[i, species_names],
                                                 ends[i, species_names],
                                                 species_names) %>%
      left_join(y = starts %>% select(-all_of(species_names)) %>% slice(i),
                by = character()) %>%
      mutate(".InterpConst" = pvals,
             ".Group" = i) %>%
      select(all_of(species_names), everything())

    # Slice data or proportions will be added as list
    interpolated_data[1:nrow(interpolated_data), ]
  }) %>% bind_rows()

  # Ensure experimental structure are specified correctly
  dotArgs <- rlang::dots_values(...)
  model <- if (!is.null(dotArgs$model)) dotArgs$model else NULL
  if(!is.null(model)){
    add_var <- check_add_var(model = model, add_var = add_var)
  }

  # Add any experimental structures
  if(length(add_var) > 0){
    plot_data <- add_add_var(add_var = add_var, data = plot_data)
  }

  # Make prediction and get marginal effect
  if(prediction){
    dots <- list(...)
    dots$data <- plot_data
    dots$interval <- if (is.null(dots$interval)) "conf" else dots$interval
    plot_data <- do.call(add_prediction, as.list(dots))

    # Calculate the marginal effect of adding a species for producing the marginal plots
    # plot_data <- plot_data %>% group_by(.data$.Group) %>%
    #   mutate(.dy = c((diff(.data$.Pred)/diff(.data$.Proportion)), 1)) %>%
    #   mutate('.Marginal' = c(.data$.dy[1:(length(.data$.dy) - 1)], .data$.dy[(length(.data$.dy) - 1)]),
    #          '.Threshold' = .data$.Proportion[abs(.data$.Marginal) == min(abs(.data$.Marginal))][1],
    #          '.MarEffect' = ifelse(!!sym(".Proportion") < .data$.Threshold, 'Negative', 'Positive')) %>%
    #   select(-.data$.dy) %>%
    #   ungroup()
  }

  cli::cli_alert_success("Finished data preparation.")
  return(plot_data)
}


#' @title
#' Visualising the change in a response variable between two points in
#' the simplex space
#'
#' @description
#' The helper function for plotting the change in a response variable over a
#' straight line between two points across the simplex space. The output of the
#' \code{\link{suggest_name_data}} function (with any desired modifications)
#' should be passed here. The generated plot will show individual curves
#' indicating the variation in the response between the points.
#' `\code{\link[PieGlyph:PieGlyph-package]{Pie-glyphs}}` are
#' used to highlight the compositions of the starting, ending and midpoint of the
#' straight line between the two points.
#'
#' @param data A data frame created using the \code{\link{suggest_name_data}} function.
#' @inheritParams visualise_effects_plot
#' @inheritParams prediction_contributions
#'
#' @inherit prediction_contributions return
#'
#' @export
#'
#' @examples
#' library(DImodels)
#'
#' ## Load data
#' data(sim2)
#'
#' ## Fit model
#' mod <- glm(response ~ (p1 + p2 + p3 + p4)^2 + 0, data = sim2)
#'
#' ## Visualise change as we move from the centroid community to each monoculture
#' plot_data <- suggest_name_data(starts = sim2[c(19, 20, 19, 20), ],
#'                                ends = sim2[c(47, 52, 55, 60), ],
#'                                prop = c("p1", "p2", "p3", "p4"),
#'                                model = mod)
#' suggest_name_plot(data = plot_data, prop = c("p1", "p2", "p3", "p4"))
#'
#' ## Show specific curves
#' suggest_name_plot(data = plot_data[plot_data$.Group %in% c(1, 4), ],
#'                   prop = c("p1", "p2", "p3", "p4"))
#'
#' ## Show uncertainty using `se = TRUE`
#' suggest_name_plot(data = plot_data[plot_data$.Group %in% c(1, 4), ],
#'                   prop = c("p1", "p2", "p3", "p4"), se = TRUE)
#'
#' ## Change colours using `colours`
#' suggest_name_plot(data = plot_data[plot_data$.Group %in% c(1, 4), ],
#'                   prop = c("p1", "p2", "p3", "p4"), se = TRUE,
#'                   colours = c("steelblue1", "steelblue4", "orange1", "orange4"))
#'
#' ## Facet plot based on specific variables
#' suggest_name_plot(data = plot_data,
#'                   prop = c("p1", "p2", "p3", "p4"), se = TRUE,
#'                   facet_var = "block",
#'                   colours = c("steelblue1", "steelblue4", "orange1", "orange4"))
#'
#' ## Simulataneously create multiple plots for additional variables
#' sim2$block <- as.numeric(sim2$block)
#' new_mod <- update(mod, ~ . + block, data = sim2)
#' plot_data <- suggest_name_data(starts = sim2[c(18), 3:6],
#'                        ends = sim2[c(48, 60), 3:6],
#'                        prop = c("p1", "p2", "p3", "p4"),
#'                        model = new_mod, conf.level = 0.95,
#'                        add_var = list("block" = c(1, 2)))
#'
#' suggest_name_plot(data = plot_data,
#'                   prop = c("p1", "p2", "p3", "p4"), se = TRUE,
#'                   colours = c("steelblue1", "steelblue4", "orange1", "orange4"))
suggest_name_plot <- function(data, prop, colours = NULL,
                              se = FALSE, facet_var = NULL,
                              nrow = 0, ncol = 0){
  # Sanity checks
  if(missing(data)){
    cli::cli_abort(c("{.var data} cannot be empty.",
                     "i" = "Specify a data-frame or tibble
                     preferably the output of {.fn suggest_name_data} or
                     a data-frame with a similar structure and column names."))
  }

  # Ensure identifiers for columns in data giving species proportions are specified
  if(missing(prop)){
    cli::cli_abort(c("{.var prop} cannot be empty.",
                     "i" = "Specify either a character or numeric vector giving
                     names/indicies of the columns containing the
                     compositional variables in {.var data}."))
  }

  sanity_checks(data = data, prop = prop,
                colours = colours,
                booleans = list("se" = se),
                numerics = list("nrow" = nrow, "ncol" = ncol),
                unit_lengths = list("nrow" = nrow, "ncol" = ncol))

  if(check_col_exists(data, ".add_str_ID")){
    ids <- unique(data$.add_str_ID)
    plots <- lapply(cli_progress_along(1:length(ids), name = "Creating plot",
                                       format = paste0(
                                         "{pb_spin} Creating plot ",
                                         "[{pb_current}/{pb_total}]   ETA:{pb_eta}"
                                       )),
                    function(i){
                      data_iter <- data %>% filter(.data$.add_str_ID == ids[i])
                      suggest_name_plot_internal(data = data_iter, prop = prop,
                                                 colours = colours, se = se,
                                                 facet_var = facet_var)+
                        labs(subtitle = ids[i])
                    })
    if(length(plots) > 1){
      plot <- new("ggmultiplot", plots = plots, nrow = nrow, ncol = ncol)
    } else {
      plot <- plots[[1]]
    }
    cli::cli_alert_success("Created all plots.")
  } else {
    plot <- suggest_name_plot_internal(data = data, prop = prop,
                              colours = colours, se = se,
                              facet_var = facet_var)
    cli::cli_alert_success("Created plot.")
  }
  plot
}

#' @title
#' Visualising the change in a response variable between two points in
#' the simplex space
#'
#' @description
#' This function will prepare the underlying data and plot the results for visualising the
#' change in a response variable between two points in the simplex space in a single
#' function call. The two points are specified by the `starts` and `ends`
#' parameters are joined by a straight line across the simplex space and the
#' response is predicted for the starting, ending and intermediate communities
#' along this line. The associated uncertainty along this prediction is also shown.
#' The generated plot will show individual curves indicating the variation in the
#' response between the points. `\code{\link[PieGlyph:PieGlyph-package]{Pie-glyphs}}`
#' are used to highlight the compositions of the starting, ending and midpoint
#' of the straight line between the two points.
#' This is a wrapper function specifically for statistical models fit using the
#' \code{\link[DImodels:DI]{DI()}} function from the
#' \code{\link[DImodels:DImodels-package]{DImodels}} R package and would implicitly
#' call \code{\link{suggest_name_data}} followed by
#' \code{\link{suggest_name_plot}}. If your model object isn't fit using
#' DImodels, consider calling these functions manually.
#'
#'
#' @inheritParams visualise_effects
#' @inheritParams suggest_name_data
#' @inheritParams suggest_name_plot
#' @inheritParams add_prediction
#'
#' @inherit prediction_contributions return
#'
#' @export
#'
#' @examples
#' library(DImodels)
#' data(sim2)
#'
#' # Fit model
#' mod <- DI(y = "response", prop = 3:6, DImodel = "AV", data = sim2)
#'
#' # Create plot
#' # Move from p3 monoculture to p4 monoculture
#' suggest_name(model = mod,
#'              starts = data.frame(p1 = 0, p2 = 0, p3 = 1, p4 = 0),
#'              ends = data.frame(p1 = 0, p2 = 0, p3 = 0, p4 = 1))
#'
#' # Move from 50:50 mixture of p1 and p2 towards their monoculture
#' suggest_name(model = mod,
#'              starts = data.frame(p1 = 0.5, p2 = 0.5, p3 = 0, p4 = 0),
#'              ends = data.frame(p1 = c(1, 0), p2 = c(0, 1), p3 = 0, p4 = 0))
#'
#' # Move from dominant mixtures to p1 monoculture
#' suggest_name(model = mod,
#'              starts = sim2[c(1, 5, 9, 13), 3:6],
#'              ends = data.frame(p1 = 1, p2 = 0, p3 = 0, p4 = 0))
#'
#' # Move from centroid community to each monoculture
#' suggest_name(model = mod,
#'              starts = sim2[c(18),],
#'              ends = sim2[c(48, 52, 56, 60), ])
#'
#' # Show change across multiple points simulataneously
#' suggest_name(model = mod,
#'              starts = sim2[c(1, 17, 22), ],
#'              ends = sim2[c(5, 14, 17), ])
#'
#' # Show confidence bands
#' suggest_name(model = mod,
#'              starts = sim2[c(1, 17, 22), ],
#'              ends = sim2[c(5, 14, 17), ], se = TRUE)
#'
#' # Change colours
#' suggest_name(model = mod,
#'              starts = sim2[c(1, 17, 22), ],
#'              ends = sim2[c(5, 14, 17), ], se = TRUE,
#'              colours = c("steelblue1", "steelblue4", "orange1", "orange4"))
#'
#' # Facet based on existing variables
#' suggest_name(model = mod,
#'              starts = sim2[c(1, 17, 22), ],
#'              ends = sim2[c(5, 14, 17), ], se = TRUE, facet_var = "block",
#'              colours = c("steelblue1", "steelblue4", "orange1", "orange4"))
#'
#' # Add additional variables and create a separate plot for each
#' suggest_name(model = mod,
#'              starts = sim2[c(1, 17, 22), 3:6],
#'              ends = sim2[c(5, 14, 17), 3:6], se = TRUE,
#'              colours = c("steelblue1", "steelblue4", "orange1", "orange4"),
#'              add_var = list("block" = factor(c(1, 3),
#'                                              levels = c(1, 2, 3, 4))))
#'
#' ## Specify `plot = FALSE` to not create the plot but return the prepared data
#' suggest_name(model = mod, plot = FALSE,
#'              starts = sim2[c(1, 17, 22), 3:6],
#'              ends = sim2[c(5, 14, 17), 3:6], se = TRUE,
#'              colours = c("steelblue1", "steelblue4", "orange1", "orange4"),
#'              add_var = list("block" = factor(c(1, 3),
#'                                              levels = c(1, 2, 3, 4))))
suggest_name <- function(model, starts, ends, add_var = list(),
                         interval = c("confidence", "prediction", "none"),
                         conf.level = 0.95,
                         se = FALSE, average = TRUE,
                         colours = NULL, FG = NULL,
                         facet_var = NULL, plot = TRUE,
                         nrow = 0, ncol = 0){
  # Sanity checks
  # Ensure model is a DImodels object
  # Ensure specified model is fit using the DI function
  if(missing(model) || !inherits(model, "DI")){
    model_not_DI(call_fn = "visualise_effects")
  }

  if(missing(starts)){
    cli::cli_abort(c("{.var starts} cannot be empty.",
                     "i" = "Specify a data frame or tibble indicating the initial variable proportions in {.var data}."))
  }

  if(missing(ends)){
    cli::cli_abort(c("{.var ends} cannot be empty.",
                     "i" = "Specify a data frame or tibble indicating the final variable proportions in {.var data}."))
  }

  # Get original data used to fit the model
  original_data <- model$original_data

  # Get all species in the model
  model_species <- eval(model$DIcall$prop)
  # If species were specified as column indices extract names
  if(is.numeric(model_species)){
    model_species <- colnames(original_data)[model_species]
  }

  interval <- match.arg(interval)

  plot_data <- suggest_name_data(model = model, starts = starts, ends = ends,
                                 prop = model_species, add_var = add_var,
                                 interval = interval, conf.level = conf.level)

  # Get functional groups
  if(is.null(FG)){
    FG <- eval(model$DIcall$FG)
  }

  # Colours for species
  if(is.null(colours)){
    colours <- get_colours(vars = model_species, FG = FG)
  }

  if(isTRUE(plot)){
    plot <- suggest_name_plot(data = plot_data, prop = model_species,
                              colours = colours, se = se,
                              facet_var = facet_var,
                              nrow = 0, ncol = 0)
    return(plot)
  } else {
    return(plot_data)
  }
}


#' @keywords internal
#' Internal function for creating a plot showing the change in response between
#' any two points in the simplex
#'
#' @usage NULL
NULL
suggest_name_plot_internal <- function(data, prop, colours = NULL,
                                       se = FALSE, facet_var = NULL){

  # Get names of columns containing species proportions
  species_names <- data %>% select(all_of(prop)) %>% colnames()

  pie_data <- data %>%
    group_by(.data$.Group) %>%
    slice(1, 51, 101) %>%
    ungroup() %>%
    # Filter out any overlapping pies to avoid overplotting
    distinct(.data$.InterpConst, .data$.Pred, .keep_all = T) %>%
    ungroup()


  # Colours for the pie-glyph slices
  if(is.null(colours)){
    colours <- get_colours(species_names)
  }

  # Create canvas for plot
  plot <- ggplot(data, aes(x = .data$.InterpConst, y = .data$.Pred))+
    theme_bw()

  # Add ribbons for uncertainty of prediction
  if(se){
    plot <- plot +
      geom_ribbon(aes(ymin = .data$.Lower, ymax = .data$.Upper,
                      group = .data$.Group),
                  colour = 'grey', alpha = 0.25)
  }

  # Add line tracing the effect of adding a particular species to the data
  plot <- plot +
    geom_line(aes(group = .data$.Group), colour = 'black', alpha = 0.75)

  # Add the pie-chart glyphs for identifying the data
  plot <- plot +
    geom_pie_glyph(data = pie_data, radius = 0.3,
                   slices = prop, colour = 'black')+
    scale_fill_manual(values = colours,
                      labels = prop)

  # Add facet if specified
  if(!is.null(facet_var)){
    plot <- add_facet(plot, data, facet_var,
                      labeller = label_both)
  }

  # Adjust plot aesthetics
  plot <- plot +
    labs(fill = "Variable",
         x = "Interpolation constant",
         y = "Predicted Response",
         caption = "The pie-glyphs on the left show the starting compositions
         while those on the right show the ending compositions.
         The pie-glyphs in the centre show the composition of the
         point midway between the starting and ending compositions.")+
    theme(legend.position = 'top')

  return(plot)
}

#' @keywords internal
#' Utility function that accepts any two points in the simplex space as
#' two numeric vectors and returns a data-frame containing `ncomms` number of
#' intermediate communities across a straight line between the two points.
#'
#' @usage NULL
NULL
interpolate_communities <- function(start, end, prop, ncomms = 101){
  ts <- seq(0, 1, length.out = ncomms)
  points <- sapply(ts, function(t) {
    start + t * (end - start)
  }) %>% t()
  points <- apply(points, 2, unlist) %>%
    as.data.frame() %>%
    `colnames<-`(prop)
  return(points)
}
