#' @title Conditional ternary diagrams
#'
#' @description
#' The helper function for preparing the underlying data for creating conditional
#' ternary diagrams, where the high dimensional simplex is sliced at various
#' values along the range of a particular variable and we visualise the change
#' in the response with respect to the remaining three variables in a ternary
#' diagram where the proportions within the ternary would sum to 1 - x, where
#' x is the value of the conditioning variable at which the simplex is sliced.
#'
#' @param prop A character vector indicating the model coefficients
#'             corresponding to variable proportions.
#' @param conditional A character vector describing the variable(s) in the
#'                    direction of which to condition the high dimensional
#'                    simplex for taking the slices.
#' @param values A vector of numbers between 0 and 1 describing the values
#'               at which to takes slices of the high dimensional simplex.
#' @inheritParams ternary_data
#' @inheritDotParams add_prediction -data
#'
#' @return A data-frame containing compositional columns with names specified
#'         in `prop` parameter along with any additional columns specified in
#'         `exp_str` parameter and the following columns appended at the end.
#'  \describe{
#'    \item{.x}{The x-projection of the points within the ternary.}
#'    \item{.y}{The y-projection of the points within the ternary.}
#'    \item{.add_str_ID}{An identifier column for grouping the cartesian product
#'                       of all additional columns specified in `exp_str`
#'                       parameter (if `exp_str` is specified).}
#'    \item{.Sp}{An identifier column specifying the variable along which the
#'               high dimensional simplex is sliced.}
#'    \item{.Value}{The value (between 0 and 1) along the direction of variable
#'                  in `.Sp` at which the high dimensional simplex is sliced.}
#'    \item{.Facet}{An identifier column formed by combining `.Sp` and `.value`
#'                  to group observations within a specific slice of the
#'                  high dimensional simplex.}
#'    \item{.Pred}{The predicted response for each community.}
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
#'
#' ## Load data
#' data(sim2)
#'
#' ## Fit model
#' mod <- glm(response ~ 0 + (p1 + p2 + p3 + p4)^2, data = sim2)
#'
#' ## Create data
#' head(conditional_ternary_data(prop = c("p1", "p2", "p3", "p4"),
#'                               conditional = "p4",
#'                               model = mod,
#'                               resolution = 1))
#'
#' ## Can also add any additional experimental structures
#' head(conditional_ternary_data(prop = c("p1", "p2", "p3", "p4"),
#'                               conditional = "p4",
#'                               exp_str = list("block" = c("1", "2")),
#'                               model = mod,
#'                               resolution = 1))
#'
#' ## It could be desirable to take the output of this function and add
#' ## additional variables to the data before making predictions
#' ## Use `prediction = FALSE` to get data without any predictions
#' head(conditional_ternary_data(prop = c("p1", "p2", "p3", "p4"),
#'                               conditional = "p4",
#'                               prediction = FALSE,
#'                               resolution = 1))
conditional_ternary_data <- function(prop,
                                     conditional = NULL, values = c(0.2,0.5,0.8),
                                     exp_str = list(),
                                     resolution = 3, prediction = TRUE, ...){

  #Sanity Checks
  if(missing(prop)){
    cli::cli_abort(c("{.var prop} cannot be empty.",
                     "i" = "Specify a character vector indicating the model
                            coefficients corresponding to variable
                            proportions."))
  }

  if(length(prop) < 3){
    cli::cli_abort(c("Ternary diagrams can only be created for models with more
                      than or equal to 3 species.",
                     "i" = "Currently only {length(prop)} are specified in
                            {.var prop}."))
  }

  # Don't condition on anything if not specified
  if(is.null(conditional)){
    conditional <- ''
  }

  # Check conditional parameter
  if (!inherits(conditional, 'character')){
    cli::cli_abort("{.var conditional} should contain the names of the
                   conditioning variable as a string vector",
                   "i" = "{.var conditional} was specified as a
                          {.cls {conditional}}.")
  }

  if (all(conditional !='') & !all(conditional %in% prop)){
    cli::cli_abort(c("The value specified in {.var conditional} should be
                      present in {.var prop}.",
                     "i" = "{conditional[!conditional %in% prop]} {?is/are}
                            not specified in {.var prop}."))
  }

  # Check values parameter
  if(!all(is.numeric(values))){
    cli::cli_abort(c("{.var values} should be a numeric vector with values
                     between 0 and 1 describing the values at which to condition
                     the n-dimensional simplex space.",
                     "i" = "{.var values} was specified as a {.cls {values}}
                            object."))
  }

  if(!all(between(values, 0 , 1))){
    cli::cli_abort(c("{.var values} should be a numeric vector with values
                     between 0 and 1 describing the values at which to condition
                     the n-dimensional simplex space.",
                     "i" = "The value{?s} specified in {.var values} was
                            {.val {values}}."))
  }

  # Get species to be shown in the ternary diagram
  tern_species <- prop[prop!=conditional]

  if (length(tern_species)!=3){
    if (length(tern_species)<3){
      cli::cli_abort(c("Ternary diagram can not be created with less than 3
                       species.",
                       "i" = "After accounting for the values specified in
                             {.var conditional}, only {length(tern_species)}
                             are left for creating the ternary diagram.",
                       "i" = "Remove any {3 - length(tern_species)} species
                             from the {.var conditional} parameter."))
    } else if (length(tern_species>3)){
      cli::cli_warn(c("After accounting for the values specified in
                      {.var conditional}, there are more than three species
                      for ternary diagram.",
                      "i" = "Creating the diagram for species
                             {.val {tern_species[1:3]}} and assuming
                             {.val {tern_species[4:length(tern_species)]}}
                             to be 0, but it might not be very informative.",
                      "i" = "Try increasing the number of values in
                             {.var conditional} or grouping species into
                             functional groups
                             (see {.var ?DImodelsVis::functional_ternary})."))
    }
  }

  # Create base ternary data to be plotted
  triangle <- ternary_data(prop = tern_species[1:3],
                           exp_str = exp_str,
                           resolution = resolution,
                           prediction = FALSE)

  # Add if there are more than three species to show in the ternary
  # add the remaining species which are assumed to have a proportion of 0
  if (length(tern_species)>3){
    other_species <- tern_species[4:length(tern_species)]
    pOther <- matrix(0, ncol=length(other_species), nrow= nrow(triangle))
    colnames(pOther) <- other_species
    triangle <- cbind(triangle, pOther)
  }

  # If there are no species to condition on then i.e. conditional = ''
  # our data is ready and we can make predictions from the model
  if (all(conditional=='')){
    if(prediction){
      cond_data <- add_prediction(data = triangle, ...)
    }
    # cond_data <- add_prediction(data = triangle,
    #                             model = model,
    #                             coefficients = coefficients,
    #                             coeff_cols = coeff_cols,
    #                             pred_name = pred_name,
    #                             interval = "none")
  } else {
    # Rescaling species to be shown in ternary for each species
    # in the conditional parameter according to the values specified
    # in the values parameter
    cond_data <- lapply(cli_progress_along(values, name = "Preparing data"), function(idx){
      x <- values[idx]
      # No need to scale if x is zero
      if(x == 0){
        scaled_data <- triangle
      } else {
        # Scale proportion of species within the ternary
        scaled_data <- triangle %>%
          mutate(!! tern_species[1] := rescale(!!sym(tern_species[1]),
                                               min = 0, max = 1-x),
                 !! tern_species[2] := rescale(!!sym(tern_species[2]),
                                               min = 0, max = 1-x),
                 !! tern_species[3] := rescale(!!sym(tern_species[3]),
                                               min = 0, max = 1-x))
      }
      # Add remaining species in the data
      remaining_species <- matrix(0, ncol=length(conditional),
                                  nrow= nrow(triangle))
      colnames(remaining_species) <- conditional
      scaled_data <- cbind(scaled_data, remaining_species)

      # Update value of species which we are conditioning on
      species_data <- lapply(conditional, function(cond_sp){
        # Add identifier for grouping data for a particular species
        sp_data <- scaled_data %>%
                      mutate(.Sp = cond_sp,
                             .Value = x,
                             .Facet = paste0(cond_sp, ' = ', x))

        # To avoid any rounding issues & ensure all species proportions sum to 1
        sp_data[, cond_sp] <- 1 - rowSums(sp_data[, tern_species])
        # Predicting the response for the communities
        if(prediction){
          sp_data <- add_prediction(data = sp_data, ...)
        }
        # Need to return data via subset of rows as bind_rows fails otherwise
        sp_data
        }) %>% bind_rows()
    }) %>% bind_rows()
  }
  # Final formatting to pretty up data
  selection <- if(all(conditional == "")) tern_species else c(tern_species, conditional)
  cond_data <- cond_data %>%
    select(all_of(c(selection, names(exp_str))), everything())
  cli::cli_alert_success("Finished data preparation.")
  return(cond_data)
}

#' @title Conditional ternary diagrams
#'
#' @description
#' The helper function for plotting conditional ternary diagrams. The output of
#' the `\code{\link{conditional_ternary_data}}` should be passed here to
#' visualise the n-dimensional simplex space as 2-d slices showing the change
#' in the response across any three variables, when a fourth is fixed at a
#' particular value.
#'
#' @importFrom ggplot2 facet_wrap
#'
#' @param data A data-frame which is the output of the
#'             `\link{conditional_ternary_data}` function.
#' @inheritParams ternary_plot
#'
#' @inherit ternary_plot return
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
#' mod <- glm(response ~ 0 + (p1 + p2 + p3 + p4)^2, data = sim2)
#'
#' ## Create data for slicing
#' plot_data <- conditional_ternary_data(prop = c("p1", "p2", "p3", "p4"),
#'                                       conditional = "p4",
#'                                       model = mod,
#'                                       resolution = 1)
#'
#' ## Create plot
#' conditional_ternary_plot(data = plot_data)
conditional_ternary_plot <- function(data,
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
  if(missing(data)){
    cli::cli_abort(c("{.var data} cannot be missing.",
                    "i" = "Specify the output of
                          {.fn conditional_ternary_data} function."))
  }

  # Create the simple ternary plot
  pl <- ternary_plot(data = data,
                     nlevels = nlevels,
                     colours = colours,
                     tern_labels = tern_labels,
                     lower_lim = lower_lim,
                     upper_lim = upper_lim,
                     contour_text = contour_text,
                     show_axis_labels = FALSE,
                     show_axis_guides = show_axis_guides,
                     axis_label_size = axis_label_size,
                     vertex_label_size = vertex_label_size)

  # Check if we need to condition on anything and get appropriate values
  if(!is.null(data$.Facet)){
    conditional <- unique(data$.Sp)
    values <- unique(data$.Value)

    # Facet and create one panel for each species to condition on
    pl <- pl +
      facet_wrap(~ .Facet, ncol = length(values))

  } else {
    conditional <- ""
    values <- 0
  }

  # Show appropriate labels for for the ternary axes
  if(show_axis_labels){
    # Labels for the ternary axes
    # (because they'll be scaled for conditional panels)
    axis_labels <- lapply(conditional, function(sp){
      lapply(values, function(x){
        positions <- tibble(x1 = seq(0.2,0.8,0.2),
                            y1 = c(0,0,0,0),
                            x2 = .data$x1/2,
                            y2 = .data$x1*sqrt(3)/2,
                            x3 = (1-.data$x1)*0.5+.data$x1,
                            y3 = sqrt(3)/2-.data$x1*sqrt(3)/2,
                            label = .data$x1*(1-x),
                            rev_label = rev(.data$label),
                            .Facet = paste0(sp, ' = ', x),
                            .Pred = 0)
      })
    }) %>% bind_rows()

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

  return(pl)
}


#' @title Conditional ternary diagrams
#'
#' @description
#' Conditional ternary diagrams are a way to visualise n-dimensional
#' compositional data residing in the n-1 dimensional space as 2-d ternary
#' diagrams. We slice the high dimensional simplex at various values along
#' the range of a particular  variable and visualise the change in the
#' response with respect to the remaining three variables in a ternary
#' diagram where the proportions within the ternary would sum to 1 - x,
#' where x is the value of the conditioning variable at which the simplex is
#' sliced. Taking multiple 2-d slices across multiple variables should allow to
#' create an approximation of how the response varies across the n-dimensional
#' simplex.
#'
#'
#' @importFrom metR geom_text_contour
#' @importFrom grDevices terrain.colors
#' @importFrom dplyr tibble between
#' @importFrom ggplot2 geom_raster scale_fill_stepsn geom_contour
#'                     theme_void guide_colorbar coord_fixed
#'
#' @param model A Diversity Interactions model object fit by using the
#'              \code{\link[DImodels:DI]{DI()}} function from the
#'              \code{\link[DImodels:DImodels-package]{DImodels}} package.
#' @inheritParams conditional_ternary_data
#' @inheritParams ternary_plot
#' @inheritParams model_diagnostics
#'
#' @inherit response_contributions return
#'
#' @export
#'
#' @examples
#' library(DImodels)
#' data(sim2)
#' m1 <- DI(y = "response", data = sim2, prop = 3:6, DImodel = "FULL")
#'
#' conditional_ternary(m1, conditional = "p4", resolution = 1)
#'
#' ## Slices for experiments for over 4 species
#' data(sim4)
#' m2 <- DI(y = "response", prop = paste0("p", 1:6),
#'          DImodel = "AV", data = sim4)
#'
#' conditional_ternary(m2, conditional = c("p4", "p5", "p6"),
#'                     resolution = 1, values = c(0.2, 0.5))
conditional_ternary <- function(model, conditional = NULL,
                                values = c(0.2, 0.5, 0.8),
                                exp_str = list(),
                                resolution = 3,
                                nlevels = 7,
                                colours = NULL,
                                lower_lim = NULL,
                                upper_lim = NULL,
                                contour_text = TRUE,
                                show_axis_labels = TRUE,
                                show_axis_guides = FALSE,
                                axis_label_size = 4,
                                vertex_label_size = 5,
                                nrow = 0, ncol = 0){

  # Ensure specified model is fit using the DI function
  if(missing(model) || !inherits(model, "DI")){
    model_not_DI(call_fn = "conditional_ternary")
  }

  # Get data used to fit the model
  og_data <- model$original_data

  # Get all species in the model
  species <- eval(model$DIcall$prop)

  if(is.numeric(species)){
    species <- colnames(og_data)[species]
  }

  # Create data in appropriate format for plotting
  plot_data <- conditional_ternary_data(prop = species,
                                        model = model,
                                        conditional = conditional,
                                        values = values,
                                        exp_str = exp_str,
                                        resolution = resolution)

  # Labels for the ternary
  tern_labels <- colnames(plot_data)[1:3]

  if(length(exp_str) > 0){
    ids <- unique(plot_data$.add_str_ID)
    plots <- lapply(cli_progress_along(1:length(ids), name = "Creating plot",
                                       format = paste0(
                                         "{pb_spin} Creating plot ",
                                         "[{pb_current}/{pb_total}]   ETA:{pb_eta}"
                                       )),
                    function(i){
                      data <- plot_data %>% filter(.data$.add_str_ID == ids[i])
                      plot <- conditional_ternary_plot(data = data,
                                                       nlevels = nlevels,
                                                       colours = colours,
                                                       tern_labels = tern_labels,
                                                       lower_lim = lower_lim,
                                                       upper_lim = upper_lim,
                                                       contour_text = contour_text,
                                                       show_axis_labels = show_axis_labels,
                                                       show_axis_guides = show_axis_guides,
                                                       axis_label_size = axis_label_size,
                                                       vertex_label_size = vertex_label_size) +
                        labs(subtitle = ids[i])
                    })
    if(length(plots) > 1){
      plot <- new("ggmultiplot", plots = plots, nrow = nrow, ncol = ncol)
    } else {
      plot <- plots[[1]]
    }
    cli::cli_alert_success("Created all plots.")
  } else {
    plot <- conditional_ternary_plot(data = plot_data,
                                     nlevels = nlevels,
                                     colours = colours,
                                     tern_labels = tern_labels,
                                     lower_lim = lower_lim,
                                     upper_lim = upper_lim,
                                     contour_text = contour_text,
                                     show_axis_labels = show_axis_labels,
                                     show_axis_guides = show_axis_guides,
                                     axis_label_size = axis_label_size,
                                     vertex_label_size = vertex_label_size)
    cli::cli_alert_success("Created plot.")
  }
  return(plot)
}
