#' @title Create data for visualising the effect of species loss or addition from a community.
#'
#' @description
#' The helper function to create the underlying data for visualising effects plots.
#' For each species community in data, several intermediary communities would be created between the initial
#' community and the monoculture of the species of interest (if looking at effect of species addition) or
#' between the original community and the community without the species of interest (if look at effect of species loss)
#' whilst maintaining the same ratio of the original species in the community.
#' The model object or coefficients vector would then be used to calculated the predicted response for each
#' intermediary community allowing us to visualise the effect of the addition/loss of the species of interest.
#' This function is particularly useful if the model object is not of class `\link[DImodels:DI]{DI()}` or only model
#' coefficients are available for prediction.
#'
#' @param data A dataframe specifying the initial communities of interest for which
#'             to visualise the effect of species loss or addition.
#'             If a model object is specified then this data should contain all the
#'             variables present in the model object including any experimental structures.
#'             If a coefficient vector is specified then data should contain same number of
#'             columns as the number of elements in the coefficient vector and a one-to-one
#'             positional mapping would be assumed between the data columns and the
#'             elements of the coefficient vector.
#' @param prop A vector of column names or indices identifying the columns containing the
#'             species proportions in the data.
#' @param species_interest A character vector specifying the species for which to visualise
#'                         the effect of addition or loss on the response. If left blank,
#'                         all species would be assumed to be of interest.
#' @param effect A character string with one value from c("addition", "loss") whether to
#'               visualise the effect of addition or loss of a species.
#' @param prediction A logical value indicating whether to pass the final data to
#'                   `\link{add_prediction}` and add predictions to the data.
#'                   Default value is \code{TRUE}, but often it would be desirable
#'                   to make additional changes to the data before making any
#'                   predictions, so the user can set this to \code{FALSE} and
#'                   manually call the `\link{add_prediction}` function.
#' @inheritDotParams add_prediction -data
#'
#' @return A data frame with the following columns appended at the end
#'  \describe{
#'    \item{.Sp}{An identifier column to discern the species of interest being modified in each curve.}
#'    \item{.Proportion}{The value of the species of interest within the community.}
#'    \item{.Group}{An identifier column to discern between the different curves.}
#'    \item{.Pred}{The predicted response for each community.}
#'    \item{.Lower}{The lower limit of the prediction/confidence interval for each observation.}
#'    \item{.Upper}{The upper limit of the prediction/confidence interval for each observation.}
#'    \item{.Marginal}{The marginal change in the response (first derivative) with respect to the
#'                     gradual change in the proportion of the species of interest.}
#'    \item{.Threshold}{A numeric value indicating the maximum proportion of the species of interest
#'                      within a particular community which results in a positive marginal effect on the response.}
#'    \item{.MarEffect}{A character string entailing whether the addition/loss of the species of interest from the
#'                      particular community would result in a positive or negative marginal effect on the response.}
#'    \item{.Effect}{An identifier column signifying whether considering the effect of species addition or species loss.}
#'  }
#'
#' @export
#'
#' @examples
#' library(DImodels)
#'
#' ## Load data
#' data(sim1)
#'
#' ## Fit model
#' mod <- glm(response ~ p1 + p2 + p3 + p4 + 0, data = sim1)
#'
#' ## Create data for visualising effect of adding species 1 to
#' ## the original communities in the data
#' head(visualise_effects_data(data = sim1, prop = c("p1", "p2", "p3", "p4"),
#'                             species_interest = "p1", effect = "addition",
#'                             model = mod))
#'
#' ## Create data for visualising effect of losing species 2 from
#' ## the original communities in the data but using model coefficients
#' ## When specifying coefficients the data should have a one-to-one
#' ## positional mapping with specified coefficients.
#' ## Note that no confidence intervals can be generated when specifying
#' ## model coefficients.
#' init_comms <- sim1[, c("p1", "p2", "p3", "p4")]
#' head(visualise_effects_data(data = init_comms, prop = 1:4,
#'                             species_interest = "p2",
#'                             effect = "loss",
#'                             coefficients = mod$coefficients))
#'
#' ## Can also create only the intermediary communities without predictions
#' ## by specifying prediction = FALSE.
#' ## Any additional columns can then be added and the `add_prediction` function
#' ## can be manually called.
#' ## Note: If calling the `add_prediction` function manually, the data would
#' ## not contain information about the marginal effect of changing the species
#' ## interest
#' effects_data <- visualise_effects_data(data = init_comms, prop = 1:4,
#'                                        species_interest = "p2",
#'                                        effect = "loss",
#'                                        prediction = FALSE)
#' head(effects_data)
#' head(add_prediction(data = effects_data, model = mod, interval = "prediction"))
#' head(add_prediction(data = effects_data, coefficients = mod$coefficients))
visualise_effects_data <- function(data, prop, species_interest = NULL,
                                   effect = c("addition", "loss"),
                                   prediction = TRUE, ...){
  if(missing(data)){
    cli::cli_abort(c("{.var data} cannot be empty.",
                     "i" = "Specify a data frame or tibble indicating the species communities
                     for which to calculate the effect of addition/loss of a species in {.var data}."))
  }

  #check_data_functions(model = model, coefficients = coefficients, species = species)

  # Get names of columns containing species proportions
  species_names <- data %>% select(all_of(prop)) %>% colnames()

  # If species_interest is not specified assume all species are of interest
  if(is.null(species_interest)){
      cli::cli_warn(c("{.var species_interest} was not specified.",
                    "i" = "Assuming all species are of interest."))
    species_interest <- species_names
  }

  # Ensure effect is only one of 'addition' or 'loss'
  effect <- match.arg(effect)

  # Progress bar
  #pro_bar <- cli_progress_bar("Preparing data", total = 100, clear = FALSE)

  # Create data to create effects plot
  plot_data <- lapply(cli_progress_along(species_interest, name = "Preparing data"), function(idx){
    x <- species_interest[idx]
    if(effect == 'addition'){
      # Interval or calculating effect of adding species
      pvals <- seq(0, 1, length.out = 101)
      # Choose all species other than the species of interest to perform extrapolation to see effect of adding species of interest
      pothers <- data %>% select(-all_of(x))
      # Filter out the monoculture of the species of interest as can't add species to that
      pothers <- pothers[rowSums(pothers[, species_names[species_names!=x]])!=0, ]
      # Rescale species proportions to ensure they sum to 1
      pothers[, species_names[species_names!=x]] <- pothers[, species_names[species_names!=x]]/rowSums(pothers[, species_names[species_names!=x]])
      pother_names <- colnames(pothers[, species_names[species_names!=x]])
      # Filter duplicate data after the rescaling
      pothers <- unique(pothers)
    } else {
      # If we are looking at species loss then we don't need to look at monocultures of species of interest or in mixture which don't contain
      # our species of interest
      pothers <- data %>% filter(!(!!sym(x) %in% c(0, 1)))
      pother_names <- colnames(pothers[, species_names[species_names!=x]])
      # Filter any duplicate data
      pothers <- unique(pothers)
    }

    # For each row in data create copies with species of interest having one value from pvals (0.01, 0.02, 0.03, ..., 0.98, 0.99, 1)
    data_expand <- lapply(seq_len(nrow(pothers)), function(i){
      # Intervals for calculating the effects of species loss
      if(effect == 'loss'){
        pvals <- seq(pothers[i, x], 0, length.out = 101)
      }
      pothers_expand <- pothers %>%
        slice(i) %>%
        slice(rep(1, each = length(pvals))) %>%
        mutate(!!sym(x) := pvals,
               ".Sp" = x,
               ".Proportion" = pvals,
               ".Group" = i) %>%
        select(all_of(species_names), everything())

      # Rescale the proportions of the other species accordingly to ensure the proportions sum to 1
      pothers_expand[, pother_names] <- pothers_expand[, pother_names]*(1 - pvals)/rowSums(pothers_expand[, pother_names])

      # Make prediction and get marginal effect
      if(prediction){
        pothers_expand <- add_prediction(pothers_expand, interval = "confidence", ...)
        # Calculate the marginal effect of adding a species for producing the marginal plots
        dy <- diff(pothers_expand$.Pred)/diff(pvals)
        pothers_expand <- pothers_expand %>%
          mutate('.Marginal' = c(dy, dy[length(dy)]),
                 '.Threshold' = pvals[abs(.data$.Marginal) == min(abs(.data$.Marginal))][1],
                 '.MarEffect' = ifelse(!!sym(x) < .data$.Threshold, 'Negative', 'Positive'))
      }
      #cli_progress_update(id = pro_bar, force = TRUE)
      # Need to return data via subset of rows as bind_rows fails otherwise
      pothers_expand[1:nrow(pothers_expand),]

    }) %>% bind_rows()
    data_expand
  }) %>% bind_rows()
  #cli_progress_done(id = pro_bar)
  # Add flag to identify whether we are looking a effect of addition or effect of loss
  plot_data$.Effect <- effect
  cli::cli_alert_success("Finished data preparation.")
  return(plot_data)
}

#' @title Visualise the effects of species loss or species addition
#'
#' @description
#' The plotting function which would create the effects plot. The output of the
#' `\code{\link{visualise_effects_data}}` function should be passed here to
#' visualise the effect of adding or removing a particular species from a community.
#'
#' @param data A data frame created using the \code{\link{visualise_effects_data}} function.
#' @param prop A vector of column names or indices identifying the columns containing the
#'             species proportions in the data.
#' @param colours A character vector indicating the colours for the slices in the pie-glyphs.
#'                If left NULL, the colour blind friednly colours will be for the pie-glyph slices.
#' @param se A boolean variable indicating whether to plot confidence intervals associated with
#'           the effect of species addition or loss
#' @param average A boolean variable indicating whether to add a line describing the "average"
#'                effect of species addition or loss
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' library(DImodels)
#'
#' ## Load data
#' data(sim1)
#'
#' ## Fit model
#' mod <- glm(response ~ p1 + p2 + p3 + p4 + 0, data = sim1)
#'
#' ## Create data for visualising effect of adding species 1 to
#' ## the original communities in the data
#' plot_data <- visualise_effects_data(data = sim1,
#'                                     prop = c("p1", "p2", "p3", "p4"),
#'                                     species_interest = "p1",
#'                                     effect = "addition", model = mod)
#'
#' ## Create plot
#' visualise_effects_plot(data = plot_data,
#'                        prop = c("p1", "p2", "p3", "p4"))
#'
#' ## Show specific curves and add prediction intervals
#' subset <- plot_data[plot_data$.Group %in% c(1,7), ]
#' visualise_effects_plot(data = subset, prop = 1:4, se = TRUE)
#'
#' ## Do not show average effect line
#' visualise_effects_plot(data = subset, prop = 1:4,
#'                        se = TRUE, average = FALSE)
#'
#' ## Change colours of the pie-glyph slices
#' visualise_effects_plot(data = subset, prop = 1:4,
#'                        colours = c("darkolivegreen", "darkolivegreen1",
#'                                    "steelblue4", "steelblue1"))
visualise_effects_plot <- function(data, prop, colours = NULL,
                                   se = FALSE, average = TRUE){
  if(missing(data)){
    cli::cli_abort(c("{.var data} cannot be empty.",
                     "i" = "Specify a data frame or tibble (preferably the output of {.fn visualise_effects_data})."))
  }

  # Ensure identifiers for columns in data giving species proportions are specified
  if(missing(prop)){
    cli::cli_abort(c("{.var prop} cannot be empty.",
                     "i" = "Specify either a character or numeric vector indicating
                     the columns containg the species proportions in {.var prop}."))
  }

  # Get names of columns containing species proportions
  species_names <- data %>% select(all_of(prop)) %>% colnames()

  # Get data for showing the pie-chart glyphs
  pie_data <- data %>%
    group_by(.data$.Sp, .data$.Group) %>%
    filter(.data$.Proportion %in% c(min(.data$.Proportion), max(.data$.Proportion))) %>%
    ungroup() %>%
    group_by(.data$.Sp) %>%
    distinct(.data$.Pred, .keep_all = T)

  # Colours for the pie-glyph slices
  if(is.null(colours)){
    colours <- get_colours(species_names)
  }

  # Create canvas for plot
  plot <- ggplot(data, aes(x = .data$.Proportion, y= .data$.Pred))+
    theme_bw()

  # Add ribbons for uncertainty of prediction
  if(se){
    plot <- plot +
      geom_ribbon(aes(ymin = .data$.Lower, ymax = .data$.Upper,
                      group = .data$.Group),
                  colour = 'grey', linetype = 'dashed',
                  linewidth = 0.7, alpha = 0.25)
  }

  # Add line tracing the effect of adding a particular species to the data
  plot <- plot +
    geom_line(aes(group = .data$.Group), colour = 'black', alpha = 0.4)

  # If user specified to chose to see average effect of adding species then add that to plot
  if(average){
    # Round species proportion to 2 decimal places and calculate average effect.
    # Rounding done because for species loss because some proportions have recurring decimals
    avg_data <- data %>%
      mutate('.Proportion' = round(.data$.Proportion, 2))

    # Add line showing average effect
    plot <- plot +
      stat_summary(group = '.Sp', fun=mean, geom="line", colour="black",
                   linewidth = 1.25, data = avg_data)
  }

  # Add the pie-chart glyphs for identifying the data
  plot <- plot +
    geom_pie_glyph(aes(group = .data$.Sp), data = pie_data, radius = 0.3,
                   slices = prop, colour = 'black')+
    scale_fill_manual(values = colours,
                      labels = prop)

  # If looking at species loss then reverse the x-axis
  if(data$.Effect[1] == 'loss'){
    plot <- plot + scale_x_reverse()
  }
  # Finally facet plot by each species of interest to ensure better visibility
  plot <- plot +
    facet_wrap(~.data$.Sp)

  # Adjust plot aesthetics
  plot <- plot +
    labs(fill = "Species",
         x = "Proportion",
         y = "Predicted Response")+
    theme(legend.position = 'top')

  return(plot)
}

#' @title Visualise effect of species addition or loss
#' @description This function helps to visualise the effect of the addition or loss
#'              of a species from a community on the response
#'
#' @importFrom ggplot2 geom_line labs scale_colour_manual geom_ribbon
#'                     stat_summary scale_x_reverse facet_wrap stat_summary
#' @importFrom dplyr %>% select filter all_of everything slice mutate bind_rows
#' @importFrom PieGlyph geom_pie_glyph
#' @importFrom rlang :=
#'
#' @param model A Diversity Interactions model object fit by using the
#'              `\code{\link[DImodels:DI]{DI()}}` function from the
#'              `\code{\link[DImodels:DImodels-package]{DImodels}}` package.
#' @param data A dataframe specifying communities of interest for which user
#'             wants visualise the effect of species loss or addition.
#'             If left blank, the communities from the original data used
#'             to fit the model would be selected.
#' @param equi A boolean variable indicating whether to include the
#'             equi-proportional community in the plot.
#' @param design A boolean variable to indicate whether to include the
#'               communities in the design in the plot.
#' @inheritParams visualise_effects_data
#' @inheritParams visualise_effects_plot
#' @inheritParams response_contributions
#'
#' @return A \code{\link[ggplot2:ggplot]{ggplot2::ggplot}} ggplot or
#'         \code{\link[ggfortify:ggmultiplot-class]{ggfortify::ggmultiplot}}
#'         object.
#'
#' @export
#'
#' @examples
#' library(DImodels)
#'
#' ## Load data
#' data(sim1)
#'
#' ## Fit model
#' mod <- DI(prop = 3:6, DImodel = "AV", data = sim1, y = "response")
#'
#' ## Get effects plot for all species in design
#' visualise_effects(model = mod)
#'
#' ## Add custom communities to plot instead of design communities
#' visualise_effects(model = mod, design = FALSE,
#'                   data = data.frame("p1" = c(0.9, 0.1),
#'                                     "p2" = c(0.1, 0.9)))
#'
#' ## Add uncertainty on plot
#' visualise_effects(model = mod,
#'                   data = data.frame("p1" = c(0.9, 0.1),
#'                                     "p2" = c(0.1, 0.9)),
#'                   se= TRUE)
#'
#' ## Visualise effect of species loss for particular species
#' visualise_effects(model = mod, effect = "loss",
#'                   species_interest = c("p1", "p3"))
visualise_effects <- function(model, data = NULL, species_interest = NULL,
                              effect = c('addition', 'loss'),
                              exp_str = list(), equi = FALSE,
                              design = TRUE, se = FALSE, average = TRUE,
                              FG = NULL, conf.level = 0.95, nrow = 0, ncol = 0){

  # Sanity checks
  # Ensure model is a DImodels object
  if (!inherits(model, 'DI')){
    cli::cli_abort(c("{.var model} should be a object creating using the {.help [{.fun DI}](DImodels::DI)}
                     function from the {.help [{.pkg DImodels}](DImodels::DImodels)}
                     package or extending the {.cls DI} class.",
                     "x" = "You specified a model of class {.cls {class(model)}}."),
                   call = call)
  }

  # Get original data used to fit the model
  original_data <- model$original_data

  # Get all species in the model
  model_species <- eval(model$DIcall$prop)
  # If species were specified as column indices extract names
  if(is.numeric(model_species)){
    model_species <- colnames(original_data)[model_species]
  }

  # Failsafe to ensure user can't specify all three options to be false
  if(is.null(data) && (!is.null(equi) && !equi) && (!is.null(design) && !design)){
    stop("Specify a dataset of communities to show effects for or set either or both of 'equi' or 'design' to TRUE.")
  }

  # If the user has not specified neither of communities, equi or model_data. Then default to plotting the design communities
  if(is.null(data) & is.null(equi) & is.null(design)){
    design_points <- nrow(unique(original_data[, model_species]))
    # By default the model will plot all data in the original data. However this might not be feasible when there are a lot
    # of data in the design. If that's the case then plot only the equi-proportional mixtures
    if(design_points <= 100 ){
      design <- TRUE
    } else {
      equi <- TRUE
    }
  }

  # Create empty dataframe to append data to plot
  plot_data <- data.frame()

  # If user has manually specified the data then add those to data to be plotted
  if(!is.null(data)){
    # If user specifies data then ignore the design and equi parameters
    equi <- FALSE
    design <- FALSE
    # Ensure data are specified as a dataframe and in proper format
    if(!inherits(data,'data.frame')){
      stop("Specify data to show effects for in the plot as a data frame with proportions of the different species")
    } else {

      sp_abs <- !model_species %in% colnames(data)
      # If the column names of the data dataframe not match the species present in the
      # model then notify user asking them to change the column names
      if(all(sp_abs)){
        stop(glue::glue("The column names of the data data frame should be same as
                         the names of the species specified when fitting the model.
                         Update the column names to ensure they match
                        {paste(model_species, collapse = ', ')}"))
      }
      # If any species present in the model is not specified in the data, then notify user and assume its proportion to be 0
      if(any(sp_abs)){
        message(glue::glue("{paste(model_species[sp_abs], collapse = ', ')} weren't
                           present in the data specified, assuming their proportions to be 0."))
      }
      data[, model_species[sp_abs]] <- 0

      # If there are any additional columns in the data data.frame other than the species proportions. Notify user and ignore those
      extra_cols <- !colnames(data) %in% model_species
      if(any(extra_cols)){
        message(glue::glue("{paste(colnames(data)[extra_cols], collapse = ', ')}
                           were additional columns present in the data dataframe, ignoring them."))
      }
      # Append the user specified data to data dataframe
      plot_data <- rbind(plot_data, data[, model_species])
    }
  }

  # If plotting design data then add those to the data
  if(!is.null(design)){
    if(!is.logical(design)){
      stop(glue::glue("design takes boolean value of {TRUE} or {FALSE} specifying
                      whether or not to plot the design data"))
    } else {
      if(design){
        # Extract the design data from data used to fit the model
        design <- original_data[, model_species] %>% unique()

        # Append design data to data dataframe
        plot_data <- rbind(plot_data, design)
      }
    }
  }

  # If plotting the equi-proportional mixtures then add those to the data
  if(!is.null(equi)){
    if(!is.logical(equi)){
      stop(glue::glue('equi takes boolean value of {TRUE} or {FALSE} specifying whether or not to plot the design data'))
    } else {
      if(equi){
        # Create the equi-proportional mixture
        nSpecies <- length(model_species)
        prop <- 1/nSpecies
        equi <-  data.frame(matrix(prop, nrow = 1, ncol = nSpecies))
        colnames(equi) <- model_species

        # Append the equi proportional mixture to data dataframe
        plot_data <- rbind(plot_data, equi)
      }
    }
  }

  # Ensure the proportions of all data in the data sum approx to 1.
  if(!all(near(rowSums(plot_data), 1))){
    stop('Proportions of data in data should sum to 1')
  }

  # If user has not specified the species of interest then assume all species from model are of interest
  if(is.null(species_interest)){
    species_interest <- model_species
  }

  # If user has specified species_interest ensure it is a character vector
  if(!is.character(species_interest)){
    stop('Specify the species_interest parameter as a character vector containing names of the species of interest to show on plot')
  }

  # If user has specified the species of interest ensure they were present in the model
  if(!all(species_interest %in% model_species)){
    stop('species_interest should be present in the DImodels object')
  }

  # Ensure effect is only one of 'addition' or 'loss'
  effect <- match.arg(effect)

  # Get all combinations of experimental structures

  # No plots can be created if user is looking for effect of species loss but none of the data specified contain
  # the species of interest
  if(effect == 'loss'){
    if(all(plot_data[, species_interest] == 0)){
      stop("Specify data which contain some positive proportion of the species of
           interest to see the effect of its loss. Can't visualize effect of loss
           of species if no data contain the species of interest")
    }
  }

  # Similarly no plots can be created if the user is trying to find the effect of adding a species to it's monoculture
  if(effect == 'addition'){
    if(length(species_interest) == 1 & all(plot_data[, species_interest] == 1)){
      stop("Specify data other than the monoculture of the species of interest.
           Can't visualize effect of addition of species if only monocultures of
           species of interest are present.")
    }
  }

  # Get functional groups
  if(is.null(FG)){
    FG <- eval(model$DIcall$FG)
  }

  # Ensure experimental structure are specified correctly
  exp_str <- check_exp_str(model = model, exp_str = exp_str)

  if(length(exp_str) > 0){
    plot_data <- add_exp_str(exp_str = exp_str, data = plot_data)
  }


  plot_data <- visualise_effects_data(model = model, data = plot_data, effect = effect,
                                 prop = model_species, species_interest = species_interest)
  # Colours for species
  colours <- get_colours(vars = model_species, FG = FG)


  if(length(exp_str) > 0){
    ids <- unique(plot_data$.add_str_ID)
    plots <- lapply(cli_progress_along(1:length(ids), name = "Creating plot",
                                       format = paste0(
                                         "{pb_spin} Creating plot ",
                                         "[{pb_current}/{pb_total}]   ETA:{pb_eta}"
                                       )),
                    function(i){
                      data <- plot_data %>% filter(.data$.add_str_ID == ids[i])
                      visualise_effects_plot(data = data, prop = model_species,
                                          colours = colours, se = se, average = average)+
                        labs(subtitle = ids[i])
                    })
    if(length(plots) > 1){
      plot <- new("ggmultiplot", plots = plots, nrow = nrow, ncol = ncol)
    } else {
      plot <- plots[[1]]
    }
    cli::cli_alert_success("Created all plots.")
  } else {
    plot <- visualise_effects_plot(data = plot_data, prop = model_species,
                                colours = colours, se = se, average = average)
    cli::cli_alert_success("Created plot.")
  }
  return(plot)
}
