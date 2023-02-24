
#' @title Visualise effect of species loss or addition
#' @description This function helps to visualise the effect of the addition or loss of a species from a community on the response
#'
#' @importFrom ggplot2 geom_line labs scale_colour_manual geom_ribbon stat_summary scale_x_reverse facet_wrap
#' @importFrom dplyr %>% select filter all_of everything slice mutate bind_rows
#' @importFrom PieGlyph geom_pie_glyph
#' @importFrom rlang :=
#'
#' @param model A Diversity Interactions model object fit by using the \code{\link[DImodels:DI]{DI()}} function from the \code{\link[DImodels:DImodels-package]{DImodels}} package.
#' @param communities A dataframe specifying communities of interest for which user wants visualise the effect of species loss or addition. If left blank, the communities from the original data used to fit the model would be selected
#' @param species_interest A character vector specifying the species for which to visualise the effect of addition or loss on the response. If left blank, depending on the number of communities in the design, either the original communities or the equi-proportional community would be selected.
#' @param effect A character string with one value from c("addition", "loss" ) whether to visualise the effect of addition or loss of a species.
#' @param .equi A boolean variable indicating whether to include the equi-proportional community in the plot
#' @param .design A boolean variable to indicate whether to include the communities in the design in the plot
#' @param .se A boolean variable indicating whether to plot confidence intervals associated with the effect of species addition or loss
#' @param .average A boolean variable indicating whether to add a line describing the "average" effect of species addtion or loss
#' @param FG A character vector specifying the functional grouping of the species in the design. Species belonging to the same functional group will be assigned with different shades of the same colour. The user can manually specify a character vector giving the functional group each species belongs to. If left empty the functional will try to get a functional grouping from the original \code{\link[DImodels:DI]{DI}} model object.
#'
#' @return A ggplot object
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
#' effects_plot(model = mod)
#'
#' ## Add custom communities to plot instead of design communities
#' effects_plot(model = mod, .design = FALSE,
#'              communities = data.frame(p1 = c(0.9, 0.1),
#'                                       p2 = c(0.1, 0.9)))
#'
#' ## Add uncertainity on plot
#' effects_plot(model = mod,
#'              communities = data.frame(p1 = c(0.9, 0.1),
#'                                       p2 = c(0.1, 0.9)),
#'              .se= TRUE)
#'
#' ## Visualise effect of species loss for particular species
#' effects_plot(model = mod, effect = "loss",
#'              species_interest = c("p1", "p3"))
effects_plot <- function(model, communities = NULL, species_interest = NULL, effect = c('addition', 'loss'), .equi = FALSE, .design = TRUE, .se = FALSE, .average = TRUE, FG = NULL){

  # Sanity checks
  # Ensure model is a DImodels object
  if(!inherits(model, 'DI')){
    stop('Model should be a DImodels object')
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
  if(is.null(communities) && (!is.null(.equi) && !.equi) && (!is.null(.design) && !.design)){
    stop("Specify a dataset of communities to show effects for or set either or both of '.equi' or '.design' to TRUE.")
  }

  # If the user has not specified neither of communities, equi or model_data. Then default to plotting the design communities
  if(is.null(communities) & is.null(.equi) & is.null(.design)){
    design_points <- nrow(unique(original_data[, model_species]))
    # By default the model will plot all communities in the original data. However this might not be feasible when there are a lot
    # of communities in the design. If that's the case then plot only the equi-proportional mixtures
    if(design_points <= 100 ){
      .design <- T
    } else {
      .equi <- T
    }
  }

  # Create empty dataframe to append communities to plot
  plot_communities <- data.frame()

  # If plotting design communities then add those to the data
  if(!is.null(.design)){
    if(!is.logical(.design)){
      stop(glue::glue('.design takes boolean value of {TRUE} or {FALSE} specifying whether or not to plot the design communities'))
    } else {
      if(.design){
        # Extract the design communities from data used to fit the model
        design <- original_data[, model_species] %>% unique()

        # Append design communities to communities dataframe
        plot_communities <- rbind(plot_communities, design)
      }
    }
  }

  # If plotting the equi-proportional mixtures then add those to the data
  if(!is.null(.equi)){
    if(!is.logical(.equi)){
      stop(glue::glue('.equi takes boolean value of {TRUE} or {FALSE} specifying whether or not to plot the design communities'))
    } else {
      if(.equi){
        # Create the equi-proportional mixture
        nSpecies <- length(model_species)
        prop <- 1/nSpecies
        equi <-  data.frame(matrix(prop, nrow = 1, ncol = nSpecies))
        colnames(equi) <- model_species

        # Append the equi proportional mixture to communities dataframe
        plot_communities <- rbind(plot_communities, equi)
      }
    }
  }


  # If user has manually specified the communities then add those to data to be plotted
  if(!is.null(communities)){
    # Ensure communities are specified as a dataframe and in proper format
    if(!inherits(communities,'data.frame')){
      stop("Specify communities to show effects for in the plot as a data frame with proportions of the different species")
    } else {

      sp_abs <- !model_species %in% colnames(communities)
      # If the column names of the communities dataframe not match the species present in the model then notify user asking them to change the column names
      if(all(sp_abs)){
        stop(glue::glue("The column names of the communities data frame should be same as the names of the species specified when fitting the model.
                  Update the column names to ensure they match {paste(model_species, collapse = ', ')}"))
      }
      # If any species present in the model is not specified in the data, then notify user and assume its proportion to be 0
      if(any(sp_abs)){
        message(glue::glue("{paste(model_species[sp_abs], collapse = ', ')} weren't present in the data specified, assuming their proportions to be 0."))
      }
      communities[, model_species[sp_abs]] <- 0

      # If there are any additional columns in the communities data.frame other than the species proportions. Notify user and ignore those
      extra_cols <- !colnames(communities) %in% model_species
      if(any(extra_cols)){
        message(glue::glue("{paste(colnames(communities)[extra_cols], collapse = ', ')} were additional columns present in the communities dataframe, ignoring them."))
      }
      # Append the user specified communities to communities dataframe
      plot_communities <- rbind(plot_communities, communities[, model_species])
    }
  }

  # Ensure the proportions of all communities in the communities sum approx to 1.
  if(!all(near(rowSums(plot_communities), 1))){
    stop('Proportions of communities in communities should sum to 1')
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
    stop('species_interest should be present in the DImodels')
  }

  # Ensure effect is only one of 'addition' or 'loss'
  effect <- match.arg(effect)

  # No plots can be created if user is looking for effect of species loss but none of the communities specified contain
  # the species of interest
  if(effect == 'loss'){
    if(all(plot_communities[, species_interest] == 0)){
      stop("Specify communities which contain some positive proportion of the species of interest to see the effect of its loss. Can't visualize effect of loss of species if no communities contain the species of interest")
    }
  }

  # Similarly no plots can be created if the user is trying to find the effect of adding a species to it's monoculture
  if(effect == 'addition'){
    if(length(species_interest) == 1 & all(plot_communities[, species_interest] == 1)){
      stop("Specify communities other than the monoculture of the species of interest. Can't visualize effect of addition of species if only monocultures of species of interest are present.")
    }
  }

  # Get functional groups
  if(is.null(FG)){
    FG <- eval(model$DIcall$FG)
  }

  # Colours for species
  colours <- get_colours(species = model_species, FG = FG)

  # Create data to create effects plot
  plot_data <- lapply(species_interest, function(x){
    if(effect == 'addition'){
      # Interval or calculating effecct of adding species
      pvals <- seq(0, 1, length.out = 101)
      # Choose all species other than the species of interest to perform extrapolation to see effect of adding species of interest
      pothers <- plot_communities %>% select(-all_of(x))
      # Filter out the monoculture of the species of interest as can't add species to that
      pothers <- pothers[rowSums(pothers)!=0, ]
      # Rescale species proportions to ensure they sum to 1
      pothers <- pothers/rowSums(pothers)
      pother_names <- colnames(pothers)
      # Filter duplicate communities after the rescaling
      pothers <- unique(pothers)
    } else {
      # If we are looking at species loss then we don't need to look at monocultures of species of interest or in mixture which don't contain
      # our species of interest
      pothers <- plot_communities %>% filter(!(!!sym(x) %in% c(0, 1)))
      pother_names <- colnames(pothers %>% select(-all_of(x)))
      # Filter any duplicate communities
      pothers <- unique(pothers)
    }

    # For each row in data create copies with species of interest having one value from pvals (0.01, 0.02, 0.03, ..., 0.98, 0.99, 1)
    data_expand <- lapply(seq_len(nrow(pothers)), function(i){
      # Intervals for calculating the effects of species loss
      if(effect == 'loss'){
        pvals <- seq(pothers[i, x], 0, length.out = 101)
      }
      pothers_expand <-  pothers %>%
        slice(i) %>%
        slice(rep(1, each = length(pvals))) %>%
        mutate(!!sym(x) := pvals,
               'Sp' = x,
               'Proportion' = pvals,
               group = i) %>%
        select(all_of(model_species), everything())

      # Rescale the proportions of the other species accordingly to ensure the proportions sum to 1
      pothers_expand[, pother_names] <- pothers_expand[, pother_names]*(1 - pvals)/rowSums(pothers_expand[, pother_names])

      # Make predictions for the communities and calculate the uncertainity for the prediction
      preds <- predict(model, newdata = pothers_expand, se.fit = T)
      pothers_expand <- pothers_expand %>%
        mutate('Pred' = preds$fit,
               'Lower' = preds$fit - 1.96*preds$se.fit,
               'Upper' = preds$fit + 1.96*preds$se.fit)

      # Calculate the marginal effect of adding a species for producing the marginal plots
      dy <- diff(pothers_expand$Pred)/diff(pvals)
      pothers_expand <- pothers_expand %>%
        mutate('marginal' = c(dy, dy[length(dy)]),
               'threshold' = pvals[abs(.data$marginal) == min(abs(.data$marginal))][1],
               'Effect' = ifelse(!!sym(x) < .data$threshold, 'Negative', 'Positive'))
      pothers_expand

    }) %>% bind_rows()
    data_expand
  }) %>% bind_rows()

  # Get the subset of the data to create the pie-charts (starting and ending points for a curve)
  pie_data <- plot_data %>%
    group_by(.data$Sp, .data$group) %>%
    filter(.data$Proportion %in% c(min(.data$Proportion), max(.data$Proportion))) %>%
    ungroup() %>%
    group_by(.data$Sp) %>%
    distinct(.data$Pred, .keep_all = T)

  # Create canvas for plot
  plot <- ggplot(plot_data, aes(x = .data$Proportion, y= .data$Pred))+
    theme_bw()

  # Add ribbons for uncertainty of prediction
  if(.se){
    plot <- plot +
      geom_ribbon(aes(ymin = .data$Lower, ymax = .data$Upper,
                      group = .data$group),
                  colour = 'grey', linetype = 'dashed',
                  linewidth = 0.7, alpha = 0.25)
  }

  # Add line tracing the effect of adding a particular species to the communities
  plot <- plot +
    geom_line(aes(group = .data$group), colour = 'black', alpha = 0.4)

  # If user specified to chose to see average effect of adding species then add that to plot
  if(.average){
    # Round species proportion to 2 decimal places and calculate average effect.
    # Rounding done because for species loss because some proportions have recurring decimals
    avg_data <- plot_data %>%
      mutate('Proportion' = round(.data$Proportion, 2))

    # Add line showing average effect
    plot <- plot +
      stat_summary(group = 'Sp', fun=mean, geom="line", colour="black", linewidth = 1.25, data = avg_data)
  }

  # Add the pie-chart glyphs for identifying the communities
  plot <- plot +
    geom_pie_glyph(aes(group = .data$Sp), data = pie_data, radius = 0.3,
                   slices = model_species, colour = 'black')+
    scale_fill_manual(values = colours,
                      labels = model_species)

  # If looking at species loss then reverse the x-axis
  if(effect == 'loss'){
    plot <- plot + scale_x_reverse()
  }
  # Finally facet plot by each species of interest to ensure better visibility
  plot <- plot +
    facet_wrap(~.data$Sp)

  # Adjust plot aesthetics
  plot <- plot +
    theme(legend.position = 'top')

  return(plot)
}
