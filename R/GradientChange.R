#' @title Visualising average change in response over a diversity gradient.
#'
#' @description
#' Helper function for creating the data for the average richness plot.
#'
#' @importFrom tools toTitleCase
#'
#' @param data A data-frame consisting of species proportions and
#'             any other necessary variables to make predictions from
#'             `model` or `coefficients`.
#' @param prop A vector identifying the column-names or indices of the
#'             columns containing the proportions in `data`.
#' @inheritParams gradient_change
#' @inheritDotParams add_prediction -data
#'
#' @return A data-frame with the following columns appended at the end
#'  \describe{
#'    \item{.Pred}{The predicted response for each community.}
#'    \item{.Lower}{The lower limit of the prediction/confidence interval
#'                  for each observation.}
#'    \item{.Upper}{The upper limit of the prediction/confidence interval
#'                  for each observation.}
#'    \item{.Richness}{The richness (number of species) within each community.}
#'    \item{.Evenness}{The evenness (metric quantifying the relative abundance
#'                     of each species) of each community.}
#'    \item{.Gradient}{An character string defining the diversity gradient used
#'                     for averaging the response.}
#'    \item{.Avg}{The averaged value of the response for each unique value of
#'                the selected diversity gradient.}
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
#' ## By default response would be averaged on the basis of richness
#' head(gradient_change_data(data = sim2,
#'                           prop = c("p1", "p2", "p3", "p4"),
#'                           model = mod))
#'
#' ## Average response with respect to evenness
#' head(gradient_change_data(data = sim2,
#'                           prop = c("p1", "p2", "p3", "p4"),
#'                           model = mod,
#'                           gradient = "evenness"))
gradient_change_data <- function(data, prop,
                                 gradient = c("richness", "evenness"), ...){
  if(missing(data)){
    cli::cli_abort(c("{.var data} cannot be empty.",
                     "i" = "Specify a data frame or tibble containing the
                            variable proportions which to use for calculating
                            the average change across a diversity gradient."))
  }

  #Sanity Checks
  if(missing(prop)){
    cli::cli_abort(c("{.var prop} cannot be empty.",
                     "i" = "Specify the column-names or indices of the columns
                           containing the variable proportions in {.var data}."))
  }

  sanity_checks(data = data, prop = prop)
  # choose variable to group by
  gradient <- match.arg(gradient)

  # Calculate predictions
  plot_data <- add_prediction(data = data, ...)

  # Get names of columns containing variable proportions
  species_names <- data %>% select(all_of(prop)) %>% colnames()

  # Choose diversity gradient to group data on
  group_var <- paste0(".", tools::toTitleCase(gradient))

  # Get name of column containing predictions
  # dotArgs <- rlang::dots_values(...)
  # pred_var <- if (!is.null(dotArgs$pred_name)) dotArgs$pred_name else ".Pred"

  # Add the average change across specified gradient
  plot_data <- plot_data %>%
    mutate(".Richness" = get_richness(.data, species_names),
           ".Evenness" = get_evenness(.data, species_names),
           ".Gradient" = group_var) %>%
    group_by(!! sym(group_var)) %>%
    mutate('.Avg' = mean(.data$.Pred)) %>%
    ungroup()

  cli::cli_alert_success("Finished data preparation")
  return(plot_data)
}


#' @title Visualising average change in response over a diversity gradient.
#'
#' @description
#' Helper for plotting average response at each level of species richness
#'
#' @param data A data-frame which is the output of the
#'             `\link{gradient_change_data}` function, consisting of the
#'             predicted response averaged over a particular diversity gradient.
#' @param pie_data A subset of data-frame specified in `data`, to visualise
#'                 the individual data-points as pie-glyphs showing the
#'                 relative proportions of the variables in the data-point.
#' @inheritParams gradient_change
#' @inheritParams gradient_change_data
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
#' data(sim4)
#'
#' ## Fit model
#' mod <- glm(response ~ 0 + (p1 + p2 + p3 + p4 + p5 + p6)^2, data = sim4)
#'
#' ## Create data
#' ## By default response would be averaged on the basis of richness
#' plot_data <- gradient_change_data(data = sim4,
#'                                   prop = c("p1", "p2", "p3",
#'                                            "p4", "p5", "p6"),
#'                                   model = mod)
#' gradient_change_plot(data = plot_data,
#'                      prop = c("p1", "p2", "p3", "p4", "p5", "p6"))
#'
#' ## If prop is not specified then the observations will be shows as points
#' gradient_change_plot(data = plot_data)
#'
#' ## Average response with respect to evenness
#' plot_data <- gradient_change_data(data = sim4,
#'                           prop = c("p1", "p2", "p3", "p4", "p5", "p6"),
#'                           model = mod,
#'                           gradient = "evenness")
#' gradient_change_plot(data = plot_data,
#'                      prop = c("p1", "p2", "p3", "p4", "p5", "p6"))
#'
#' ## Change colours of the pie-slices
#' gradient_change_plot(data = plot_data,
#'                      prop = c("p1", "p2", "p3", "p4", "p5", "p6"),
#'                      pie_colours = c("darkolivegreen1", "darkolivegreen4",
#'                                      "orange1", "orange4",
#'                                      "steelblue1", "steelblue4"))
#'
#'
#' ## Manually specify only specific communities to be shown as pie-chart
#' ## glyphs using `pie_data`.
#' gradient_change_plot(data = plot_data,
#'                      prop = c("p1", "p2", "p3", "p4", "p5", "p6"),
#'                      pie_data = plot_data %>% filter(.Richness %in% c(1, 6)))
#' ## Note: It is important for the data specified in
#' ## `pie_data` to have the .Pred and .Gradient columns.
#' ## So the best use case for this parameter is to accept
#' ## a subset of the data specified in `data`.
gradient_change_plot <- function(data, prop = NULL,
                                 pie_data = NULL,
                                 pie_colours = NULL){

  # Ensure inputs are appropriate
  sanity_checks(data = data, prop = prop, colours = pie_colours)

  # If pie_data is specified then prop cannot be NULL
  if(!is.null(pie_data) && is.null(prop)){
    cli::cli_abort(c("{.var prop} cannot be {.pkg NULL} when {.var pie_data}
                     is specified.",
                   "i" = "Specify the column-names or indices of the columns
                         containing the variable proportions in {.var pie_data}
                         to show as pie-glyphs on the plot. "))
  }

  # If prop is specified without pie_data, assume data is pie_data
  if(!is.null(prop) && is.null(pie_data)){
      pie_data <- data
  }

  # Check prop is present in pie_data too
  sanity_checks(data = pie_data, prop = prop)

  # Deciding whether to dodge points
  gradient_var <- data$.Gradient[1]
  position <- if(gradient_var == ".Richness") position_dodge(0.5) else position_identity()

  # Create plot
  plot <- ggplot(data = data, aes(x = !! sym(gradient_var), y = .data$.Pred))+
    geom_line(aes(y = .data$.Avg), linewidth = 1)+
    geom_point(aes(group = .data$.Evenness),
               colour = '#555555',
               position = position)+
    labs(x = strsplit(gradient_var, split = "")[[1]][-1] %>%
               paste0(collapse = ""),
         y = 'Predicted Response',
         fill = 'Species')+
    theme_DI()

  if(!is.null(pie_data)){
    # Colours for the pie-glyph slices
    if(is.null(pie_colours)){
      pie_colours <- get_colours(prop)
    } else {
      # Ensure there are same number of colours as variables specified in prop
      if(length(pie_colours) != length(prop)){
        cli::cli_abort(c("The number of {.var pie_colours} specified should
                         be same as the number of variables specified in
                         {.var prop}.",
                         "i" = "There are {length(prop)} values in {.var prop}
                               but only {length(pie_colours)} colour{?s}
                               {?was/were} specified."))
      }
    }

    plot <- plot +
      geom_pie_glyph(aes(group = .data$.Evenness),
                     slices = prop,
                     colour = 'black', data = pie_data,
                     position = position)+
      scale_fill_manual(values = pie_colours, name = "Species")
  }
  return(plot)
}



#' @title Visualise average response at each level of richness
#'
#' @importFrom dplyr mutate %>% group_by distinct across all_of
#' @importFrom ggplot2 ggplot geom_line aes geom_point position_dodge
#'                     position_identity labs theme theme_bw
#' @importFrom PieGlyph geom_pie_glyph
#' @importFrom rlang !!! !! sym .data
#' @importFrom DImodels DI_data_E_AV
#'
#' @param gradient Diversity gradient to show on the X-axis, one of
#'                 "richness" or "evenness". Defaults to "richness".
#' @param communities Communities which are used to calculate the average
#'                    response. Accepts one of "original" or "equi".
#'                    "original" (the default) calculates the average using
#'                    the original communities from the raw data used to fit
#'                    the \code{\link[DImodels:DI]{DI}} model while "equi"
#'                    uses all possible equi-proportional communities at each
#'                    level of richness to calculate the average.
#' @param threshold For high levels of richness (>20) it might not be feasible
#'                  to compute and store all equi-proportional communities.
#'                  This parameter defines the number of communities randomly
#'                  selected at each level of richness.
#' @inheritParams response_contributions
#' @inheritParams model_diagnostics
#'
#' @inherit response_contributions return
#'
#' @export
#'
#' @examples
#' ## Load DImodels package to fit the model
#' library(DImodels)
#'
#' ## Load data
#' data(sim4)
#'
#' ## Fit DI model
#' mod <- DI(prop = 3:8, DImodel = 'AV', data = sim4, y = 'response')
#'
#' ## Create visualisation
#'
#' ## By default, 'richness' is the gradient and commmunities from the
#' ## raw data are used to calcualte average response
#' gradient_change(model = mod)
#'
#' ## Can also calculate average response using all possible
#' ## equi-proportional communities at each level of richness
#' gradient_change(model = mod, communities = 'equi')
#'
#' ## Can also plot across evenness
#' data(sim5)
#' mod1 <- DI(prop = paste0('p', 1:9), DImodel = 'AV',
#'           data = sim5, y = 'response', estimate_theta = TRUE)
#' gradient_change(model = mod1, gradient = 'evenness')
gradient_change <- function(model, gradient = c('richness', 'evenness'),
                            communities = c('original', 'equi'),
                            exp_str = list(),
                            threshold = NULL,
                            pie_colours = NULL,
                            nrow = 0, ncol = 0){
  # Ensure specified model is fit using the DI function
  if(missing(model) || !inherits(model, "DI")){
    model_not_DI(call_fn = "gradient_change")
  }

  # choose variable to plot on x-axis
  gradient <- match.arg(gradient)

  # choose communities to calculate average from
  communities <- match.arg(communities)

  # Get original data used to fit the model
  original_data <- model$original_data

  # Get all species in the model
  model_species <- eval(model$DIcall$prop)
  # If species were specified as column indices extract names
  if(is.numeric(model_species)){
    model_species <- colnames(original_data)[model_species]
  }

  # If threshold is missing set default
  if(is.null(threshold)){
    threshold <- round(1000000/length(model_species))
  }

  # Sanity checks for threshold
  sanity_checks(colours = pie_colours,
                numerics = list("threshold" = threshold),
                unit_lengths = list("threshold" = threshold))


  if(communities == 'equi'){
    plot_communities <- get_equi_comms(nSpecies = length(model_species),
                                       species = model_species,
                                       threshold = threshold)
  } else {
    plot_communities <- original_data
  }

  # Ensure experimental structure are specified correctly
  exp_str <- check_exp_str(model = model, exp_str = exp_str)

  if(length(exp_str) > 0){
    plot_communities <- add_exp_str(data = plot_communities,
                                    exp_str = exp_str)
  }

  # Create data for plotting
  plot_data <- gradient_change_data(data = plot_communities,
                                     model = model,
                                     prop = model_species,
                                     gradient = gradient)
  # If exp_str is specified then unique plot for each combination of exp_str
  if(length(exp_str) > 0){
    ids <- unique(plot_data$.add_str_ID)
    plots <- lapply(cli_progress_along(1:length(ids), name = "Creating plot",
                                       format = paste0(
                                         "{pb_spin} Creating plot ",
                                         "[{pb_current}/{pb_total}]  ETA:{pb_eta}"
                                       )),
                    function(i){
                      data <- plot_data %>% filter(.data$.add_str_ID == ids[i])
                      prop <- if(length(model_species) > 10 || nrow(data) > 1024) NULL else model_species
                      gradient_change_plot(data = data, prop = prop,
                                           pie_colours = pie_colours)+
                        labs(subtitle = ids[i])
                    })
    if(length(plots) > 1){
      plot <- new("ggmultiplot", plots = plots, nrow = nrow, ncol = ncol)
    } else {
      plot <- plots[[1]]
    }
    cli::cli_alert_success("Created all plots.")
  # Else return a single ggplot
  } else {
    prop <- if(length(model_species) > 10 || nrow(plot_data) > 1024) NULL else model_species
    plot <- gradient_change_plot(data = plot_data, prop = prop,
                                 pie_colours = pie_colours)
    cli::cli_alert_success("Created plot.")
  }
  return(plot)
}
