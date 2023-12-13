#' @title Average response over diversity gradient
#'
#' @description
#' Helper function for creating the data to visualise the average
#' response over a diversity gradient. The average can be calculated
#' either over the "richness" or "evenness" gradients. The average is
#' calculated from all communities present at a given level of the
#' chosen diversity gradient in `data`. The output of this
#' function can be passed to the \code{\link{gradient_change_plot}} function
#' to visualise the results.
#'
#' @importFrom tools toTitleCase
#'
#' @param data A data-frame consisting of variable proportions and
#'             any other necessary variables to make predictions from
#'             `model` or `coefficients`.
#' @param prop A vector identifying the column-names or indices of the
#'             columns containing the variable proportions in `data`.
#' @inheritParams gradient_change
#' @inheritDotParams add_prediction -data
#'
#' @return The data-frame with the following columns appended at the end
#'  \describe{
#'    \item{.Pred}{The predicted response for each obvervation.}
#'    \item{.Lower}{The lower limit of the prediction/confidence interval
#'                  for each observation.}
#'    \item{.Upper}{The upper limit of the prediction/confidence interval
#'                  for each observation.}
#'    \item{.Richness}{The richness (number of non-zero compositional variables)
#'                     within each observation.}
#'    \item{.Evenness}{The evenness (metric quantifying the relative abundance
#'                     of each compositional variable) within each observation.}
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
#' library(dplyr)
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
#'
#' ## Additional variables can also be added to the data by either specifying
#' ## them directly in the `data` or by using the `add_var` argument
#' ## Refit model
#' sim2$block <- as.numeric(sim2$block)
#' new_mod <- update(mod, ~. + block, data = sim2)
#' ## This model has block so we can either specify block in the data
#' subset_data <- sim2[c(1,5,9,11), 2:6]
#' subset_data
#' head(gradient_change_data(data = subset_data,
#'                           prop = c("p1", "p2", "p3", "p4"),
#'                           model = mod,
#'                           gradient = "evenness"))
#' ## Or we could add the variable using `add_var`
#' subset_data <- sim2[c(1,5,9,11), 3:6]
#' subset_data
#' head(gradient_change_data(data = subset_data,
#'                           prop = c("p1", "p2", "p3", "p4"),
#'                           model = new_mod,
#'                           gradient = "evenness",
#'                           add_var = list(block = c(1, 2))))
#' ## The benefit of specifying the variable this way is we have an ID
#' ## columns now called `.add_str_ID` which could be used to create a
#' ## separate plot for each value of the additional variable
#'
#'
#' ## Model coefficients can also be used, but then user would have
#' ## to specify the data with all columns corresponding to each coefficient
#' coef_data <- sim2 %>%
#'                mutate(`p1:p2` = p1*p2, `p1:p3` = p1*p2, `p1:p4` = p1*p4,
#'                       `p2:p3` = p2*p3, `p2:p4` = p2*p4, `p3:p4` = p3*p4) %>%
#'                select(p1, p2, p3, p4,
#'                       `p1:p2`, `p1:p3`, `p1:p4`,
#'                       `p2:p3`, `p2:p4`, `p3:p4`) %>%
#'                slice(1,5,9,11)
#' print(coef_data)
#' print(mod$coefficients)
#' gradient_change_data(data = coef_data,
#'                      prop = c("p1", "p2", "p3", "p4"),
#'                      gradient = "evenness",
#'                      coefficients = mod$coefficients,
#'                      interval = "none")
gradient_change_data <- function(data, prop, add_var = list(),
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

  # Ensure experimental structure are specified correctly
  dotArgs <- rlang::dots_values(...)
  model <- if (!is.null(dotArgs$model)) dotArgs$model else NULL
  if(!is.null(model)){
    add_var <- check_add_var(model = model, add_var = add_var)
  }

  # Add any experimental structures
  if(length(add_var) > 0){
    data <- add_add_var(add_var = add_var, data = data)
  }

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


#' @title Average response over diversity gradient
#'
#' @description
#' Helper function for plotting the average response at each level of the
#' diversity gradient. The output of the \code{\link{gradient_change_data}}
#' function should be passed here to visualise the average response. The
#' proportions of the compositional variables present in the data would be
#' shown on the plot using `\code{\link[PieGlyph:PieGlyph-package]{pie-glyphs}}`.
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
#'                                   prop = c("p1", "p2", "p3",
#'                                            "p4", "p5", "p6"),
#'                                   model = mod,
#'                                   gradient = "evenness")
#' gradient_change_plot(data = plot_data,
#'                      prop = c("p1", "p2", "p3", "p4", "p5", "p6"))
#'
#' ## Change colours of the pie-slices using `pie_colours`
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
#'
#' ## Use `facet_var` to facet the plot on an additional variable
#' gradient_change_plot(data = plot_data,
#'                      prop = c("p1", "p2", "p3", "p4", "p5", "p6"),
#'                      pie_data = plot_data %>% filter(.Richness %in% c(1, 6)),
#'                      facet_var = "treatment")
#'
#' ## If `add_var` was used during the data preparation step then
#' ## multiple plots will be produced and can be arranged using nrow and ncol
#' new_mod <- update(mod, ~. + treatment, data = sim4)
#' plot_data <- gradient_change_data(data = sim4[c(seq(1, 18, 3), 19:47), -2],
#'                                   prop = c("p1", "p2", "p3",
#'                                            "p4", "p5", "p6"),
#'                                   model = new_mod,
#'                                   add_var = list("treatment" = c(50, 250)))
#' ## Create plot arranged in 2 rows
#' gradient_change_plot(data = plot_data,
#'                      prop = c("p1", "p2", "p3", "p4", "p5", "p6"),
#'                      pie_data = plot_data %>% filter(.Richness %in% c(1, 6)),
#'                      nrow = 2)
gradient_change_plot <- function(data, prop = NULL,
                                 pie_data = NULL,
                                 pie_colours = NULL,
                                 facet_var = NULL,
                                 nrow = 0, ncol = 0){

  # If add_var is specified then unique plot for each combination of add_var
  if(check_col_exists(data, ".add_str_ID")){
    ids <- unique(data$.add_str_ID)
    plots <- lapply(cli_progress_along(1:length(ids), name = "Creating plot",
                                       format = paste0(
                                         "{pb_spin} Creating plot ",
                                         "[{pb_current}/{pb_total}]  ETA:{pb_eta}"
                                       )),
                    function(i){
                      data_iter <- data %>% filter(.data$.add_str_ID == ids[i])
                      gradient_change_plot_internal(data = data_iter, prop = prop,
                                                    pie_colours = pie_colours,
                                                    facet_var = facet_var,
                                                    pie_data = pie_data)+
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
    plot <- gradient_change_plot_internal(data = data, prop = prop,
                                          pie_colours = pie_colours,
                                          facet_var = facet_var,
                                          pie_data = pie_data)
    cli::cli_alert_success("Created plot.")
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
#' @param data A dataframe specifying communities of interest for which user
#'             wants to visualise the gradient. If left blank, the value specified
#'             in `communities` will be used to create data for the user.
#' @param pie_data A subset of data-frame specified in `data`, to visualise
#'                 the individual data-points as pie-glyphs showing the
#'                 relative proportions of the variables in the data-point.
#' @param gradient Diversity gradient to show on the X-axis, one of
#'                 "richness" or "evenness". Defaults to "richness".
#' @inheritParams prediction_contributions
#' @inheritParams model_diagnostics
#'
#' @inherit prediction_contributions return
#'
#' @export
#'
#' @examples
#' ## Load DImodels package to fit the model
#' library(DImodels)
#' library(dplyr)
#'
#' ## Load data
#' data(sim4)
#'
#' ## Fit DI model
#' mod <- DI(prop = 3:8, DImodel = 'AV', data = sim4, y = 'response') %>%
#'          suppressWarnings()
#'
#' ## Create visualisation
#'
#' ## By default, 'richness' is the gradient and communities from the
#' ## raw data are used to calculate average response
#' gradient_change(model = mod)
#'
#' ## Specify custom data
#' gradient_change(model = mod, data = sim4 %>% filter(richness <= 4))
#'
#' ## Create plot for all equi-proportional communities at a
#' ## given level of richness
#' plot_data <- get_equi_comms(6, variables = paste0("p", 1:6))
#' gradient_change(model = mod, data = plot_data)
#'
#' ## Can also plot average response across evenness
#' gradient_change(model = mod, gradient = 'evenness')
#'
#' ## Change colours of the pie-slices using `pie_colours`
#' gradient_change(model = mod,
#'                 pie_colours = c("darkolivegreen1", "darkolivegreen4",
#'                                 "orange1", "orange4",
#'                                 "steelblue1", "steelblue4"))
#'
#' ## Manually specify only specific communities to be shown as pie-chart
#' ## glyphs using `pie_data`.
#' gradient_change(model = mod,
#'                 pie_data = sim4 %>% filter(richness %in% c(1, 6)))
#'
#' ## Use `facet_var` to facet the plot on an additional variable
#' gradient_change(model = mod,
#'                 pie_data = sim4 %>% filter(richness %in% c(1, 6)),
#'                 facet_var = "treatment")
#'
#' ## Use `add_var` to add additional variables independent of the compositions
#' ## Multiple plots will be produced and can be arranged using nrow and ncol
#' ## Create plot arranged in 2 rows
#' gradient_change(model = mod,
#'                 data = sim4[, -2],
#'                 add_var = list("treatment" = c(50, 250)),
#'                 pie_data = sim4[, -2] %>% filter(richness %in% c(1, 6)),
#'                 nrow = 2)
#'
#' ## Specify `plot = FALSE` to not create the plot but return the prepared data
#' gradient_change(model = mod, plot = FALSE,
#'                 pie_data = sim4 %>% filter(richness %in% c(1, 6)),
#'                 facet_var = "treatment")
gradient_change <- function(model, data = NULL,
                            gradient = c('richness', 'evenness'),
                            add_var = list(),
                            plot = TRUE,
                            pie_data = NULL,
                            pie_colours = NULL,
                            facet_var = NULL,
                            nrow = 0, ncol = 0){
  # Ensure specified model is fit using the DI function
  if(missing(model) || !inherits(model, "DI")){
    model_not_DI(call_fn = "gradient_change")
  }

  # choose variable to plot on x-axis
  gradient <- match.arg(gradient)

  # Get original data used to fit the model
  original_data <- model$original_data

  # Get all species in the model
  model_species <- eval(model$DIcall$prop)
  # If species were specified as column indices extract names
  if(is.numeric(model_species)){
    model_species <- colnames(original_data)[model_species]
  }


  # If data is not specified then create data
  if(is.null(data)){
      data <- original_data
  }

  # Sanity checks for colours and data
  sanity_checks(data = data, colours = pie_colours)


  # Ensure experimental structure are specified correctly
  add_var <- check_add_var(model = model, add_var = add_var)
  if(length(add_var) > 0){
    data <- add_add_var(data = data,
                        add_var = add_var)

    if(!is.null(pie_data)){
      pie_data <- add_add_var(data = pie_data,
                              add_var = add_var)
    }
  }

  # Create data for plotting
  plot_data <- rlang::try_fetch(gradient_change_data(data = data,
                                                     model = model,
                                                     prop = model_species,
                                                     gradient = gradient),
                                error = function(cnd)
                                  rlang::abort("The following error was encountered while preparing `data` for plotting.",
                                               parent = cnd))
  if(!is.null(pie_data)){
    pie_data <- rlang::try_fetch(suppressMessages(gradient_change_data(
                                                      data = pie_data,
                                                      model = model,
                                                      prop = model_species,
                                                      gradient = gradient,
                                                      )),
                                 error = function(cnd)
                                   rlang::abort("The following error was encountered while preparing `pie_data` for plotting.",
                                                parent = cnd))
  }

  prop <- if(length(model_species) > 10 || nrow(plot_data) > 1024) NULL else model_species

  if(isTRUE(plot)){
    plot <- gradient_change_plot(data = plot_data,
                                 prop = prop,
                                 pie_data = pie_data,
                                 pie_colours = pie_colours,
                                 facet_var = facet_var,
                                 nrow = nrow, ncol = ncol)
    return(plot)
  } else {
    return(plot_data)
  }
}

#' @keywords internal
#' Internal function for the gradient change plot
#'
#' @usage NULL
NULL
gradient_change_plot_internal <- function(data, prop = NULL,
                                 pie_data = NULL,
                                 pie_colours = NULL,
                                 facet_var = NULL){

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

  # Choose a proper scale for x-axis
  xlab <- strsplit(gradient_var, split = "")[[1]][-1] %>% paste0(collapse = "")
  pie_max <- ifelse(is.null(pie_data) || !check_col_exists(pie_data, ".Richness"),
                    -Inf, max(pie_data$.Richness))
  xmax <- max(max(data$.Richness), pie_max)
  x_scale <- if(gradient_var == ".Richness")
    ggplot2::scale_x_continuous(breaks = seq(1:xmax))
  else
    ggplot2::scale_x_continuous()

  # Caption
  if(gradient_var == ".Richness"){
    caption <- paste0("The solid black line connects the average predicted
                      response of all points shown at a given level of ", tolower(xlab),
                      ".\n The points at a given level of richness are jittered according
                      to the relative abundances of compositional variables in an observation.")
  } else {
    caption <- paste0("The solid black line connects the average predicted
                      response of all points shown at a given level of ", tolower(xlab), ".")
  }

  # Create plot
  plot <- ggplot(data = data, aes(x = !! sym(gradient_var), y = .data$.Pred))+
    geom_line(aes(y = .data$.Avg), linewidth = 1)+
    geom_point(aes(group = .data$.Evenness),
               colour = '#555555',
               position = position)+
    labs(x = xlab,
         y = "Predicted Response",
         fill = "Species",
         caption = caption)+
    theme_DI()+
    x_scale

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

  if(!is.null(facet_var)){
    plot <- add_facet(plot, data, facet_var)
  }

  return(plot)
}
