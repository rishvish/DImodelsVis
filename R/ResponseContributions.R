#' @title Splitting predicted response into contributions by
#'        individual model coefficients
#'
#' @description
#' The data preparation function for splitting the predicted response
#' from a regression model into contributions by the coefficients in the model.
#'
#' @importFrom stats coef formula model.matrix
#' @importFrom dplyr bind_cols ends_with row_number
#' @importFrom rlang :=
#'
#' @inheritParams response_contributions
#' @inheritParams add_prediction
#'
#' @return A data-frame with the following columns
#'  \describe{
#'    \item{.Community}{An identifier column to discern each
#'                      observation in the data.}
#'    \item{.Pred}{The predicted repsonse for each observation.}
#'    \item{.Lower}{The lower limit of the prediction interval for
#'                  each observation.}
#'    \item{.Upper}{The lower limit of the prediction interval for
#'                  each observation.}
#'    \item{.Contributions}{An identifier describing the name of the
#'                          coefficient contributing to the response.}
#'    \item{.Value}{The contributed value of the respective coefficient/group
#'                  to the total prediction.}
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
#' response_contributions_data(data = sim2[c(1,5,9,11), ],
#'                             model = mod)
#'
#' ## Specific coefficients can also be grouped together
#' ## Either by their indices in the model coefficient vector
#' response_contributions_data(data = sim2[c(1,5,9,11), ],
#'                             model = mod,
#'                             groups = list("Interactions" = 5:10))
#' ## Or by specifying the coefficient names as character strings
#' response_contributions_data(data = sim2[c(1,5,9,11), ],
#'                             model = mod,
#'                             groups = list("p1_Ints" = c("p1:p2",
#'                                                         "p1:p3",
#'                                                         "p1:p4")))
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
#' response_contributions_data(data = coef_data,
#'                             coefficients = mod$coefficients)
#' # Remember that it won't be possible to get prediction intervals
#' # when using the coefficients parameter
response_contributions_data <- function(data, model = NULL, coefficients = NULL,
                                        groups = list(), conf.level = 0.95,
                                        interval = "prediction"){
  if(missing(data)){
    cli::cli_abort(c("{.var data} cannot be empty.",
                     "i" = "Specify a data frame or tibble indicating the species communties which to use for calculating the average change across a diversity gradient."))
  }

  # Add predictions from model/coefficients to check things are fine
  pred_data <- add_prediction(data = data, model = model,
                              coefficients = coefficients,
                              interval = interval,
                              conf.level = conf.level)
  # Get names of columns containing species proportions
  # species_names <- data %>% select(all_of(species)) %>% colnames(.)

  # Branch here if model is specified
  if(!is.null(model)){
    form <- formula(model)
    form[2] <- NULL
    X_matrix <- model.matrix(form, data)
    coefficients <- coef(model)
  }
  # Branch here if regression coefficients are specified
  else if(!is.null(coefficients)){
    X_matrix <- as.matrix(sapply(data[, colnames(data)], as.numeric))
  }
  # Express the prediction into contribution from each coefficient in the model
  contr_data <- as.data.frame(sweep(X_matrix, MARGIN = 2, coefficients, `*`))

  # Add identifier for each row
  contr_data$.Community <- paste0("Obs", 1:nrow(contr_data))

  # Check coefficients can be properly grouped if groups are specified
  check_coeff_groupings(coefficients, groups)

  # Group contributions
  grouped_data <- data.frame(".Community" = contr_data$.Community)
  for(group in names(groups)){
    vars <- groups[[group]]
    grouped_data <- grouped_data %>%
      mutate(!!group := rowSums(contr_data %>% select(all_of(vars)), na.rm = TRUE))
  }
  # Add any variables left over after grouping the contributions
  grouped_data <- bind_cols(contr_data %>%
                              select(-.data$.Community) %>%
                              select(-all_of(unlist(groups))),
                            grouped_data)

  # Name of the contributions prediction is split into
  contributions <- grouped_data %>%
    select(-.data$.Community) %>%
    colnames()

  # If the data has an identifier for exp str then add that in and
  # reset observation numbers
  if(!is.null(data$.add_str_ID)){
    grouped_data$.add_str_ID <- data$.add_str_ID
    grouped_data <- grouped_data %>%
                      group_by(.data$.add_str_ID) %>%
                      mutate(".Community" = paste0("Obs", row_number())) %>%
                      ungroup()
  }

  # Add total predictions and uncertainity
  grouped_data <- grouped_data %>%
    mutate(.Community = fct_inorder(.data$.Community),
           .Pred = pred_data$.Pred,
           .Lower = pred_data$.Lower,
           .Upper = pred_data$.Upper) %>%
    # Converting to long format so it's easier to plot
    pivot_longer(cols = all_of(contributions),
                 names_to = '.Contributions', values_to = '.Value')

  cli::cli_alert_success("Finished data preparation.")
  return(grouped_data)
}

#' @rdname response_contributions_data
#' @title Splitting predicted response into contributions by
#'        individual model coefficients
#'
#' @description
#' The plotting function to visualise the predicted response from a
#' regression model as a stacked bar-chart showing the contributions
#' of each regression coefficient to the total predicted value. Requires
#' the output of the `\link{response_contributions_data}` as an input.
#'
#' @param data A data-frame which is the output of the
#'             `\link{response_contributions_data}` function, consisting of
#'             the predicted response split into the contributions by the
#'             different coefficients.
#' @param colours A character vector specifying the colours for the
#'                contributions of the different coefficients.
#'                If not specified, a default colour-scheme would be chosen,
#'                however it could be uninformative and it is recommended to
#'                manually choose contrasting colours for each coefficient
#'                group to render a more useful plot.
#' @param se A logical value indicating whether to show prediction intervals
#'           for predictions in the plot.
#'
#' @return A ggplot object
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
#' ## Create data for plotting
#' plot_data <- response_contributions_data(data = sim2[c(1,5,9,11,15,19,23), ],
#'                             model = mod)
#' ## Create plot
#' response_contributions_plot(data = plot_data)
#'
#' ## Choose distinct colours for groups of coefficients
#' ID_cols <- get_colours(4)
#' int_cols <- get_shades("#808080", 6)[[1]]
#' cols <- c(ID_cols, int_cols)
#'
#' response_contributions_plot(data = plot_data, colours = cols)
#'
#' ## Show prediction intervals
#' response_contributions_plot(data = plot_data, colours = cols, se = TRUE)
response_contributions_plot <- function(data, colours = NULL, se = FALSE){
  # Ensure data is specified
  if(missing(data)){
    cli::cli_abort(c("{.var data} cannot be empty.",
                     "i" = "Specify a data frame or tibble (preferably the
                            output of {.fn response_contributions_data})."))
  }

  # Colours for the bar segments
  if(is.null(colours)){
    rlang::warn(c("No colours were specified for the response contributions.",
                  "i" = "The default colours might not result in an informative
                         plot, consider choosing specific colours to contrast
                         the contributions of different groups in the response."),
                .frequency = "regularly", .frequency_id = "2")
    colours <- get_colours(levels(factor(data$.Contributions)))
  }

  plot <- ggplot(data, aes(x = .data$.Community, y = .data$.Value,
                           fill = fct_rev(fct_inorder(.data$.Contributions))))+
            geom_col(colour = "black")+
            theme_DI()+
            scale_fill_manual(values = rev(colours))+
            guides(fill = guide_legend(reverse = TRUE)) +
            labs(y = 'Predicted Response',
                 x = 'Community',
                 fill = 'Contributions')

  if(se){
    plot <- plot +
      geom_errorbar(aes(y = .data$.Pred,
                        ymin = .data$.Lower, ymax = .data$.Upper),
                    colour = "black")
  }
  # if(FALSE){
  #   plot <- plot +
  #     facet_grid(~.data$Richness, labeller = label_both,
  #                scales = 'free_x', space = 'free') +
  #     theme(strip.text.x = element_text(size=12),
  #           panel.border = element_blank(),
  #           panel.spacing = unit(0.25, "cm"))
  # }
  return(plot)
}

#' @title Split response into contributions by species identities,
#'        species interactions and any additional experimental structures
#'
#' @description
#' Helper function to simplify the creation of the response contributions plot
#' for models fit using the \code{\link[DImodels:DI]{DI()}} function from the
#' \code{\link[DImodels:DImodels-package]{DImodels}} package.
#'
#'
#' @importFrom dplyr mutate %>% group_by distinct across
#'                   all_of slice_sample ungroup near
#' @importFrom tidyr pivot_longer
#' @importFrom forcats fct_rev fct_inorder
#' @importFrom ggplot2 ggplot element_text aes geom_col labs geom_errorbar
#'                     theme theme_bw facet_grid guides guide_legend
#'                     scale_fill_manual unit element_blank label_both
#' @importFrom PieGlyph geom_pie_glyph
#' @importFrom rlang !!! !! sym .data
#' @importFrom stats predict
#' @importFrom cli cli_progress_along
#' @importFrom DImodels DI_data_FG
#'
#' @param model A Diversity Interactions model object fit by using the
#'              \code{\link[DImodels:DI]{DI()}} function from the
#'              \code{\link[DImodels:DImodels-package]{DImodels}} package.
#' @param data A dataframe specifying communities of interest for which user
#'             wants to compare response. If left blank, a random selection
#'             of communities from the original data used to fit the model
#'             would be selected.
#' @param no_random Number of random communities to select from each level of
#'                  richness from original data if communities are not specified.
#' @param FG The functional grouping of the species in the design. Species
#'           belonging to the same functional group will be assigned with
#'           different shades of the same colour. The user can manually specify
#'           a character vector giving the functional group each species
#'           belongs to. If left empty the function will try to get a functional
#'           grouping from the original \code{\link[DImodels:DI]{DI}} model
#'           object.
#' @param exp_str A list specifying values for additional experimental
#'                structures in the model other than the proportions.
#'                This would be useful to compare the predictions across
#'                different values for a categorical variable. One plot will
#'                be generated for each unique combination of values specified
#'                here and they will be arranged in a grid according to the
#'                value specified in `nrow` and `ncol`.
#' @param groups A list specifying groupings to arrange coefficients into.
#'               The coefficients within a group will be added together and
#'               shown as a single entity on the plot. This could be useful for
#'               grouping the species interactions into a single term for
#'               better visibility.
#' @param conf.level The confidence level for calculating confidence or
#'                   prediction intervals.
#' @param se A logical value indicating whether to show prediction intervals
#'           for predictions in the plot.
#' @param nrow Number of rows in which to arrange the final plot
#'             (when `exp_str` is specified).
#' @param ncol Number of columns in which to arrange the final plot
#'             (when `exp_str` is specified).
#'
#' @return A ggmultiplot (ggplot if single plot is returned) class object
#' @export
#'
#' @examples
#' #' ## Load DImodels package to fit the model
#' library(DImodels)
#'
#' ## Load data
#' data(sim2)
#'
#' ## Fit DI model
#' model1 <- DI(prop = 3:6, DImodel = 'FULL', data = sim2, y = 'response')
#'
#' ## Create visualisation
#' ## If no communities are specified a random selection
#' ## of communities from the original data is used
#' response_contributions(model1)
#'
#' ## Three communities are selected by default but user can modify
#' ## Select 5 communities from each level of richness
#' response_contributions(model1, no_random = 4)
#'
#' ## Can also manually specify communities of interest
#' my_comms <- data.frame(p1 = c(1, 0, 0,   0.5, 1/3, 0.25),
#'                        p2 = c(0, 0, 0.5, 0,   1/3, 0.25),
#'                        p3 = c(0, 1, 0.5, 0,   1/3, 0.25),
#'                        p4 = c(0, 0, 0,   0.5, 0,   0.25))
#'
#' response_contributions(model1, data = my_comms)
#'
#' ## Group species by functional groups
#' response_contributions(model1, FG = c("G", "G", "H", "H"))
response_contributions <- function(model, data = NULL, no_random = 3,
                                   exp_str = list(), groups = list(),
                                   conf.level = 0.95,
                                   se = FALSE, FG = NULL,
                                   nrow = 0, ncol = 0){
  # Ensure specified model is fit using the DI function
  if(missing(model) || !inherits(model, "DI")){
    model_not_DI(call_fn = "response_contributions")
  }

  # Get original data used to fit the model
  original_data <- model$original_data

  # Get all species in the model
  model_species <- eval(model$DIcall$prop)
  # If species were specified as column indices extract names
  if(is.numeric(model_species)){
    model_species <- colnames(original_data)[model_species]
  }

  # If data is missing then use the original data used to fit the model
  # and choose a random subset of the communities
  if(is.null(data)){
    data <- original_data %>% distinct(across(all_of(model_species)), .keep_all = T) %>%
      add_ID_terms(model) %>%
      mutate(.Richness = get_richness(.data, model_species)) %>%
      group_by(.data$.Richness) %>%
      slice_sample(n = no_random, replace = TRUE) %>%
      ungroup() %>%
      distinct(across(all_of(model_species)), .keep_all = T) %>%
      mutate(.Community = as.character(1:length(.data$.Richness))) %>%
      as.data.frame()

  } else {

    sanity_checks(data = data, prop = model_species)
    # if(!inherits(data, "data.frame")){
    #   stop("'data' should be a dataframe specifying the proportions of different species in data")
    # }
    # if(! all(model_species %in% colnames(data))){
    #   stop(glue::glue("The column names of the proportions of species in 'data' should be same as those specified when fitting the model.
    #                   Rename the column names to {paste0(model_species, collapse = ',')}"))
    # }
    # if(! all(near(rowSums(data[, model_species]), 1, tol = (.Machine$double.eps)^0.25))){
    #   stop("The species proportions should sum to 1")
    # }

    data <- data %>%
      add_ID_terms(model) %>%
      mutate(.Richness = get_richness(.data, model_species),
             .Community = as.character(1:length(.data$.Richness)))
  }

  if(is.null(FG)){
    FG <- eval(model$DIcall$FG)
  }

  # model_coeffs <- model$coefficients
  # iden_coeffs <- model_coeffs[grepl(x = names(model_coeffs),
  #                                   pattern = paste0('^', model_species, '$', collapse = '|'))]
  #
  # # Predictions
  # data$Pred <- predict(model, data)
  #
  # # Hadamard product of species matrix with identity coeff vector
  # data[, paste(model_species, 'cont', sep = '_')] <- data[, model_species] * tcrossprod(rep.int(1L, nrow(data[, model_species])), iden_coeffs)
  #
  # data$Int <- data$Pred - rowSums(data[, paste(model_species, 'cont', sep = '_')])
  #
  # plot_data <- data %>%
  #   pivot_longer(cols = all_of(c(paste(model_species, 'cont', sep = '_'), 'Int')),
  #                names_to = 'Contributions', values_to = 'Value')
  #

  # Prepare data for plotting
  # Add interaction terms
  old <- ncol(data) # To find number of colours for interaction terms
  data <- add_interaction_terms(model = model, data = data)
  new <- ncol(data)
  # Add any experimental structures specified by user
  # Ensure experimental structure are specified correctly
  exp_str <- check_exp_str(model = model, exp_str = exp_str)

  if(length(exp_str) > 0){
    data <- add_exp_str(data = data,
                        exp_str = exp_str)
  }

  # Split the predicted response into respective contributions
  plot_data <- response_contributions_data(data = data, model = model,
                                           groups = groups,
                                           conf.level = conf.level)

  # Colours for species

  ## Colours for ID effects
  mod_ids <- if (is.null(model$DIcall$ID)) model_species else eval(model$DIcall$ID)
  ID_cols <- get_colours(vars = mod_ids, FG = FG)

  ## Colours for interaction effects
  ## Number of interaction terms
  DImodel_tag <- eval(model$DIcall$DImodel)
  if (is.null(DImodel_tag)) {
    DImodel_tag <- "CUSTOM"
  }

  if(DImodel_tag == "FG"){
    n_ints <- ncol(DI_data_FG(prop = model_species, data = data,
                              FG = eval(model$DIcall$FG))$FG)
  } else {
    n_ints <- new - old
  }
  if(n_ints == 0){
    int_cols <-  NULL
  } else {
    int_cols <- get_shades(shades = n_ints)[[1]]
  }

  ## Colours for experimental structures
  ## Number of experimental structures
  coeffs <- coef(model)
  n_coeff <- ifelse(is.na(coeffs["theta"]), length(coeffs), length(coeffs) - 1)
  n_exp <- n_coeff - length(ID_cols) - n_ints
  if(n_exp == 0){
    exp_cols <- NULL
  } else {
    exp_cols <- get_shades(colours = "red", shades = n_exp)[[1]]
  }

  # Final colours
  colours <- c(ID_cols, int_cols, exp_cols)

  if(length(exp_str) > 0){
    ids <- unique(plot_data$.add_str_ID)
    plots <- lapply(cli_progress_along(1:length(ids), name = "Creating plot",
                                       format = paste0(
                                         "{pb_spin} Creating plot ",
                                         "[{pb_current}/{pb_total}]   ETA:{pb_eta}"
                                       )),
                    function(i){
                      data <- plot_data %>% filter(.data$.add_str_ID == ids[i])
                      plot <- response_contributions_plot(data = data,
                                                          colours = colours,
                                                          se = se) +
                                  labs(subtitle = ids[i])
                    })
    if(length(plots) > 1){
      plot <- new("ggmultiplot", plots = plots, nrow = nrow, ncol = ncol)
    } else {
      plot <- plots[[1]]
    }
    cli::cli_alert_success("Created all plots.")
  } else {
    plot <- response_contributions_plot(data = plot_data, colours = colours,
                                        se = se)
    cli::cli_alert_success("Created plot.")
  }
  return(plot)
}
