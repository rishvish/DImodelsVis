#' @title Visualising model selection
#' @description This function helps to visualise model selection by showing a visual comparison between models using information criteria
#'
#' @importFrom ggplot2 geom_line geom_point geom_label aes scale_colour_manual geom_hline geom_segment arrow
#' @importFrom dplyr %>% mutate rename arrange
#'
#' @param models List of Diversity Interaction (DI) models fit by using the \code{\link[DImodels:DI]{DI()}} function from the \code{\link[DImodels:DImodels-package]{DImodels}} package.
#' @param metric Metric used for comparisons between models. Takes values from c("AIC", "BIC", "AICc", "BICc", "logLik"), Can choose a single or multiple metrics for comparing the different models
#' @param sort A boolean value indicating whether to sort the model from highest to lowest value of chosen metric
#' @param model_names A character string describing the names to display on X-axis for each model in order they appear in the models parameter
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' library(DImodels)
#'
#' ## Load data
#' data(sim2)
#'
#' ## Fit different DI models
#' mod_AV <- DI(prop = 3:6, DImodel = "AV", data = sim2, y = "response")
#' mod_FULL <- DI(prop = 3:6, DImodel = "FULL", data = sim2, y = "response")
#' mod_FG <- DI(prop = 3:6, DImodel = "FG", FG = c("G","G","L","L"),
#'              data = sim2, y = "response")
#' mod_AV_theta <- DI(prop = 3:6, DImodel = "AV", data = sim2,
#'                    y = "response", estimate_theta = TRUE)
#' mod_FULL_theta <- DI(prop = 3:6, DImodel = "FULL", data = sim2,
#'                      y = "response", estimate_theta = TRUE)
#' mod_FG_theta <- DI(prop = 3:6, DImodel = "FG", FG = c("G","G","L","L"),
#'                    data = sim2, y = "response", estimate_theta = TRUE)
#'
#' models_list <- list(mod_AV, mod_FULL, mod_FG,
#'                     mod_AV_theta, mod_FULL_theta, mod_FG_theta)
#'
#' ## Compare different information criteria for the different models
#' model_selection(models = models_list,
#'                 metric = c("AIC", "AICc"))
#'
#' ## If single metric is specified then breakup of metric
#' ## between loglik and penalty will be shown
#' model_selection(models = models_list,
#'                 metric = c("AIC"))
#'
#' ## Sort models
#' model_selection(models = models_list,
#'                 metric = c("AIC"), sort = TRUE)
#'
#' ## If multiple metrics are specified then sorting will be done on first metric
#' model_selection(models = models_list,
#'                 metric = c("AIC", "logLik", "BIC"), sort = TRUE)
#'
#' ## By default the labels on the X-axis for the models will be
#' ## created by assigning a unique ID to each model sequentially
#' ## in the order they appear in the models object
#' model_selection(models = models_list,
#'                 metric = c("AIC", "logLik", "BIC"), sort = TRUE)
#'
#' ## When possible the variables names of objects contatining the
#' ## individual models would be used as axis labels
#' model_selection(models = list(mod_AV, mod_FULL, mod_FG,
#'                               mod_AV_theta, mod_FULL_theta, mod_FG_theta),
#'                 metric = c("AIC", "logLik", "BIC"), sort = TRUE)
#'
#' ## If neither of these two situations are desirable custom labels
#' ## for each model can be specified using the model_names parameter
#' model_selection(models = list(mod_AV, mod_FULL, mod_FG,
#'                               mod_AV_theta, mod_FULL_theta, mod_FG_theta),
#'                 metric = c("AIC", "logLik", "BIC"), sort = TRUE,
#'                 model_names = c("AV model", "Full model", "FG model",
#'                                 "AV theta", "Full theta", "FG theta"))
model_selection <- function(models, metric = c("AIC", "BIC", "AICc", "BICc", "logLik"), sort = FALSE,
                            model_names = names(models)){
  # Ensure all objects in model list are DI models
  if(!all(sapply(models, function(x) {inherits(x, "DI")}))){
    stop("Models should be a list of DImodels objects")
  }

  # Ensure model_names is same length as models list
  if(!is.null(model_names) & (length(model_names) != length(models))){
    stop("'model_names' should be a character string with same length as models list.")
  }

  # If the user has not named the list of models, give them same names as the variable names
  if(is.null(model_names)){
    p_tree <- substitute(models)
    # If possible then use variables names as X-axis labels
    if(length(p_tree) > 1) {
      message("The models are given the same names as the variables they were stored in. If this is not desirable consider providing names for the models in the models_name 'parameter'.")
      names(models) <- sapply(p_tree[-1], deparse)
    # Otherwise give a unique ID to models sequentially
    } else {
      message("The models are given numeric identifiers in the order they appear in the models parameter. If this is not desirable consider providing names for the models in the models_name 'parameter'.")
      names(models) <- paste("Model", 1:length(models))
    }
  } else {
    names(models) <- model_names
  }


  # Ensure all metrics specified by user match the choices
  # Manual check because partial mapping not working due to choices being too similar
  if(!all(metric %in% c("AIC", "BIC", "AICc", "BICc", "logLik"))){
    stop("'metric' should be one of c(\"AIC\", \"BIC\", \"AICc\", \"BICc\", \"logLik\")")
  }

  # If user didn't choose anything then AIC is default
  if(missing(metric)){
    metric <- "AIC"
  }
  metric <- match.arg(metric, choices = c("AIC", "BIC", "AICc", "BICc", "logLik"),
                      several.ok = TRUE)

  # Create dataframe for plotting
  plot_data <- data.frame(model_name = names(models))
  # Add the relevant metrics
  plot_data <- plot_data %>%
                  mutate("AIC" = AIC_vec(models) %>% round(2),
                         "BIC" = BIC_vec(models) %>% round(2),
                         "AICc" = AICc_vec(models) %>% round(2),
                         "BICc" = BICc_vec(models) %>% round(2),
                         "logLik" = -2*logLik_vec(models) %>% round(2))
  # Sort models
  if(sort){
    plot_data <- plot_data %>%
      arrange(desc(!!sym(metric[1]))) %>%
      mutate(model_name = forcats::fct_inorder(.data$model_name))
  }

  # Convert logLik to -2*logLik
  if (any("logLik" %in% metric)){
    metric[metric == "logLik"] <- "-2 * logLik"
    plot_data <- plot_data %>% rename("-2 * logLik" = logLik)
  }

  # Show multiple metrics
  if(length(metric) > 1){
    plot_data <- plot_data %>%
      pivot_longer(cols = all_of(metric), values_to = "Value", names_to = "Metric") %>%
      mutate('Metric' = forcats::fct_inorder(.data$Metric))

    pl <- ggplot(plot_data, aes(x = .data$model_name, y = .data$Value,
                                colour = .data$Metric, group = .data$Metric))+
      geom_line(linewidth = 1.5)+
      geom_point(size = 3)+
      geom_label(aes(label = .data$Value), colour = 'black', nudge_y = 1.5)+
      labs(x = "Model")+
      theme_bw()+
      theme(legend.position = 'top')+
      scale_colour_manual(values = colour_blind_friendly_cols(length(metric)))
    # Show a single metric
  } else {
    # Dashed lines for rule of 2
    rule_of_2 <- min(plot_data[, metric]) + c(-2, 2)

    # If logLik is metric there will be no penalty
    if(metric == "-2 * logLik"){
      pl <- ggplot(plot_data, aes(x = .data$model_name, y = !!sym(metric), group = 1))+
        geom_line(colour = "#CC79A7", linewidth = 1.5)+
        geom_point(colour = "black",
                   fill = "#CC79A7", size = 5, pch = 21, stroke = 1)+
        geom_hline(yintercept = rule_of_2, linetype = 3, linewidth = 0.75)+
        geom_label(aes(label = !!sym(metric)), nudge_y = 1.5)+
        labs(x = "Model", y = "-2 * LogLik")+
        theme_bw()
    # For all other metrics there will be a penalty
    } else {
      plot_data <- plot_data %>%
        mutate("Goodness of fit" = -2*logLik_vec(models) %>% round(2),
               "Penalty" = (!!sym(metric) - .data$`Goodness of fit`) %>% round(2))

      label_data <- plot_data %>%
        pivot_longer(cols = all_of(c("Goodness of fit", "Penalty", metric))) %>%
        group_by(.data$model_name) %>%
        mutate(position = ifelse(.data$name == "Penalty", cumsum(.data$value) - .data$value/2, .data$value)) %>%
        ungroup()

      pl <- ggplot(plot_data, aes(x = .data$model_name, group = 1))+
        geom_line(aes(y = .data$`Goodness of fit`), colour = "#f59424", linewidth = 1.5)+
        geom_line(aes(y = !!sym(metric)), colour = "#CC79A7", linewidth = 1.5)+
        geom_point(aes(y = .data$`Goodness of fit`), colour = "black",
                   fill = "#f59424", size = 5, pch = 21, stroke = 1)+
        geom_point(aes(y = !!sym(metric)), colour = "black",
                   fill = "#CC79A7", size = 5, pch = 21, stroke = 1)+
        geom_segment(aes(xend = .data$model_name, y = .data$`Goodness of fit` + 1, yend = !!sym(metric) - 1),
                     arrow = arrow(ends = "both"), linewidth = 1)+
        geom_hline(yintercept = rule_of_2, linetype = 3, linewidth = 0.75)+
        geom_label(data = label_data,
                   aes(y = .data$position, label = .data$value, fill = forcats::fct_inorder(.data$name)),
                   nudge_y = c(-1.75, 0, 1.75),
                   nudge_x = c(0, 0.1, 0))+
        labs(x = "Model", y = "Value", fill = "Measure")+
        theme_bw()+
        theme(legend.position = 'top')+
        scale_fill_manual(values = c("#f59424", "#0a77bd", "#CC79A7"))
    }
  }
  return(pl)
}
