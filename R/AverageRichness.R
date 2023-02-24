#' @title Visualise average response at each level of richness
#'
#' @importFrom dplyr mutate %>% group_by distinct across all_of
#' @importFrom ggplot2 ggplot geom_line aes geom_point position_dodge position_identity labs theme theme_bw
#' @importFrom PieGlyph geom_pie_glyph
#' @importFrom rlang !!! !! sym .data
#' @importFrom DImodels DI_data_E_AV
#' @importFrom stats predict
#' @importFrom modelr add_predictions
#'
#' @param model A Diversity Interactions model object fit by using the \code{\link[DImodels:DI]{DI()}} function from the \code{\link[DImodels:DImodels-package]{DImodels}} package.
#' @param gradient Diversity gradient to show on the X-axis, one of 'richness' or 'evenness'. Defaults to 'richness'.
#' @param communities Communities which are used to calculate the average response. Accepts one of 'original' or 'equi'. 'original' (the default) calculates the average using the original communities from the raw data used to fit the \code{\link[DImodels:DI]{DI}} model while 'equi' uses all possible equi-proportional communities at each level of richness to calculate the average.
#' @param threshold For high levels of richness (>20) it might not be feasible to compute and store all equi-proportional communities. This parameter defines the number of communities randomly selected at each level of richness.
#'
#' @return A ggplot object
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
#' average_richness(model = mod)
#'
#' ## Can also calculate average response using all possible
#' ## equi-proportional communities at each level of richness
#' average_richness(model = mod, communities = 'equi')
#'
#' ## Can also plot across evenness
#' data(sim5)
#' mod1 <- DI(prop = paste0('p', 1:9), DImodel = 'AV',
#'           data = sim5, y = 'response', estimate_theta = TRUE)
#' average_richness(model = mod1, gradient = 'evenness')
average_richness <- function(model, gradient = c('richness', 'evenness'),
                             communities = c('original', 'equi'),
                             threshold = NULL){
  # Sanity checks
  # Ensure model is a DImodels object
  if(!inherits(model, 'DI')){
    stop("Model should be a DImodels object")
  }

  # Ensure no_random is numeric
  if(!is.null(threshold) && (length(threshold) != 1 || !is.numeric(threshold))){
    stop("'threshold' should be number describing the number of random communities to select at each level of richness.")
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


  if(communities == 'equi'){
    communities <- get_equi_comms(length(model_species),
                                  model_species,
                                  threshold)
  } else {
    communities <- original_data
  }

  # Create data for plotting average line
  plot_data <- communities %>%
    mutate('richness' = get_richness(.data, model_species),
           'evenness' = get_evenness(.data, model_species)) %>%
    add_predictions(model = model, var = "Pred") %>%
    group_by(!! sym(gradient)) %>%
    mutate('Avg' = mean(.data$Pred))

  # Create data from plotting pie-charts
  pie_data <- original_data %>% distinct(across(all_of(model_species)), .keep_all = T) %>%
    mutate('richness' = get_richness(.data, model_species),
           'evenness' = get_evenness(.data, model_species)) %>%
    add_predictions(model = model, var = "Pred") %>%
    group_by(!! sym(gradient)) %>%
    mutate('Avg' = mean(.data$Pred))

  # Create plot
  ggplot(data = plot_data, aes(x = !! sym(gradient), y = .data$Pred))+
    geom_line(aes(y = .data$Avg), linewidth = 1)+
    geom_pie_glyph(aes(group = .data$evenness),
                   slices = model_species,
                   colour = 'black', data = pie_data,
                   position = if(gradient == 'richness') position_dodge(0.5) else position_identity())+
    labs(x = gradient,
         y = 'Predicted Response',
         fill = 'Species')+
    theme_bw()+
    theme(legend.position = 'top')
}
