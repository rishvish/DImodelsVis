#' @title Split response into contributions by species identities and interactions
#'
#' @importFrom dplyr mutate %>% group_by distinct across all_of slice_sample ungroup near
#' @importFrom tidyr pivot_longer
#' @importFrom forcats fct_rev fct_inorder
#' @importFrom ggplot2 ggplot element_text aes geom_col labs theme theme_bw facet_grid guides guide_legend scale_fill_manual unit element_blank label_both
#' @importFrom PieGlyph geom_pie_glyph
#' @importFrom rlang !!! !! sym .data
#' @importFrom stats predict
#'
#' @param model A Diversity Interactions model object fit by using the \code{\link[DImodels:DI]{DI()}} function from the \code{\link[DImodels:DImodels-package]{DImodels}} package.
#' @param communities A dataframe specifying communities of interest for which user wants to compare response. If left blank, a random selection of communities from the original data used to fit the model would be selected
#' @param no_random Number of random communities to select from each level of richness from original data if communities are not specified
#' @param FG The functional grouping of the species in the design. Species belonging to the same functional group will be assigned with different shades of the same colour. The user can manually specify a character vector giving the functional group each species belongs to. If left empty the functional will try to get a functional grouping from the original \code{\link[DImodels:DI]{DI}} model object.
#'
#' @return A ggplot object
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
#' response_contributions(model1, communities = my_comms)
#'
#' ## Group species by functional groups
#' response_contributions(model1, FG = c("G", "G", "H", "H"))
response_contributions <- function(model, communities = NULL, no_random = 3, FG = NULL){
  # Sanity checks
  # Ensure model is a DImodels object
  if(!inherits(model, 'DI')){
    stop("Model should be a DImodels object")
  }

  # Get original data used to fit the model
  original_data <- model$original_data

  # Get all species in the model
  model_species <- eval(model$DIcall$prop)
  # If species were specified as column indices extract names
  if(is.numeric(model_species)){
    model_species <- colnames(original_data)[model_species]
  }

  if(is.null(communities)){
    communities <- original_data %>% distinct(across(all_of(model_species)), .keep_all = T) %>%
      mutate(Richness = get_richness(.data, model_species)) %>%
      group_by(.data$Richness) %>%
      slice_sample(n = no_random, replace = TRUE) %>%
      ungroup() %>%
      distinct(across(all_of(model_species)), .keep_all = T) %>%
      mutate(community = as.character(1:length(.data$Richness))) %>%
      as.data.frame()

  } else {

    if(!inherits(communities, "data.frame")){
      stop("'communities' should be a dataframe specifying the proportions of different species in communities")
    }
    if(! all(model_species %in% colnames(communities))){
      stop(glue::glue("The column names of the proportions of species in 'communities' should be same as those specified when fitting the model.
                      Rename the column names to {paste0(model_species, collapse = ',')}"))
    }
    if(! all(near(rowSums(communities[, model_species]), 1, tol = (.Machine$double.eps)^0.25))){
      stop("The species proportions should sum to 1")
    }

    communities <- communities %>%
      mutate(Richness = get_richness(.data, model_species),
             community = as.character(1:length(.data$Richness)))
  }

  if(is.null(FG)){
    FG <- eval(model$DIcall$FG)
  }

  model_coeffs <- model$coefficients
  iden_coeffs <- model_coeffs[grepl(x = names(model_coeffs),
        pattern = paste0('^', model_species, '$', collapse = '|'))]

  # Predictions
  communities$Pred <- predict(model, communities)

  # Hadamard product of species matrix with identity coeff vector
  communities[, paste(model_species, 'cont', sep = '_')] <- communities[, model_species] * tcrossprod(rep.int(1L, nrow(communities[, model_species])), iden_coeffs)

  communities$Int <- communities$Pred - rowSums(communities[, paste(model_species, 'cont', sep = '_')])

  plot_data <- communities %>%
    pivot_longer(cols = all_of(c(paste(model_species, 'cont', sep = '_'), 'Int')),
                 names_to = 'Contributions', values_to = 'Value')

  # Colours for species
  colours <- get_colours(species = model_species, FG = FG)

  ggplot(plot_data, aes(x = .data$community, y = .data$Value, fill = fct_rev(fct_inorder(.data$Contributions))))+
    geom_col()+
    facet_grid(~.data$Richness, labeller = label_both,
               scales = 'free_x', space = 'free') +
    theme_bw()+
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          legend.position="top",
          legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          strip.text.x = element_text(size=12),
          panel.border = element_blank(),
          panel.spacing = unit(0.25, "cm"))+
    scale_fill_manual(values = c('#808080', rev(colours)),
                      labels = c('Interactions', rev(model_species)))+
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(y = 'Predicted Response',
         x = 'Community',
         fill = 'Contributions')
}
