#' @usage NULL
NULL
create_conditional_data <- function(model, conditional=NULL, values=c(0.2,0.5,0.8),
                                    resolution = 3){

  #Sanity Checks
  if (!inherits(model, 'DI')){
    stop('Specify a DImodels object in model')
  }

  if(is.null(conditional)){
    conditional <- ''
  }

  if (!inherits(conditional, 'character')){
    stop('Specify the names of the conditional species as a string vector')
  }

  if(!all(between(values, 0 , 1))){
    stop(glue::glue('Values should be a numeric vector with values between 0 and 1'))
  }

  # Get data used to fit the model
  og_data <- model$original_data

  # Get all species in the model
  species <- eval(model$DIcall$prop)

  if(is.numeric(species)){
    species <- colnames(og_data)[species]
  }

  if(length(species) < 3){
    stop(glue::glue("Can't create ternary diagram for models with less than 3 species"))
  }

  if (all(conditional !='') & !all(conditional %in% species)){
    stop('The conditioning species should be present in the model.
         \nPlease check if you have entered the name correctly.')
  }

  # Get species to be shown in the ternary diagram
  tern_species <- species[species!=conditional]

  if (length(tern_species)!=3){
    if (length(tern_species)<3){
      stop(glue::glue("Can't create ternary diagram with less than {length(tern_species)} species.
                Remove any {3 - length(tern_species)} species from the conditional parameter"))
    } else if (length(tern_species>3)){
      warning(glue::glue("More than three species for ternary diagram;
Creating the diagram for species `{tern_species[1]}`, `{tern_species[2]}`, and `{tern_species[3]}` and assuming {paste0(tern_species[4:length(tern_species)], collapse = ',')} to be 0, but it might not be very informative.
Try increasing the number of conditioning species or grouping species into functional groups"))
    }
  }

  # Ensure resolution isn't to high or low
  if(!(is.numeric(resolution) & length(resolution) == 1 & between(resolution, 1, 5))){
    warning(glue::glue('Resolution should be a numeric value between 1 and 10, reverting back to the default value of 3'))
  }

  # Prepare data for creating conditional proportions
  base <- seq(0,1,l=100*2*resolution)
  high <- seq(0,sin(pi/3),l=87*2*resolution)
  triangle <- expand.grid(base = base, high = high)
  triangle <- subset(triangle, (((base*sin(pi/3)*2) > high) & (((1-base)*sin(pi/3)*2) > high)))

  # Extrapolate 2-d simplex coordinates to represent proportions
  # of the three species shown in the simplex
  triangle <- triangle %>%
                mutate(!! tern_species[1] := .data$high*2/sqrt(3),
                       !! tern_species[3] := .data$base - .data$high/sqrt(3),
                       !! tern_species[2] := 1 - .data$high*2/sqrt(3) - (.data$base - .data$high/sqrt(3)))

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
    cond_data <- add_predictions(data = triangle,
                                 model = model,
                                 var = 'Pred')
  } else {
    # Rescaling species to be shown in ternary for each species
    # in the conditional parameter according to the values specified
    # in the values parameter
    cond_data <- lapply(values, function(x){
      # Scale proportion of species within the ternary
      scaled_data <- triangle %>%
                        mutate(!! tern_species[1] := rescale(!!sym(tern_species[1]), min = 0, max = 1-x),
                               !! tern_species[2] := rescale(!!sym(tern_species[2]), min = 0, max = 1-x),
                               !! tern_species[3] := rescale(!!sym(tern_species[3]), min = 0, max = 1-x))
      # Add remaining species in the data
      remaining_species <- matrix(0, ncol=length(conditional), nrow= nrow(triangle))
      colnames(remaining_species) <- conditional
      scaled_data <- cbind(scaled_data, remaining_species)

      # Update value of species which we are conditioning on
      species_data <- lapply(conditional, function(cond_sp){
        # Add identifier for grouping data for a particular species
        sp_data <- scaled_data %>%
                      mutate(Sp = cond_sp,
                             Value = paste0(cond_sp, ' = ', x))

        # To avoid any rounding issues and ensure all species prportions sum to 1
        sp_data[, cond_sp] <- 1 - rowSums(sp_data[, tern_species])
        sp_data
        }) %>% bind_rows()
    }) %>% bind_rows()

    # Predicting the response for the communities
    cond_data <- add_predictions(data = cond_data,
                                 model = model,
                                 var = 'Pred')
  }
  return(cond_data)
}


#' @title Conditional ternary diagrams
#'
#' @importFrom metR geom_text_contour
#' @importFrom grDevices terrain.colors
#' @importFrom dplyr tibble between
#' @importFrom ggplot2 geom_raster scale_fill_stepsn geom_contour theme_void guide_colorbar coord_fixed
#'
#' @param model A Diversity Interactions model object fit by using the \code{\link[DImodels:DI]{DI()}} function from the \code{\link[DImodels:DImodels-package]{DImodels}} package.
#' @param conditional A character string describing the species to condition the ternary slices on.
#' @param values A vector of numbers between 0 and 1 describing the values at which to takes slices of the ternary
#' @param resolution A number between 1 and 5 describing the resolution of the resultant graph. A high value would result in a higher definition figure but at the cost of being computationlly expernsive.
#' @param nlevels The number of levels for the contour map
#' @param colours The colours for the contour map. The default colours scheme is the terrain.colors from hcl.pals()
#' @param lower_lim A number to set a custom lower limit for the contour. The default is minimum of the prediction
#' @param upper_lim A number to set a custom upper limit for the contour. The default is maximum of the prediction
#' @param .contour_text A boolean value indicating whether to include labels on the contour lines showing their values
#' @param .axis_guides A boolean value indicating whether to show axis guides in the ternary
#' @param axis_label_size A numeric value to adjust the size of the axis labels in the contour
#' @param vertex_label_size A numeric value to adjust the size of the vertex labels in the contour
#'
#' @return A ggplot object
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
#'                     resolution = 1)
conditional_ternary <- function(model, conditional = NULL,
                                values = c(0.2, 0.5, 0.8),
                                resolution = 3,
                                nlevels = 7,
                                colours = NULL,
                                lower_lim = NULL,
                                upper_lim = NULL,
                                .contour_text = TRUE,
                                .axis_guides = FALSE,
                                axis_label_size = 4,
                                vertex_label_size = 5
                                ){
  # Create data in appropriate format for plotting
  plot_data <- create_conditional_data(model = model, conditional = conditional,
                                       values = values, resolution = resolution)

  # Create colour-scale (legend) for plot
  # Lower limit of legend
  if(!is.null(lower_lim)){
    # If user specified lower limit ensure it is numeric and of length 1
    if(!is.numeric(lower_lim) | length(lower_lim) != 1){
      stop('lower_lim should be a single numeric value')
    }
  } else {
    # If user didn't specify lower limit assume it to be min of predicted response
    lower_lim <- round(min(plot_data$Pred), 2)
  }

  # Upper limit of legend
  if(!is.null(upper_lim)){
    # If user specified upper limit ensure it is numeric and of length 1
    if(!is.numeric(upper_lim) | length(upper_lim) != 1){
      stop('upper_lim should be a single numeric value')
    }
  } else {
    # If user didn't specify upper limit assume it to be max of predicted response
    upper_lim <- round(max(plot_data$Pred),2)
  }

  # Create breaks between range of legend
  if(!is.numeric(nlevels) | length(nlevels) != 1){
    stop('nlevels should be a single numeric value')
  }

  size <- nlevels + 1
  breaks <- round(seq(lower_lim, upper_lim, length.out= size), 2)

  # Choose colours for contour
  if(!is.null(colours)){
    # If user has specified colours ensure they have same length as nlevels
    if(length(colours)!= nlevels){
      stop(glue::glue('Colours should be specified as a vector of character strings with same length as value of nlevels ({nlevels} in this case)'))
    }

    # If user has specified colours ensure they are valid
    cols <- areColours(colours)
    if(!all(cols)){
      stop(glue::glue('{paste(names(cols)[which(cols == F)], collapse = ", ")} are not valid colours'))
    }
  } else {
    # If user didn't specify colours then use the default terrain colours
    colours <- terrain.colors(nlevels, rev = T)
  }

  # Labels for the ternary
  tern_labels <- colnames(plot_data)[3:5]

  # Labels for the ternary axes
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
                        Value = paste0(sp, ' = ', x),
                        Pred = 0)
    })
  }) %>% bind_rows()

  # Ensure vertex and axis label sizes are numeric
  if(!is.numeric(axis_label_size) | length(axis_label_size) != 1){
    stop('axis_label_size should be a single numeric value')
  }

  if(!is.numeric(vertex_label_size) | length(vertex_label_size) != 1){
    stop('vertex_label_size should be a single numeric value')
  }

  # Create plot
  pl <- ggplot(plot_data, aes(x = .data$base, y = .data$high,
                              z = .data$Pred))+
    geom_raster(aes(fill = .data$Pred))+
    scale_fill_stepsn(colours = colours, breaks = breaks,
                      labels = function(x){
                        x
                      },
                      limits = c(lower_lim, upper_lim),
                      show.limits = T)+
    geom_contour(breaks = breaks, colour = 'black')+
    geom_text(data = data.frame(x = c(0.5, 1, 0), y = c(sqrt(3)/2, 0,  0),
                                label = tern_labels,
                                Pred = 0),
              aes(x= .data$x, y= .data$y, label = .data$label),
              size = vertex_label_size, fontface='plain',
              nudge_x = c(0, 0.05, -0.05),
              nudge_y = c(0.05, 0, 0))+
    geom_segment(data = data.frame(x = c(0, 0, 1), y = c(0,0,0),
                                   xend = c(1, 0.5, 0.5),
                                   yend = c(0, sqrt(3)/2, sqrt(3)/2),
                                   Pred = 0),
                 aes(x=.data$x, y=.data$y, xend=.data$xend, yend=.data$yend), linewidth = 1)+
    geom_text(data = axis_labels,
              aes(x=.data$x1, y=.data$y1, label=.data$label), nudge_y=-0.055, size = axis_label_size)+
    geom_text(data = axis_labels,
              aes(x=.data$x2, y=.data$y2, label=.data$rev_label),  nudge_x=-0.055, nudge_y=0.055, size = axis_label_size)+
    geom_text(data = axis_labels,
              aes(x=.data$x3, y=.data$y3, label=.data$rev_label),  nudge_x=0.055, nudge_y=0.055, size = axis_label_size)+
    theme_void()+
    guides(fill = guide_colorbar(frame.colour = 'black',
                                 ticks.colour = 'black',
                                 title = 'Predicted\nResponse',
                                 show.limits = T))+
    facet_wrap(~ Value, ncol = length(values))+
    coord_fixed()+
    theme(legend.key.size = unit(0.1, 'npc'),
          legend.key.height = unit(0.04, 'npc'),
          legend.title = element_text(size = 14, vjust = 0.75),
          plot.subtitle = element_text(hjust=0.5, size=14),
          strip.text = element_text(size =14, vjust = 0.5),
          legend.text = element_text(size = 12, angle = 45, vjust = 1.2, hjust = 1.2),
          legend.position = 'bottom')



  if(.axis_guides){
    pl <- pl +
      geom_segment(data = axis_labels,
                   aes(x = .data$x1, y = .data$y1,
                       xend = .data$x2, yend = .data$y2), colour='grey',
                   linetype='dashed', linewidth=1, alpha = .75)+
      geom_segment(data = axis_labels,
                   aes(x = .data$x1, y = .data$y1,
                       xend = .data$x3, yend = .data$y3), colour='grey',
                   linetype='dashed', linewidth=1, alpha = .75)+
      geom_segment(data = axis_labels,
                   aes(x = .data$x2, y = .data$y2,
                       xend = rev(.data$x3), yend = rev(.data$y3)), colour='grey',
                   linetype='dashed', linewidth=1, alpha = .75)
  }


  if(.contour_text){
    pl <- pl +
      geom_text_contour(skip=0, breaks = breaks,
                        label.placer = metR::label_placer_fraction(0.15),
                        size=3.5, nudge_x = 0.015, nudge_y = 0.015)

  }

  # if(.pies){
  #   pl <- pl +
  #     # scale_fill_manual(values = colours,
  #     #                   guide = guide_legend(order = 2))+
  #     ggnewscale::new_scale_fill()+
  #     PieGlyph::geom_pie_glyph(data = data.frame(Sp = conditional,
  #                                                values = values,
  #                                                Value = apply(arrange(expand.grid(conditional, values), Var1), 1, paste, collapse=" = "),
  #                                                Others = 1 - values,
  #                                                Pred = 0),
  #                              slices = c('values', 'Others'),
  #                              aes(x = 0.5, y = 1))+
  #     scale_fill_manual(values = c('green', 'grey'),
  #                       guide = guide_legend(order = 2))
  #
  #     # guides(# fill = guide_colorbar(frame.colour = 'black',
  #     #         #                     ticks.colour = 'black',
  #     #         #                     title = 'Predicted\nResponse'),
  #     #        fill = guide_legend(values = c('green', 'grey')))
  #
  # } else {
  #   pl <- pl +
  #     guides(fill = guide_colorbar(frame.colour = 'black',
  #                                  ticks.colour = 'black',
  #                                  title = 'Predicted\nResponse'))
  # }






  return(pl)
}

