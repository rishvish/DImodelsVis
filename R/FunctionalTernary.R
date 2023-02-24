#' @usage NULL
NULL
get_FG_values <- function(FG){
  nFG <- length(FG)
  counts <- table(FG)
  values <- setNames(vector(mode = 'numeric', length = nFG), FG)
  for(prop in 1:nFG){
    values[prop] <- 1/counts[names(counts) == FG[prop]]
  }
  return(values)
}

#' @usage NULL
NULL
check_FG_values <- function(FG, values){
  groups <- unique(FG)
  sp_in_FG_sums <- sapply(groups, function(group){
    sum(values[FG == group])
  })
  names(sp_in_FG_sums) <- groups
  sum_to_1 <- ifelse(near(sp_in_FG_sums, 1), T, F)
  return(list('sums' = sp_in_FG_sums, 'bool' = sum_to_1))
}

#' @usage NULL
NULL
get_FG_value_mapping <- function(FG, focus = unique(FG)){

  FG_mapping <- lapply(focus, function(group){
    which(FG == group)
  })
  names(FG_mapping) <- focus
  return(FG_mapping)
}

#' @usage NULL
NULL
create_functional_data <- function(model, FG,
                                   values = NULL,
                                   cond_values = NULL,
                                   cond_FG = NULL,
                                   resolution = 3){

  #Sanity Checks
  if (!inherits(model, 'DI')){
    stop('Specify a DImodels object in model')
  }

  if(missing(FG)){
    if(!is.null(eval(model$DIcall$FG))){
      FG <- eval(model$DIcall$FG)
    } else {
      stop('The functional groups should be specified as a character vector of same length as the number of species in the experiment specifying the functional group to which each species belongs.')
    }
  }

  if(!is.character(FG)){
    stop('The functional groups should be specified as a character vector of same length as the number of species in the experiment specifying the functional group to which each species belongs.')
  }

  # Get data used to fit the model
  og_data <- model$original_data

  # Get all species in the model
  species <- eval(model$DIcall$prop)

  if(is.numeric(species)){
    species <- colnames(og_data)[species]
  }

  if(length(species) != length(FG)){
    stop(glue::glue('The functional groups vector should have the same length as the number of species in the experiment ({length(species)} in this case)'))
  }

  if(is.null(values)){
    values <- get_FG_values(FG)
  } else {
    if(!is.numeric(values)){
      stop(glue::glue('The proportions of the species within a function group should be specfied a numeric vector'))
    }

    if(length(values) != length(FG)){
      stop(glue::glue('The values vector should have same length as the FG vector ({length(FG)} in this case)'))
    }

    sp_props_in_FG <- check_FG_values(FG, values)

    if(!all(sp_props_in_FG$bool)){
      faults <- sp_props_in_FG$sums[sp_props_in_FG$bool == F]
      stop(glue::glue("The species proportions within a functional group should sum to 1. The proportions for, {paste(names(faults), collapse = ', ')} equal {paste(faults, collapse = ', ')} respectively."))
    }
  }

  # Get species to be shown in the ternary diagram
  all_FGs <- unique(FG)

  if (length(all_FGs)<3){
    stop(glue::glue("Can't create ternary diagram with less than {length(all_FGs)} functional groups."))
  }

  # Check for conditional FG and conditioning values
  if(length(all_FGs) > 3){
    if(!is.null(cond_FG)){
      if(!is.character(cond_FG) | !all(cond_FG %in% all_FGs)){
        stop(glue::glue("cond_FG should be a character vector specifying the functional group to condtion the ternary diagram when there are more than 3 functional groups in the model."))
      }
    } else {
      message(glue::glue('More than three functional groups in model and no function group specified to condition on. Choosing `{all_FGs[1]}`, `{all_FGs[2]}`, and `{all_FGs[3]}` to plot in ternary and conditioning on {paste0(all_FGs[4:length(all_FGs)], collapse = ', ')}'))
      cond_FG <- all_FGs[4:length(all_FGs)]
    }

    if(!is.null(cond_values)){
      if(!is.numeric(cond_values) | !all(between(cond_values, 0 , 1))){
        stop(glue::glue('cond_values should be a numeric vector with values between 0 and 1'))
      }
    } else {
      message(glue::glue('Choosing default values of 0.2, 0.5 and 0.8 for conditioning.'))
      cond_values <- c(0.2, 0.5, 0.8)
    }
  }

  tern_FGs <- setdiff(unique(FG), cond_FG)
  if(length(tern_FGs) > 3){
    warning(glue::glue("More than three functional groups to show on ternary diagram;
Creating the diagram for species `{tern_FGs[1]}`, `{tern_FGs[2]}`, and `{tern_FGs[3]}` and assuming all of {paste0(tern_FGs[4:length(tern_FGs)], collapse = ', ')} to be 0, but it might not be very informative."))
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
    mutate(!! tern_FGs[1] := .data$high*2/sqrt(3),
           !! tern_FGs[3] := .data$base - .data$high/sqrt(3),
           !! tern_FGs[2] := 1 - .data$high*2/sqrt(3) - (.data$base - .data$high/sqrt(3)))

  # If there are only three functional groups in the model
  # our data is ready and we can make predictions from the model
  if (length(all_FGs) == 3){
    triangle[, species] <- 0
    focus <- tern_FGs
    FG_mapping <- get_FG_value_mapping(FG, focus = focus)

    for (FGs in focus){
      idx <- FG_mapping[[FGs]]
      triangle[, species[idx]] <- (triangle[, FGs] %o% values[idx])
    }

    cond_data <- add_predictions(data = triangle,
                                 model = model,
                                 var = 'Pred')

    # Update the conditional FG and condtioning values to be NULL as
    # there is nothing to condition on
    cond_FG <- ''
    cond_values <- 0
  } else {
    # Rescaling species to be shown in ternary for each species
    # in the cond_FG parameter according to the values specified
    # in the cond_values parameter

    cond_data <- lapply(cond_values, function(x){
      # Scale proportion of species within the ternary
      scaled_data <- triangle %>%
        mutate(!! tern_FGs[1] := rescale(!!sym(tern_FGs[1]), min = 0, max = 1-x),
               !! tern_FGs[2] := rescale(!!sym(tern_FGs[2]), min = 0, max = 1-x),
               !! tern_FGs[3] := rescale(!!sym(tern_FGs[3]), min = 0, max = 1-x))

      # Add remaining species in the data
      remaining_FG <- matrix(0, ncol=length(cond_FG), nrow= nrow(triangle))
      colnames(remaining_FG) <- cond_FG
      scaled_data <- cbind(scaled_data, remaining_FG)

      # Update value of species which we are conditioning on
      species_data <- lapply(cond_FG, function(cond_sp){
        # Add identifier for grouping data for a particular species
        sp_data <- scaled_data %>%
          mutate(FG = cond_sp,
                 Value = paste0(cond_sp, ' = ', x))

        # To avoid any rounding issues and ensure all species prportions sum to 1
        sp_data[, cond_sp] <- 1 - rowSums(sp_data[, tern_FGs[1:3]])

        sp_data[, species] <- 0
        focus <- c(tern_FGs[1:3], cond_sp)
        FG_mapping <- get_FG_value_mapping(FG, focus = focus)

        for (FGs in focus){
          idx <- FG_mapping[[FGs]]
          sp_data[, species[idx]] <- (sp_data[, FGs] %o% values[idx])
        }

        sp_data
      }) %>% bind_rows()
    }) %>% bind_rows()

    # Predicting the response for the communities
    cond_data <- add_predictions(data = cond_data,
                                 model = model,
                                 var = 'Pred')
  }
  return(list(data = cond_data, cond_FG = cond_FG, cond_values = cond_values))
}

#' @title Conditional ternary diagrams at functional group level
#'
#' @importFrom stats setNames
#'
#' @param model A Diversity Interactions model object fit by using the \code{\link[DImodels:DI]{DI()}} function from the \code{\link[DImodels:DImodels-package]{DImodels}} package.
#' @param FG A character vector specifying the functional grouping of the species in the design. If left empty the functional will try to get a functional grouping from the original \code{\link[DImodels:DI]{DI}} model object.
#' @param values A vector specifying the proportional split between the species within a functional group. The default is to split the functional group proportional equally between each species.
#' @param cond_values A vector of numbers between 0 and 1 describing the values at which to takes slices of the ternary. Necessary only when there are more than 3 funcitonal groups.
#' @param cond_FG A character vector specifying the functional groups to slice the ternary on if there are more than three functional groups
#' @param resolution A number between 1 and 5 describing the resolution of the resultant graph. A high value would result in a higher definition figure but at the cost of being computationlly expernsive.
#' @param nlevels The number of levels for the contour map.
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
#' data(sim4)
#' m1 <- DI(y = "response", prop = paste0("p", 1:6),
#'          DImodel = "AV", data = sim4)
#'
#' ## If values is not specified then each the proportion is
#' ## split equally between species in an FG
#' functional_ternary(model = m1, FG = c("G","G","L","L","H","H"),
#'                    resolution = 1)
#'
#' ## Can also manually specify the split to see different
#' ## slice of ternary
#' functional_ternary(m1, FG = c("G","G","L","L","H","H"),
#'                    resolution = 1,
#'                    values = c(0, 1, 0.5, 0.5, 0, 1))
#'
#' ## If there are more than three functional groups
#' ## then the ternary can be further sliced on three FGs
#' functional_ternary(m1, FG = c("G","G","L","H","H","M"),
#'                    resolution = 1, cond_FG = "M")
#'
#' ## Can further split functional group proportion for species
#' functional_ternary(m1, FG = c('G","G","L","H","H","M"),
#'                    resolution = 1, cond_FG = "M",
#'                    values = c(0, 1, 1, 0.5, 0.5, 1))
functional_ternary <- function(model,
                               FG,
                               values = NULL,
                               cond_values = NULL,
                               cond_FG = NULL,
                               resolution = 3,
                               nlevels = 7,
                               colours = NULL,
                               lower_lim = NULL,
                               upper_lim = NULL,
                               .contour_text = T,
                               .axis_guides = F,
                               axis_label_size = 4,
                               vertex_label_size = 5
){

  # Create data in appropriate format for plotting
  data_list <- create_functional_data(model = model, FG = FG,
                                      values = values, resolution = resolution,
                                      cond_FG = cond_FG,
                                      cond_values = cond_values)

  plot_data <- data_list$data
  cond_FG <- data_list$cond_FG
  cond_values <- data_list$cond_values

  # Create colour-scale (legend) for plot
  # Lower limit of legend
  if(!is.null(lower_lim)){
    # If user specified lower limit ensure it is numeric and of length 1
    if(!is.numeric(lower_lim) | length(lower_lim) != 1){
      stop('lower_lim should be a single numeric value')
    }
  } else {
    # If user didn't specify lower limit assume it to be min of predicted response
    lower_lim <- round(min(plot_data$Pred),2)
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
  axis_labels <- lapply(cond_FG, function(sp){
    lapply(cond_values, function(x){
      positions <- tibble(x1 = seq(0.2,0.8,0.2),
                          y1 = c(0,0,0,0),
                          x2 = .data$x1/2,
                          y2 = .data$x1*sqrt(3)/2,
                          x3 = (1-.data$x1)*0.5+.data$x1,
                          y3 = sqrt(3)/2-.data$x1*sqrt(3)/2,
                          label = .data$x1*(1-x),
                          rev_label = rev(.data$label),
                          Value = ifelse(sp == '', '', paste0(sp, ' = ', x)),
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
  pl <- ggplot(plot_data, aes(x = .data$base, y = .data$high, z = .data$Pred))+
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
              aes(x = .data$x, y = .data$y, label = .data$label),
              size = vertex_label_size, fontface='plain',
              nudge_x = c(0, 0.05, -0.05),
              nudge_y = c(0.05, 0, 0))+
    geom_segment(data = data.frame(x = c(0, 0, 1), y = c(0,0,0),
                                   xend = c(1, 0.5, 0.5),
                                   yend = c(0, sqrt(3)/2, sqrt(3)/2),
                                   Pred = 0),
                 aes(x=.data$x, y=.data$y,
                     xend=.data$xend, yend=.data$yend), linewidth = 1)+
    geom_text(data = axis_labels,
              aes(x=.data$x1, y=.data$y1, label=.data$label, fontface = 'plain'), nudge_y=-0.055, size = axis_label_size)+
    geom_text(data = axis_labels,
              aes(x=.data$x2, y=.data$y2, label=.data$rev_label, fontface = 'plain'),  nudge_x=-0.055, nudge_y=0.055, size = axis_label_size)+
    geom_text(data = axis_labels,
              aes(x=.data$x3, y=.data$y3, label=.data$rev_label, fontface = 'plain'),  nudge_x=0.055, nudge_y=0.055, size = axis_label_size)+
    theme_void()+
    guides(fill = guide_colorbar(frame.colour = 'black',
                                 ticks.colour = 'black',
                                 title = 'Predicted\nResponse',
                                 show.limits = T))+
    facet_wrap(~ Value, ncol = length(cond_values))+
    coord_fixed()+
    theme(legend.key.size = unit(0.1, 'npc'),
          legend.key.height = unit(0.04, 'npc'),
          legend.title = element_text(face='plain', size = 14, vjust = 0.75),
          plot.subtitle = element_text(face='plain', hjust=0.5, size=14),
          strip.text = element_text(face = 'plain', size =14, vjust = 0.5),
          legend.text = element_text(face = 'plain', size = 12, angle = 45, vjust = 1.2, hjust = 1.2),
          legend.position = 'bottom')

  if(.axis_guides){
    pl <- pl +
      geom_segment(data = axis_labels,
                   aes(x = .data$x1, y = .data$y1, xend = .data$x2, yend = .data$y2), colour='grey',
                   linetype='dashed', linewidth=1, alpha = .75)+
      geom_segment(data = axis_labels,
                   aes(x = .data$x1, y = .data$y1, xend = .data$x3, yend = .data$y3), colour='grey',
                   linetype='dashed', linewidth=1, alpha = .75)+
      geom_segment(data = axis_labels,
                   aes(x = .data$x2, y = .data$y2, xend = rev(.data$x3), yend = rev(.data$y3)), colour='grey',
                   linetype='dashed', linewidth=1, alpha = .75)
  }

  if(.contour_text){
    pl <- pl +
      geom_text_contour(skip=0, breaks = breaks,
                        label.placer = metR::label_placer_fraction(0.15),
                        size=3.5, nudge_x = 0.015, nudge_y = 0.015)

  }

  return(pl)
}
