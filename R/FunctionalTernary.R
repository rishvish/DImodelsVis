#' @title Conditional ternary diagrams at functional group level
#'
#' @description
#' `Need to think of a simple description. Will be updated soon`
#'
#'
#' @inheritParams conditional_ternary_data
#' @param FG A character vector specifying the functional grouping of the species in the design.
#' @param values A vector specifying the proportional split between the species within a functional group.
#'               The default is to split the functional group proportional equally between each species.
#' @param cond_values A vector of numbers between 0 and 1 describing the values at which to takes slices of the ternary.
#'                    Necessary only when there are more than 3 functional groups.
#' @param cond_FG A character vector specifying the functional groups to slice the ternary on if there
#'                are more than three functional groups. If there are more than three functional groups,
#'                by default, the ternary diagrams will be conditioned on `add more stuff here`
#' @inheritDotParams add_prediction -data
#'
#' @return A data-frame containing compositional columns with names specified
#'         in `FG` and `prop` parameters along with any additional columns
#'         specified in `exp_str` parameter and the following columns appended
#'         at the end.
#'  \describe{
#'    \item{.x}{The x-projection of the points within the ternary.}
#'    \item{.y}{The y-projection of the points within the ternary.}
#'    \item{.add_str_ID}{An identifier column for grouping the cartesian product
#'                       of all additional columns specified in `exp_str`
#'                       parameter (if `exp_str` is specified).}
#'    \item{.Sp}{An identifier column specifying the functional group along
#'               which the high dimensional simplex is sliced (if there are
#'               more than 3 functional groups).}
#'    \item{.Value}{The value (between 0 and 1) along the direction of functional
#'                  group in `.Sp` at which the high dimensional simplex is sliced.}
#'    \item{.Facet}{An identifier column formed by combining `.Sp` and `.value`
#'                  to group observations within a specific slice of the
#'                  high dimensional simplex.}
#'    \item{.Pred}{The predicted response for each community.}
#'    \item{.Lower}{The lower limit of the prediction/confidence interval
#'                  for each observation.}
#'    \item{.Upper}{The upper limit of the prediction/confidence interval
#'                  for each observation.}
#'  }
#'
#' @export
#'
#' @examples
#' library(DImodels)
#'
#' ## Load data
#' data(sim4)
#'
#' ## Fit model
#' mod <- glm(response ~ 0 + (p1 + p2 + p3 + p4 + p5 + p6)^2, data = sim4)
#'
#' ## Create data
#' ## By default an equal split of the FG will be assumed for each species
#' head(functional_ternary_data(prop = paste0("p", 1:6),
#'                              FG = c("G1", "G1", "G2", "G2", "G3", "G3"),
#'                              resolution = 1,
#'                              model = mod))
#'
#' ## Change the proportional split of species within an FG by using `values`
#' ## `values` takes a numeric vector where the position of each element
#' ## describes the proportion of the corresponding species within the
#' ## corresponding FG
#' ## For examples this vector describes, 0% and 100% of p1 and p2, respectively
#' ## in G1, 100% and 0% of p3 and p4, respectively in G2 and 70% and 30% of
#' ## p5 and p6, respectively in G3.
#' vals <- c(0, 1, 1, 0, 0.7, 0.3)
#' head(functional_ternary_data(prop = paste0("p", 1:6),
#'                              FG = c("G1", "G1", "G2", "G2", "G3", "G3"),
#'                              values = vals,
#'                              resolution = 1,
#'                              model = mod))
#'
#' ## Can also add any additional experimental structures
#' head(functional_ternary_data(prop = paste0("p", 1:6),
#'                              FG = c("G1", "G1", "G2", "G2", "G3", "G3"),
#'                              exp_str = list("treatment" = c("50", "150")),
#'                              values = vals,
#'                              model = mod,
#'                              resolution = 1))
#'
#' ## It could be desirable to take the output of this function and add
#' ## additional variables to the data before making predictions
#' ## Use `prediction = FALSE` to get data without any predictions
#' head(functional_ternary_data(prop = paste0("p", 1:6),
#'                              FG = c("G1", "G1", "G2", "G2", "G3", "G3"),
#'                              values = vals,
#'                              resolution = 1,
#'                              prediction = FALSE))
#'
#' ## When there are more than three functional groups, we can
#' ## fix the value of one of the functional groups using `cond_FG`
#' ## and `cond_values` like conditional_ternary
#' head(functional_ternary_data(prop = paste0("p", 1:6),
#'                              FG = c("G1", "G1", "G2", "G2", "G3", "G4"),
#'                              resolution = 1,
#'                              cond_FG = "G4",
#'                              cond_values = c(0.2, 0.5),
#'                              model = mod,
#'                              prediction = TRUE))
functional_ternary_data <- function(prop, FG,
                                    values = NULL,
                                    cond_values = NULL,
                                    cond_FG = NULL,
                                    exp_str = list(),
                                    resolution = 3,
                                    prediction = TRUE, ...){

  #Sanity Checks
  if(missing(prop)){
    cli::cli_abort(c("{.var prop} cannot be empty.",
                     "i" = "Specify a character vector indicating the model coefficients corresponding to variable proportions."))
  }

  # Ensure FG is specified
  if(missing(FG)){
    # if(!is.null(eval(model$DIcall$FG))){
    #   FG <- eval(model$DIcall$FG)
    # } else {
      cli::cli_abort(c("The {.var FG} argument cannot be empty.",
                       "i" = "The {.var FG} argument should be specified as a character vector of same length as the {.var prop} argument, specifying the functional group to which each species in prop belongs."))
    # }
  }

  # Check inputs of function are appropriate and return default values
  def_vals <- FG_sanity_checks(prop = prop, FG = FG,
                   values = values,
                   cond_values = cond_values,
                   cond_FG = cond_FG,
                   exp_str = exp_str,
                   resolution = resolution,
                   prediction = prediction, ...)

  values <- def_vals$values
  all_FGs <- def_vals$all_FGs
  tern_FGs <- def_vals$tern_FGs
  cond_FG <- def_vals$cond_FG
  cond_values <- def_vals$cond_values

  # Create base ternary data to be plotted
  triangle <- ternary_data(prop = tern_FGs[1:3],
                           exp_str = exp_str,
                           resolution = resolution,
                           prediction = FALSE)

  # If there are only three functional groups in the model
  # our data is ready and we can make predictions from the model after accounting for species proportions
  if (length(all_FGs) == 3){
    triangle[, prop] <- 0
    focus <- tern_FGs
    FG_mapping <- get_FG_value_mapping(FG, focus = focus)

    for (FGs in focus){
      ids <- FG_mapping[[FGs]]
      triangle[, prop[ids]] <- (triangle[, FGs] %o% values[ids])
    }

    cond_data <- triangle

    if(prediction){
      cond_data <- add_prediction(data = cond_data, ...)
    }
    # cond_data <- add_prediction(data = triangle,
    #                             model = model,
    #                             coefficients = coefficients,
    #                             coeff_cols = coeff_cols,
    #                             pred_name = pred_name,
    #                             conf.level = conf.level)
  } else {
    # Rescaling species to be shown in ternary for each species
    # in the cond_FG parameter according to the values specified
    # in the cond_values parameter

    cond_data <- lapply(cli_progress_along(cond_values, name = "Preparing data"), function(idx){
      x <- cond_values[idx]
      if(x == 0){
        scaled_data <- triangle
      } else {
        # Scale proportion of species within the ternary
        scaled_data <- triangle %>%
          mutate(!! tern_FGs[1] := rescale(!!sym(tern_FGs[1]),
                                           min = 0, max = 1-x),
                 !! tern_FGs[2] := rescale(!!sym(tern_FGs[2]),
                                           min = 0, max = 1-x),
                 !! tern_FGs[3] := rescale(!!sym(tern_FGs[3]),
                                           min = 0, max = 1-x))
      }

      # Add remaining species in the data
      remaining_FG <- matrix(0, ncol=length(cond_FG), nrow= nrow(triangle))
      colnames(remaining_FG) <- cond_FG
      scaled_data <- cbind(scaled_data, remaining_FG)

      # Update value of species which we are conditioning on
      species_data <- lapply(cond_FG, function(cond_sp){
        # Add identifier for grouping data for a particular species
        sp_data <- scaled_data %>%
          mutate(.Sp = cond_sp,
                 .Value =  x,
                 .Facet = paste0(cond_sp, ' = ', x))

        # To avoid any rounding issues and ensure all species proportions sum to 1
        sp_data[, cond_sp] <- 1 - rowSums(sp_data[, tern_FGs[1:3]])

        sp_data[, prop] <- 0
        focus <- c(tern_FGs[1:3], cond_sp)
        FG_mapping <- get_FG_value_mapping(FG, focus = focus)

        for (FGs in focus){
          ids <- FG_mapping[[FGs]]
          # Converting to df because when there is only 1 prop, %o% returns
          # a vector which causes conflicts
          sp_data[, prop[ids]] <- (sp_data[, FGs] %o% values[ids]) %>%
                                    as.data.frame() %>%
                                    `colnames<-`(prop[ids])
        }

        # Predicting the response for the communities
        if(prediction){
          sp_data <- add_prediction(data = sp_data, ...)
        }
        # sp_data <- add_prediction(data = sp_data,
        #                           model = model,
        #                           coefficients = coefficients,
        #                           coeff_cols = coeff_cols,
        #                           pred_name = pred_name,
        #                           conf.level = conf.level,
        #                           interval = "none")

        # Need to return data via subset of rows as bind_rows fails otherwise
        sp_data[1:nrow(sp_data),]
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  cli::cli_alert_success("Finished data preparation.")
  return(cond_data %>% select(all_of(c(FG, prop, names(exp_str))), everything()))
  #return(list(data = cond_data, cond_FG = cond_FG, cond_values = cond_values))
}

#' @title Conditional ternary diagrams at functional group level
#'
#' @description
#' `Will be updated soon`
#'
#' @inheritParams conditional_ternary_plot
#'
#' @inherit ternary_plot return
#'
#' @export
#'
#' @examples
#' library(DImodels)
#'
#' ## Load data
#' data(sim4)
#'
#' ## Fit model
#' mod <- glm(response ~ 0 + (p1 + p2 + p3 + p4 + p5 + p6)^2, data = sim4)
#'
#' ## Create data
#' ## By default an equal split of the FG will be assumed for each species
#' plot_data <- functional_ternary_data(prop = paste0("p", 1:6),
#'                              FG = c("G1", "G1", "G2", "G2", "G3", "G3"),
#'                              resolution = 1,
#'                              model = mod)
#' functional_ternary_plot(plot_data,
#'                         tern_labels = c("G1", "G2", "G3"))
#'
#' ## Change the proportional split of species within an FG by using `values`
#' ## `values` takes a numeric vector where the position of each element
#' ## describes the proportion of the corresponding species within the
#' ## corresponding FG. This is equvialent to looking at a different slice
#' ## of the high dimensional simplex.
#' ## For examples this vector describes, 0% and 100% of p1 and p2, respectively
#' ## in G1, 100% and 0% of p3 and p4, respectively in G2 and 70% and 30% of
#' ## p5 and p6, respectively in G3.
#' vals <- c(0, 1, 1, 0, 0.7, 0.3)
#' plot_data <- functional_ternary_data(prop = paste0("p", 1:6),
#'                              FG = c("G1", "G1", "G2", "G2", "G3", "G3"),
#'                              values = vals,
#'                              resolution = 1,
#'                              model = mod)
#'
#' functional_ternary_plot(plot_data,
#'                         tern_labels = c("G1", "G2", "G3"))
#'
#' ## When there are more than three functional groups, we can
#' ## fix the value of one of the functional groups using `cond_FG`
#' ## and `cond_values` like conditional_ternary
#' plot_data <- functional_ternary_data(prop = paste0("p", 1:6),
#'                              FG = c("G1", "G1", "G2", "G2", "G3", "G4"),
#'                              resolution = 1,
#'                              cond_FG = "G4",
#'                              cond_values = c(0.2, 0.5),
#'                              model = mod,
#'                              prediction = TRUE)
#'
#' functional_ternary_plot(plot_data,
#'                         tern_labels = c("G1", "G2", "G3"))
functional_ternary_plot <- function(data,
                                    nlevels = 7,
                                    colours = NULL,
                                    lower_lim = NULL,
                                    upper_lim = NULL,
                                    tern_labels = c("FG1", "FG2", "FG3"),
                                    contour_text = TRUE,
                                    show_axis_labels = TRUE,
                                    show_axis_guides = FALSE,
                                    axis_label_size = 4,
                                    vertex_label_size = 5){
  if(missing(data)){
    cli::cli_abort(c("{.var data} cannot be missing.",
                    "i" = "Specify the output of {.fn conditional_ternary_data} function."))
  }

  pl <- conditional_ternary_plot(data = data,
                                 nlevels = nlevels,
                                 colours = colours,
                                 tern_labels = tern_labels,
                                 lower_lim = lower_lim,
                                 upper_lim = upper_lim,
                                 contour_text = contour_text,
                                 show_axis_labels = show_axis_labels,
                                 show_axis_guides = show_axis_guides,
                                 axis_label_size = axis_label_size,
                                 vertex_label_size = vertex_label_size)
  return(pl)
}

#' @title Conditional ternary diagrams at functional group level
#'
#' @importFrom stats setNames
#'
#' @param model A Diversity Interactions model object fit by using the \code{\link[DImodels:DI]{DI()}} function from the \code{\link[DImodels:DImodels-package]{DImodels}} package.
#' @inheritParams functional_ternary_data
#' @inheritParams ternary_plot
#' @inheritParams model_diagnostics
#'
#' @inherit response_contributions return
#'
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
#' functional_ternary(m1, FG = c("G","G","L","H","H","M"),
#'                    resolution = 1, cond_FG = "M",
#'                    values = c(0, 1, 1, 0.5, 0.5, 1))
functional_ternary <- function(model,
                               FG,
                               values = NULL,
                               cond_values = NULL,
                               cond_FG = NULL,
                               exp_str = list(),
                               resolution = 3,
                               nlevels = 7,
                               colours = NULL,
                               lower_lim = NULL,
                               upper_lim = NULL,
                               contour_text = TRUE,
                               show_axis_labels = TRUE,
                               show_axis_guides = FALSE,
                               axis_label_size = 4,
                               vertex_label_size = 5,
                               nrow = 0, ncol = 0){

  # Ensure specified model is fit using the DI function
  if(missing(model) || !inherits(model, "DI")){
    model_not_DI(call_fn = "functional_ternary")
  }

  # Get data used to fit the model
  og_data <- model$original_data

  # Get all species in the model
  species <- eval(model$DIcall$prop)
  if(is.numeric(species)){
    species <- colnames(og_data)[species]
  }

  # Ensure FG is specified
  if(missing(FG)){
    if(!is.null(eval(model$DIcall$FG))){
      FG <- eval(model$DIcall$FG)
    } else {
      cli::cli_abort(c("The {.var FG} argument cannot be empty.",
                       "i" = "The {.var FG} argument should be specified as a character vector of same length as the {.var prop} argument, specifying the functional group to which each species in prop belongs."))
    }
  }

  # Create data in appropriate format for plotting
  plot_data <- functional_ternary_data(prop = species, FG = FG,
                                       model = model,
                                       values = values,
                                       cond_values = cond_values,
                                       cond_FG = cond_FG,
                                       exp_str = exp_str,
                                       resolution = resolution)

  all_FGs <- unique(FG)
  if(length(all_FGs > 3)){
    if(is.null(cond_FG)){
      cond_FG <- all_FGs[4:length(all_FGs)]
    }
    if(is.null(cond_values)){
      cond_values <- c(0.2, 0.5, 0.8)
    }
  }

  # Labels for the ternary
  tern_labels <- colnames(plot_data)[1:3]

  if(length(exp_str) > 0){
    ids <- unique(plot_data$.add_str_ID)
    plots <- lapply(cli_progress_along(1:length(ids), name = "Creating plot",
                                       format = paste0(
                                         "{pb_spin} Creating plot ",
                                         "[{pb_current}/{pb_total}]   ETA:{pb_eta}"
                                       )),
                    function(i){
                      data <- plot_data %>% filter(.data$.add_str_ID == ids[i])
                      plot <- functional_ternary_plot(data = data,
                                                      nlevels = nlevels,
                                                      colours = colours,
                                                      tern_labels = tern_labels,
                                                      lower_lim = lower_lim,
                                                      upper_lim = upper_lim,
                                                      contour_text = contour_text,
                                                      show_axis_labels = show_axis_labels,
                                                      show_axis_guides = show_axis_guides,
                                                      axis_label_size = axis_label_size,
                                                      vertex_label_size = vertex_label_size) +
                        labs(subtitle = ids[i])
                    })
    if(length(plots) > 1){
      plot <- new("ggmultiplot", plots = plots, nrow = nrow, ncol = ncol)
    } else {
      plot <- plots[[1]]
    }
    cli::cli_alert_success("Created all plots.")
  } else {
    plot <- functional_ternary_plot(data = plot_data,
                                    nlevels = nlevels,
                                    colours = colours,
                                    tern_labels = tern_labels,
                                    lower_lim = lower_lim,
                                    upper_lim = upper_lim,
                                    contour_text = contour_text,
                                    show_axis_labels = show_axis_labels,
                                    show_axis_guides = show_axis_guides,
                                    axis_label_size = axis_label_size,
                                    vertex_label_size = vertex_label_size)
    cli::cli_alert_success("Created plot.")
  }
  return(plot)
}

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

FG_sanity_checks <- function(prop, FG,
                             values = NULL,
                             cond_values = NULL,
                             cond_FG = NULL,
                             exp_str = list(),
                             resolution = 3,
                             prediction = TRUE, ...){

  if(!is.character(FG)){
    cli::cli_abort(c("The {.var FG} argument should be specified as a character vector of same length as the {.var prop} argument, specifying the functional group to which each species in prop belongs.",
                     "i" = "{.var FG} was specified as a {.cls {class(FG)}} object."))
  }

  if(length(prop) != length(FG)){
    cli::cli_abort(c("The {.var FG} argument should be specified as a character vector of same length as the {.var prop} argument, specifying the functional group to which each species in prop belongs.",
                     "i" = "{.var FG} has length {length(FG)} while {.var prop} has length {length(prop)}."))
  }

  if(is.null(values)){
    values <- get_FG_values(FG)
    cli::cli_warn("The proportional split of species in the functional groups was not specified in {.var values},
                  assuming an equal split for species in each functional group.")
  } else {
    if(!is.numeric(values)){
      cli::cli_abort(c("{.var values} should be a numeric vector with values between 0 and 1 specifying the proportion of each species within a function group.",
                       "i" = "{.var values} was specified as a {.cls {class(values)}} object."))
    }

    if(!all(between(values, 0, 1))){
      cli::cli_abort(c("{.var values} should be a numeric vector with values between 0 and 1 specifying the proportion of each species within a function group.",
                       "i" = "{.var values} was specified with value{?s} {.val {as.character(values)}}."))
    }

    if(length(values) != length(FG)){
      cli::cli_abort(c("{.var values} should have the same length as the {.var FG} argument.",
                       "i" = "{.var values} has length {length(values)} while {.var FG} has length {length(FG)}."))
    }

    sp_props_in_FG <- check_FG_values(FG, values)

    if(!all(sp_props_in_FG$bool)){
      faults <- sp_props_in_FG$sums[sp_props_in_FG$bool == F]
      cli::cli_abort(c("The species proportions within a functional group should sum to 1.",
                       "i" = "The proportions for {names(faults)} equal {faults}. respectively."))
    }
  }

  # Get species to be shown in the ternary diagram
  all_FGs <- unique(FG)

  if(length(all_FGs) < 3){
    cli::cli_abort(c("Ternary diagrams cannot be created for less than 3 unique functional groups",
                     "i" = "Currently only {length(all_FGs)} unique functional groups are specified in {.var FG}."))
  }

  # Check for conditional FG and conditioning values
  if(length(all_FGs) > 3){
    if(!is.null(cond_FG)){

      if(!is.character(cond_FG)){
        cli::cli_abort(c("{.var cond_FG} should be a character vector specifying the functional group to condition the ternary diagram when there are more than 3 functional groups in the model.",
                         "i" = "{.var cond_FG} was specified as a {.cls {cond_FG}} object."))
      }

      if(!all(cond_FG %in% all_FGs)){
        cli::cli_abort(c("The value specified in {.var cond_FG} was not found in {.var FG}.",
                         "i" = "{.var cond_FG} has value {.val {cond_FG}} while {.var FG} has values {.val {all_FGs}}."))
      }
    } else {
      cli::cli_warn(c("More than three functional groups in model and no function group specified to condition on.",
                      "i" = "Choosing {.val {all_FGs[1:3]}} to show in ternary and conditioning on {.val {all_FGs[4:length(all_FGs)]}}"))
      cond_FG <- all_FGs[4:length(all_FGs)]
    }

    if(!is.null(cond_values)){
      if(!all(is.numeric(cond_values))){
        cli::cli_abort(c("{.var cond_values} should be a numeric vector with values between 0 and 1 describing the values at which to condition the n-dimensional simplex space.",
                         "i" = "{.var cond_values} was specified as a {.cls {cond_values}} object."))
      }

      if(!all(between(cond_values, 0, 1))){
        cli::cli_abort(c("{.var cond_values} should be a numeric vector with values between 0 and 1 describing the values at which to condition the n-dimensional simplex space.",
                         "i" = "The value{?s} specified in {.var cond_values} {?was/were} {.val {as.character(cond_values)}}."))
      }
    } else {
      cli::cli_warn(c("Choosing default values of 0.2, 0.5 and 0.8 for conditioning the n-dimensional simplex space."))
      cond_values <- c(0.2, 0.5, 0.8)
    }
  }

  tern_FGs <- setdiff(unique(FG), cond_FG)
  if(length(tern_FGs) > 3){
    cli::cli_warn(c("After accounting for the values specified in {.var cond_FG}, there are more than three functional groups to show on the ternary diagram.",
                    "i" = "Creating the diagram for functional groups `{tern_FGs[1]}`, `{tern_FGs[2]}`, and `{tern_FGs[3]}` and assuming all of {.var {tern_FGs[4:length(tern_FGs)]}} to be 0, but it might not be very informative."))
  }
  return(list("values" = values, "all_FGs" = all_FGs, "tern_FGs" = tern_FGs,
              "cond_FG" = cond_FG, "cond_values" = cond_values))
}
