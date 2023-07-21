#' @keywords internal
#' Utility function for adding incorporating any additional experimental
#' structures into the final data before making predictions
#'
#' @usage NULL
NULL
add_exp_str <- function(data, exp_str){
  # Add experimental structures to communities for predictions
  exp_str_data <- expand.grid(exp_str)

  # Add identifier for each unique combination of experimental structures
  exp_str_data$.add_str_ID <- do.call(paste, c(Map(function(x, y) {paste(x, y, sep=": ")}, names(exp_str_data), exp_str_data), sep="\t"))

  # If any clashes with common columns in data and exp_str
  # Rename columns in data
  common <- intersect(colnames(data), colnames(exp_str_data))
  if(length(common)>0){
    cli::cli_warn(c("Certain names specified in {.var exp_str} are already present in the data.",
                       "i" = "The column{?s} {common} in the data will be given a {.val .data} suffix to avoid conflicts."))
    idx <- which(colnames(data) == common)
    colnames(data)[idx] <- paste0(colnames(data)[idx], ".data")
  }
  # Merge experimental structures with species communities
  data <- left_join(data, exp_str_data, by = character())
  return(data)
}

#' @keywords internal
#' Utility function to return colour-blind friendly colours/palette
#'
#' @usage NULL
NULL
colour_blind_friendly_cols <- function(n){
  # A safe colour-blind palette created by combining Okabe-Ito colours and safe colour pallete from rcartocolor
  safe_colorblind_palette <- c("#009E73", "#D55E00", "#AA4499",  "#0072B2", "#F0E442", "#661100",
                               "#332288", "#CC6677", "#E69F00", "#88CCEE", "#882255", "#DDCC77",
                               "#56B4E9", "#117733", "#CC79A7", "#44AA99", "#484848", "#999933",
                               "#6699CC", "#000000")
  if(n<=20){
    cols <- safe_colorblind_palette[1:n]
  } else {
    message('There are too many colours and they might not all be colour-blind friendly and distiguishable from one another')
    cols <- grDevices::hcl.colors(palette = 'Spectral', n = n)
  }
}

#' @title Returns shades of colours
#'
#' @param colours A character vector of colours recognizable by R, to produces shades of
#' @param shades A numeric vector giving the number of shades for each colour
#'
#' @return A list consisting of hex codes describing the shades of each colour
#' @export
#'
#' @examples
#' ## Shades for a single colour
#' get_shades(c("red"))
#'
#' ## Shades for multiple colours
#' get_shades(c("red", "blue" ,"#A5F8E3", "#808080"), shades = c(2, 3, 4, 5))
#'
#' ## A single value for shade would imply all colours get the same number of shades
#' get_shades(c("red", "blue" ,"#A5F8E3", "#808080"), shades = 2)
get_shades <- function(colours = c("#808080"), shades = 3){
  if(length(shades) != 1 & length(shades) != length(colours)){
    stop(glue::glue("'shades' should be of a vector of length 1 or same as the number of colours ({length(colours)})"))
  }

  if(length(shades) == 1){
    shades <- rep(shades, each = length(colours))
  }

  # Get the L levels to adjust the shades of the colours in HSL format
  light_levels <- lapply(shades, function(x){
    light_level <- switch (x,
                           0,
                           c(-.1, .1),
                           c(-.15, 0, .15),
                           c(-.21, -.07, .07, .21),
                           c(-.2, -.1, 0, .1, .2),
                           c(-0.25, -.15, -0.05, 0.05, 0.15, 0.25),
                           c(-.3, -.2,-.1, 0, .1, .2, .3)
    )

    # If no of shades is greater than 7 then space the colours uniformly along the spectrum
    if(is.null(light_level)){
      light_level <- seq(-.35, .35, length.out = x)
    }
    return(light_level)
  })

  # Convert colours to HSL format
  hsl_cols <- t(plotwidgets::col2hsl(colours))

  # Shades of each colour will be creating by adjusting the L parameter in the HSL respresentation of each colour
  new_cols <- lapply(1:length(colours), function(idx){
    col <- as.vector(hsl_cols[idx, ])
    light_level <- light_levels[[idx]]

    # Adjust for situations when L can become greater than 0 or 1
    threshold <- max(abs(light_level))
    if((col[3] + threshold) > 1){
      col[3] <- 1 - threshold - 0.05    # -0.05 for adding some wiggle room
    }

    if((col[3] - threshold) < 0){
      col[3] <- 0 + threshold + 0.05    # 0.05 for adding some wiggle room
    }

    shades_list <- sapply(light_level, function(x){
      new_col <- col + c(0, 0, x)
      plotwidgets::hsl2col(matrix(new_col, nrow = 3))
    })

    # # Return vector if only one colour is specified
    # if(shades == 1)
    #   shades_list <- shades_list[[1]]

    return(shades_list)
  })

  # Map the shades of colours to the original colour
  names(new_cols) <- colours

  return(new_cols)
}

#' @title Return colour-blind friendly colours
#'
#' @description
#' Utility function to return either a distinct colour-blind friendly colour
#' for each variable or if a functional grouping is specified, then shades of
#' the same colour for species within a functional group
#'
#' @param vars Either a numeric value `n` to get n colours, or a character
#'             vector of values where each value will be mapped to a colour.
#' @param FG A character vector describing the functional grouping to which each
#'           variable belongs. Variables within the same group will have
#'           different shades of the same colour.
#'
#' @return A named vector containing the hex codes of colours
#'
#' @export
#'
#' @examples
#' ## Get n colours
#' get_colours(vars = 4)
#'
#' # Get a color-map for each value specified in vars
#' get_colours(vars = c("p1", "p2", "p3", "p4"))
#'
#' # Group values in vars to have same shades of a colour
#' get_colours(vars = 4, FG = c("G1", "G1", "G2", "G2"))
get_colours <- function(vars, FG = NULL){
  if(missing(vars)){
    cli::cli_abort("{.var vars} cannot be missing.",
                   "i" = "Specify a character vector of variables names to map specific colours to each variable.",
                   "i" = "Or if a specific mapping is not desired, specify a number {.val n} to {.var n} colours.")
  }

  if(is.numeric(vars) && length(vars) == 1){
    vars <- as.character(1:vars)
  }

  nSpecies <- length(vars)
  if(!is.null(FG) & (!is.character(FG) | length(FG) != nSpecies)){
    stop(glue::glue("'FG' should be a character vector having the same length as species vector ({nSpecies}) giving the functional group each species belongs to"))
  }

  if(is.null(FG)){
    cols <- colour_blind_friendly_cols(n = nSpecies)
  } else {
    FG_counts <- table(FG)
    base_cols <- colour_blind_friendly_cols(n = length(FG_counts))

    shades <- get_shades(colours = base_cols, shades = FG_counts)
    names(shades) <- unique(FG)

    cols <- stats::setNames(vector(mode = 'character', length = length(FG)), FG)

    for (sp in names(shades)){
      cols[which(names(cols) == sp)] <- shades[[sp]]
    }
    cols <- unname(cols)
  }
  return(cols)
}

#' @title Get all equi-proportional communities at each level of richness
#'
#' @importFrom dplyr  %>% left_join
#'
#' @param nSpecies Number of species in the design
#' @param species Names for the species. Will be used as column names for the final result. Default is "Sp" followed by column number.
#' @param threshold The maximum number of communities to select for each level of richness for situations when there are too many equi-proportional communities.
#'
#' @return A dataframe consisting all or a random selection of equi-proportional communities at each level of richness
#' @export
#'
#' @examples
#' ## Get all equi-proportional communities for each level of richnness upto 10
#' get_equi_comms(10)
#'
#' ## Change species names
#' get_equi_comms(4, species = c("Lollium perenne", "Chicorum intybus",
#'                               "Trifolium repens", "Trifolium pratense"))
get_equi_comms <- function(nSpecies,
                           species = paste0('Sp', 1:nSpecies),
                           threshold = round(1000000/nSpecies)){

  # Possible levels of richness
  choices <- 1:nSpecies
  # Communities to be used as seed for combinations
  seed_comms <- lapply(choices, function(rich_level){
    c(rep(1, rich_level), rep(0, nSpecies - rich_level))
  })

  chosen_comms <- lapply(seed_comms, function(comm){
    rich_level <- sum(comm)
    nrows <- choose(nSpecies, rich_level)
    # If nCr is less than threshold parameter then get
    # all combinations
    if(nrows <= threshold){
      combs <- matrix(0, nrow = nrows, ncol = nSpecies)
      idx <- utils::combn(1:nSpecies, sum(comm)) %>% t()
      for (i in 1:nrows){
        combs[i, idx[i, ]] <- 1
      }
      combs
      # Other wise choose random combinations of
      # communities
    } else {
      return(t(replicate(threshold, sample(comm))))
    }
  })

  # Bind everything together into a data.frame
  comms <- do.call(rbind, chosen_comms) %>%
    as.data.frame()
  # Add column names and richness variable
  colnames(comms) <- species
  comms$Richness <- rowSums(comms)
  # Convert to proportions
  comms[, 1:nSpecies] <- comms[, 1:nSpecies]/comms$Richness
  return(comms)
}

#' @keywords internal
#' Utility function to return the richness of each community in the data.
#' Particularly useful in dplyr pipeline.
#'
#' @usage NULL
NULL
get_richness <- function(data, species){
  # Special situation if funtion is called in dplyr pipeline as .data doesn't allow to select over multiple columns
  if (inherits(data, "rlang_data_pronoun")) {
    # Add filler column to create data frame
    my_data <- data.frame('filler', data[[species[1]]])
    for (var in species){
      my_data[, var] <- data[[var]]
    }
    data <- my_data
  }
  return(rowSums(data[, species] != 0))
}

#' @keywords internal
#' Utility function to return the evenness of each community in the data.
#' Particularly useful in dplyr pipeline.
#'
#' @usage NULL
NULL
get_evenness <- function(data, species){
  # Special situation if funtion is called in dplyr pipeline as .data doesn't allow to select over multiple columns
  if (inherits(data, "rlang_data_pronoun")) {
    # Add filler column to create data frame
    my_data <- data.frame('filler', data[[species[1]]])
    for (var in species){
      my_data[, var] <- data[[var]]
    }
    data <- my_data
  }
  return(DImodels::DI_data_E_AV(prop = species, data = data)$E)
}

#' @keywords internal
#' Utility function to return the sum of species proportions in each community
#' in the data. Particularly useful in dplyr pipeline.
#'
#' @usage NULL
NULL
get_comm_sum <- function(data, species){
  # Special situation if funtion is called in dplyr pipeline as .data doesn't allow to select over multiple columns
  if (inherits(data, "rlang_data_pronoun")) {
    # Add filler column to create data frame
    my_data <- data.frame('filler', data[[species[1]]])
    for (var in species){
      my_data[, var] <- data[[var]]
    }
    data <- my_data
  }
  return(rowSums(data[, species], na.rm = TRUE))
}

#' @importFrom stats AIC BIC logLik
#' @usage NULL
NULL
AICc <- utils::getFromNamespace("AICc.default", "DImodels") #function(mod) {UseMethod("AICc", mod)} #
#' @usage NULL
NULL
BICc <- utils::getFromNamespace("BICc.default", "DImodels")

# Create vectorized versions of the functions of information criteria
#' @usage NULL
NULL
AIC_vec <- Vectorize(AIC)
#' @usage NULL
NULL
BIC_vec <- Vectorize(BIC)
#' @usage NULL
NULL
#' @usage NULL
NULL
AICc_vec <- Vectorize(AICc)
#' @usage NULL
NULL
BICc_vec <- Vectorize(BICc)
#' @usage NULL
NULL
logLik_vec <- Vectorize(logLik)

#' @usage NULL
NULL
dropInf <- function(x, h) {
  if (any(isInf <- h >= 1)) {
    warning(gettextf("not plotting observations with leverage one:\n  %s",
                     paste(which(isInf), collapse = ", ")), call. = FALSE,
            domain = NA)
    x[isInf] <- NaN
  }
  x
}

#' @keywords internal
#' Function for return a smoothed curve over data points in a plot.
#' Useful for diagnostics plots
#'
#' @usage NULL
NULL
smoothing_fun <- function(x, y) {
  as.data.frame(stats::lowess(x, y, f = 2/3, iter = 3))
}

#' @keywords internal
#' Utility function for rescaling a vector between a specified min and max value
#'
#' @usage NULL
NULL
rescale <- function(x, min = 0, max = 1){
  ((x - min(x))/(max(x) - min(x)) * (max - min)) + min
}

#' @keywords internal
#' Function for checking for the presence of certain columns in the data
#' and print appropriate error message. Useful for plotting function to
#' ensure necessary columns are present in the data.
#'
#' @usage NULL
NULL
check_presence <- function(data, col, message = NULL){
  # Check if any columns are not present
  data_cols <- colnames(data)
  if(is.numeric(col)){
    not_present <- col[!col %in% 1:length(data_cols)]
    if(is.null(message)){
      message <- paste0("Location ", not_present, " is not present in the data.")
    }
  } else {
    not_present <- col[!col %in% data_cols]
    if(is.null(message)){
      message <- paste0("Cannot find `", not_present, "` column in the data.")
    }
  }

  # If column not present then print appropriate error message
  # specified in message
  if(length(not_present)>0){
    cli::cli_abort(message, call = caller_env())
  }
  return(TRUE)
}

#' Default theme for DImodelsVis
#'
#' @importFrom ggplot2 theme_bw theme element_text element_rect %+replace%
#' @importFrom dplyr between
#' @importFrom cli cli_abort
#'
#' @param font_size Base font size for text across the plot
#' @param font_family Font family for text across the plot
#' @param legend One of c("top", "bottom", "left", "right", "none") specifying the position of the legend. The legend position can also be specified as a numeric vector of form c(x, y) with x and y having values between 0 and 1. If specified as a numeric vector the legend within the plotting region where c(0,0) corresponds to the "bottom left" and c(1,1) corresponds to the "top right" position. The default position is "top".
#'
#' @return
#' A ggplot theme object
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#'
#' ggplot(data = iris,
#'        aes(x = Sepal.Length, y = Sepal.Width, colour = Species))+
#'    geom_point(size = 3)+
#'    theme_DI()
#'
theme_DI <- function (font_size = 14, font_family = "",
                      legend = c("top", "bottom", "left", "right", "none")) {
  # Ensure input is appropriate
  sanity_checks(numerics = list("font_size" = font_size),
                characters = list("font_family" = font_family))

  # Legend position
  if(all(is.character(legend))){
    legend <- match.arg(legend)
  } else if(all(is.numeric(legend))){
    if(!all(between(legend, 0, 1))){
      cli_abort(c("If specified as a numeric vector the values for the legend should be between 0 and 1.",
                  "i" = "You specified the values as {.var x = {x}} and {.var {y}}"))
    }
  } else {
    cli_abort(c("{.var legend} should be of type {.cls character} or {.cls numeric}.",
                "i" = "You specified an object of type {.cls {class(legend)}}."))
  }

  # Scale font size of axis text relative to axis titles
  axis_text_size <- font_size * 0.85

  .theme <- theme_bw(base_size = font_size, base_family = font_family) %+replace%
    theme(axis.text = element_text(color = "black", size = axis_text_size),
          strip.background = element_rect(fill = "#465849", colour = "black"),
          legend.position = legend,
          legend.text = element_text(size = axis_text_size),
          complete = TRUE)
  .theme
}


#' @title Add interaction terms used in a Diversity-Interactions (DI) model to new data
#'
#' @importFrom DImodels DI_data
#' @importFrom stats coef
#'
#' @param data A dataframe with species proportions that sum to 1 to create the appropriate interaction structures
#' @param model A Diversity Interactions model object fit by using the \code{\link[DImodels:DI]{DI()}} function from the \code{\link[DImodels:DImodels-package]{DImodels}} package.
#'
#' @description
#' Utility function that accepts a fitted Diversity-Interactions (DI) model
#' object and a data frame and adds the necessary interaction structures to
#' the data for making predictions.
#'
#' @return
#' A data-frame with additional columns describing the interactions terms present in the model object.
#'
#' @export
#'
#' @examples
#' library(DImodels)
#' data(sim1)
#'
#' # Fit different DI models
#' mod1 <- DI(y = "response", prop = 3:6, data = sim1, DImodel = "AV")
#' mod2 <- DI(y = "response", prop = 3:6, data = sim1, DImodel = "FULL")
#' mod3 <- DI(y = "response", prop = 3:6, data = sim1, DImodel = "ADD")
#' mod4 <- DI(y = "response", prop = 3:6, data = sim1,
#'            FG = c("G", "G", "H", "H"), DImodel = "FG")
#'
#' # Create new data for adding interaction terms
#' newdata <- sim1[sim1$block == 1, 3:6]
#' print(head(newdata))
#'
#' add_interaction_terms(data = newdata, model = mod1)
#' add_interaction_terms(data = newdata, model = mod2)
#' add_interaction_terms(data = newdata, model = mod3)
#' add_interaction_terms(data = newdata, model = mod4)
#'
add_interaction_terms <- function(data, model){
  # Ensure model is a DImodels object
  sanity_checks(DImodel = model)

  DImodel_tag <- model$DIcall$DImodel
  if (is.null(DImodel_tag)) {
    DImodel_tag <- "CUSTOM"
  }
  if (DImodel_tag == "CUSTOM") {
    cli::cli_alert_warning("Cannot add interaction terms for a custom DI model.")
  }
  original_data <- model$original_data
  model_data <- eval(model$model)
  prop_cols <- eval(model$DIcall$prop)
  if (is.numeric(prop_cols)) {
    prop <- names(original_data[, prop_cols])
  }
  else {
    prop <- prop_cols
  }

  only_one_row <- nrow(data) == 1
  if (only_one_row) {
    data <- rbind(data, rep(0, ncol(data)))
  }
  theta_flag <- model$coefficients["theta"]
  betas <- coef(model)
  if (!DImodel_tag %in% c("ID", "STR")) {
    if (!is.na(theta_flag)) {
      theta_value <- coef(model)["theta"]
    }
    else {
      theta_value <- 1
    }
    extra_variables <- DI_data(prop = prop, FG = eval(model$DIcall$FG),
                               data = data, theta = theta_value, what = DImodel_tag)
    if (DImodel_tag == "E") {
      updated_data <- data.frame(data, E = extra_variables)
    }
    if (DImodel_tag == "AV") {
      updated_data <- data.frame(data, AV = extra_variables)
    }
    if (DImodel_tag == "ADD") {
      updated_data <- data.frame(data, extra_variables)
    }
    if (DImodel_tag == "FG") {
      data[, "FG_"] <- extra_variables
      updated_data <- data
    }
    if (DImodel_tag == "FULL") {
      updated_data <- data.frame(data, extra_variables,
                                 check.names = FALSE)
    }
  }
  else {
    updated_data <- data
  }
  if (only_one_row) {
    updated_data <- updated_data[1, ]
  }
  return(updated_data)
}

#' @title Add identity effect groups used in a Diversity-Interactions (DI) model to new data
#'
#' @param data A dataframe with species proportions that sum to 1 to create the identity effect groupings.
#' @param model A Diversity Interactions model object fit by using the \code{\link[DImodels:DI]{DI()}} function from the \code{\link[DImodels:DImodels-package]{DImodels}} package.
#'
#' @description
#' Utility function that accepts a fitted Diversity-Interactions (DI) model
#' object along with a data frame and adds the appropriate species identity
#' effect groupings to the data for making predictions.
#'
#' @return
#' A data-frame with additional columns describing the interactions terms present in the model object.
#'
#' @export
#'
#' @examples
#' library(DImodels)
#' data(sim1)
#'
#' # Fit DI models with different ID effect groupings
#' mod1 <- DI(y = "response", prop = 3:6,
#'            data = sim1, DImodel = "AV") # No ID grouping
#' mod2 <- DI(y = "response", prop = 3:6,
#'            data = sim1, DImodel = "AV",
#'            ID = c("ID1", "ID1", "ID2", "ID2"))
#' mod3 <- DI(y = "response", prop = 3:6,
#'            data = sim1, DImodel = "AV",
#'            ID = c("ID1", "ID1", "ID1", "ID1"))
#'
#' # Create new data for adding interaction terms
#' newdata <- sim1[sim1$block == 1, 3:6]
#' print(head(newdata))
#'
#' add_ID_terms(data = newdata, model = mod1)
#' add_ID_terms(data = newdata, model = mod2)
#' add_ID_terms(data = newdata, model = mod3)
#'
add_ID_terms <- function(data, model){
  # Ensure model is a DImodels object
  sanity_checks(DImodel = model)

  # Get original species proportion columns
  prop_cols <- eval(model$DIcall$prop)
  if (is.numeric(prop_cols)) {
    prop <- names(model$original_data[, prop_cols])
  }
  else {
    prop <- prop_cols
  }

  # Check species proportions in the data sum to 1
  sanity_checks(data = data, prop = prop)

  ID <- eval(model$DIcall$ID)
  if(is.null(ID)){
    ID <- paste0(prop, "_ID")
  }

  IDs <-  group_IDs(data = data, prop = prop, ID = ID)

  return(cbind(data, IDs))
}


#' @title Add predictions and confidence interval to data using either statistical model or model coefficients
#'
#' @importFrom insight get_predicted
#'
#' @param data A data-frame containing appropriate values for all the terms in the model
#' @param model A regression model object which will be used to make predictions for the observations in `data`.
#'              Will override `coefficients` if specified.
#' @param coefficients If a regression model is not available (or can't be fit in R), the regression coefficients from a model fit in some other language can be used to calculate predictions.
#'                     However, the user would have to ensure there's an appropriate one-to-one positional mapping between the data columns and the coefficient values.
#'                     Further, it would not be possible to calculate confidence/prediction intervals using this method.
#' @param coeff_cols If `coefficients` are specified and there isn't a one-to-one positional mapping between
#' @param conf.level The confidence level for calculating confidence/prediction intervals. Default is 0.95.
#' @param interval Type of interval to calculate:
#'  \describe{
#'    \item{"none" (default)}{No interval to be calculated.}
#'    \item{"confidence"}{Calculate a confidence interval.}
#'    \item{"prediction"}{Calculate a prediction interval.}
#'  }
#'
#' @return
#' A data-frame with the following additional columns
#'  \describe{
#'    \item{.Pred}{The predicted response for each observation.}
#'    \item{.Lower}{The lower limit of the confidence/prediction interval (if possible) for each observation.}
#'    \item{.Upper}{The lower limit of the confidence/prediction interval (if possible) for each observation.}
#'  }
#'
#' @export
#'
#' @examples
#' library(DImodels)
#' data(sim1)
#'
#' # Fit DI models with different ID effect groupings
#' mod <- lm(response ~ 0 + p1 + p2 + p3 + p4 + p1:p2 + p3:p4, data = sim1)
#'
#' # Create new data for adding predictions
#' newdata <- head(sim1[sim1$block == 1,])
#' print(newdata)
#'
#' # Adding predictions to data with confidence interval
#' add_prediction(data = newdata, model = mod)
#'
#' # Calculate prediction intervals instead
#' add_prediction(data = newdata, model = mod, interval = "prediction")
#'
#' # Change confidence level
#' add_prediction(data = newdata, model = mod, interval = "prediction",
#'                conf.level = 0.99)
#'
#' # Use model coefficients for prediction
#' coeffs <- mod$coefficients
#'
#' # Would now have to add columns corresponding to each coefficient in the
#' # data and ensure there is a desired mapping between data columns and
#' # the coefficients. Also confidence intervals cannot be calculated.
#' newdata$`p1:p2` = newdata$p1 * newdata$p2
#' newdata$`p3:p4` = newdata$p3 * newdata$p4
#'
#' # If the coefficients are named then the function will try to
#' # perform matching between data columns and the coefficients
#' add_prediction(data = newdata, coefficients = coeffs)
#'
#' # However, if the coefficients are not named
#' coeffs <- unname(coeffs)
#'
#' # The user would have to manually specify the subset
#' # of data columnsvarranged according to the coefficients
#' subset_data <- newdata[, c(3:6, 8,9)]
#' add_prediction(data = subset_data, coefficients = coeffs)
#'
#' # Or specify a selection (either by name or index) in coeff_cols
#' add_prediction(data = newdata, coefficients = coeffs,
#'                coeff_cols = c("p1", "p2", "p3", "p4", "p1:p2", "p3:p4"))
#'
add_prediction <- function(data, model = NULL,
                           coefficients = NULL, coeff_cols = NULL,
                           interval = c("none","confidence", "prediction"),
                           conf.level = 0.95){

  # Ensure model or coefficients are specified correctly
  check_data_functions(model = model, coefficients = coefficients)

  # Choose appropriate interval for prediction CI
  interval <- match.arg(interval)

  # Branch here if model is specified
  if(!is.null(model)){
    # Prediction with SE
    if(inherits(model, "DI")){
      # TODO --------------------------------------------------------------------
      # Remember to add functionality for getting both, confidence and prediction interval
      if(interval != "none"){
        preds <- suppressWarnings(predict(model, newdata = data, se.fit = TRUE))
        data <- data %>%
          mutate(".Pred" := preds$fit,
                 ".Lower" = preds$fit - round(qnorm((1 - conf.level)/2, lower.tail = F), 2) * preds$se.fit,
                 ".Upper" = preds$fit - round(qnorm((1 + conf.level)/2, lower.tail = F), 2) * preds$se.fit)
      } else{
        preds <- suppressWarnings(predict(model, newdata = data, se.fit = FALSE))
        data <- data %>%
          mutate(".Pred" := preds)
      }
    } else {
      preds <- suppressMessages(insight::get_predicted(x = model, data = data))

      data <- data %>%
        mutate(".Pred" := as.numeric(preds))

      if(interval != "none"){
        CIs <- suppressMessages(insight::get_predicted_ci(x = model,
                                                          predictions = preds,
                                                          data = data,
                                                          ci = conf.level,
                                                          ci_type = interval))

        if(ncol(CIs) > 1){
          data <- data %>% mutate(".Lower" = CIs[, "CI_low"],
                                  ".Upper" = CIs[, "CI_high"])
        }
      }
    }
  }
  # Branch here if regression coefficients are specified
  else if(!is.null(coefficients)){
    # If coefficients are named then order data columns accordingly
    if(!is.null(names(coefficients))){
      if(!all(names(coefficients) %in% colnames(data))){
        cli::cli_abort(c("All coefficient names should be present in the data as columns.",
                         "i" = "{.var {names(coefficients)[!names(coefficients) %in% colnames(data)]}} {?was/were} not present in the data."))
      }
      # Create X matrix
      X_matrix <- as.matrix(sapply(data[, names(coefficients)], as.numeric))
    }
    # If a selection is specified from the data
    else if (!is.null(coeff_cols)){
      if(length(coefficients)!=length(coeff_cols)){
        cli::cli_abort(c("The number of columns specified for selection in {.var coeff_cols} should be same as the number of coefficients.",
                         "i" = "The were {length(coefficients)} coefficients while {.var coeff_cols} specified {length(coeff_cols)} column{?s} to select."))
      }
      # select the specified columns
      X_cols <- data %>% select(all_of(coeff_cols))
      # Created X_matrix
      X_matrix <- as.matrix(sapply(X_cols, as.numeric))
    }
    # If neither coefficients are named nor a selection provided then
    # assume the user has specified everything correctly and calculate predictions
    else {
      if(ncol(data)!=length(coefficients)){
        cli::cli_abort(c("The number of columns in data should be the same as the number of coefficients.",
                         "i" = "The were {length(coefficients)} coefficients while data had {ncol(data)} columns.",
                         "i" = "Consider giving names to the coefficient vector specified in {.var coefficients} corresponding to the respective data columns or providing a selection of columns in {.var coeff_cols} corresponding (in sequential order) to each coefficient."))
      }
      # Create X_matrix
      X_matrix <- as.matrix(sapply(data, as.numeric))
    }
    # Calculate predictions
    preds <- X_matrix %*% coefficients
    data <- data %>%
      mutate(".Pred" := preds,
             ".Upper" = preds,
             ".Lower" = preds)
  }
  return(data)
}

# Function to print generic message if model is not DI
model_not_DI <- function(call_fn){
  data_fn <- paste0(call_fn, "_data")
  plot_fn <- paste0(call_fn, "_plot")
  cli::cli_abort(c("{.var model} should be a regression model fit using the {.help [{.fun DI}](DImodels::DI)} function from
                     the {.help [{.pkg DImodels}](DImodels::DImodels)} package or an object extending the {.cls DI} class.",
                   "i" = "If your model cannot be fit using the {.help [{.fun DI}](DImodels::DI)} function,
                     manually call the {.help [{.fn {data_fn}}]} function followed by the
                     {.help [{.fn {plot_fn}}]} function to create the plot."))
}

# Import function from DImodels
group_IDs <- getFromNamespace("group_IDs", "DImodels")
