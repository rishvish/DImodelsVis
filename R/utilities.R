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

    sapply(light_level, function(x){
      new_col <- col + c(0, 0, x)
      plotwidgets::hsl2col(matrix(new_col, nrow = 3))
    })
  })

  # Map the shades of colours to the original colour
  names(new_cols) <- colours

  return(new_cols)
}

#' @importFrom grDevices col2rgb
#' @usage NULL
NULL
areColours <- function(colours) {
  sapply(colours, function(colour) {
    tryCatch(is.matrix(col2rgb(colour)),
             error = function(e) FALSE)
  })
}

#' @usage NULL
NULL
get_colours <- function(species, FG = NULL){
  if(missing(species) | !is.character(species)){
    stop("'species' should be a character vector giving the names of the species in the model")
  }

  nSpecies <- length(species)
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

#' @importFrom stats AIC BIC logLik
#' @usage NULL
NULL
AICc <- utils::getFromNamespace("AICc.default", "DImodels")
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

#' @usage NULL
NULL
smoothing_fun <- function(x, y) {
  as.data.frame(stats::lowess(x, y, f = 2/3, iter = 3))
}

#' @usage NULL
NULL
rescale <- function(x, min = 0, max = 1){
  ((x - min(x))/(max(x) - min(x)) * (max - min)) + min
}
