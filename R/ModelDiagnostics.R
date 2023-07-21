#' @title Regression diagnostics plots with pie-glyphs
#' @description
#' This function returns regression diagnostics plots for a model with points
#' replaced by pie-chart glyphs making it easier to track various data points
#' in the plots.
#'
#' @importFrom PieGlyph geom_pie_glyph pieGrob
#' @importFrom dplyr bind_rows desc filter arrange slice
#'                   select all_of mutate .data
#' @importFrom stats family residuals qqnorm weights quantile qnorm
#' @importFrom ggplot2 fortify xlim ylim xlab ylab ggtitle geom_hline
#'                     geom_text geom_abline geom_linerange geom_segment
#'                     aes ggplot scale_fill_manual unit labs guides
#'                     scale_y_continuous expansion
#' @importFrom methods new
#' @importFrom ggtext geom_richtext
#' @importFrom cli cli_progress_bar cli_progress_update cli_process_done
#' @importClassesFrom ggfortify ggmultiplot
#'
#' @param model A statistical regression model object.
#' @param which A subset of the numbers 1 to 6, by default 1, 2, 3, and 5,
#'              referring to
#'                  1. - "Residuals vs Fitted", aka ‘Tukey-Anscombe’ plot
#'                  2. - "Normal Q-Q" plot, an enhanced qqnorm(resid(.))
#'                  3. - "Scale-Location"
#'                  4. - "Cook's distance"
#'                  5. - "Residuals vs Leverage"
#'                  6. - "Cook's dist vs Lev./(1-Lev.)"
#' @param npoints Number of points to be labelled in each plot, starting
#'                with the most extreme.
#' @param cook.levels A numeric vector specifying levels of Cook's distance
#'                    at which to draw contours.
#' @param nrow Number of rows in which to arrange the final plot
#' @param ncol Number of columns in which to arrange the final plot
#' @param prop A character vector giving names of columns containing
#'             proportions to show in the pie-chart glyphs. If not specified,
#'             black points (geom_point) will be shown for each observation in
#'             the model. Note: `\code{prop}` can be left blank and will be
#'             interpreted if model is a \code{Diversity-Interactions (DI)}
#'             model object fit using the\code{\link[DImodels:DI]{DI()}}
#'             function from the \code{\link[DImodels:DImodels-package]{DImodels}}
#'             package.
#' @param pie_radius A numeric value specifying the radius (in cm) for the
#'                   pie-chart glyphs.
#' @param pie_colours A character vector specifying the colours for the slices
#'                    within the pie-chart glyphs.
#' @param label_size A numeric value specifying the size of the labels
#'                   identifying extreme observations.
#' @param points_size A numeric value specifying the size of points (when
#'                    pie-glyphs not shown) shown in the plots.
#' @param only_extremes A logical value indicating whether to show pie-glyphs
#'                      only for extreme observations.
#'
#' @return A ggmultiplot (ggplot if single plot is returned) class object
#' @export
#'
#' @examples
#' library(DImodels)
#'
#' ## Load data
#' data(sim1)
#'
#' ## Fit model
#' mod1 <- lm(response ~ 0 + (p1 + p2 + p3 + p4)^2, data = sim1)
#'
#' \dontrun{
#' ## Get diagnostics plot
#' ## Recommend to store plot in a variable, to access individual plots later
#' diagnostics <- model_diagnostics(mod1, prop = c("p1", "p2", "p3", "p4"))
#' print(diagnostics)
#'
#' ## Access individual plots
#' print(diagnostics[[1]])
#' print(diagnostics[[4]])
#'
#' ## Change plot arrangement
#' model_diagnostics(mod1, prop = c("p1", "p2", "p3", "p4"),
#'                   which = c(1, 3), nrow = 2, ncol = 1)
#'
#' ## Show only extreme points as pie-glyphs to avoid overplotting
#' model_diagnostics(mod1, prop = c("p1", "p2", "p3", "p4"),
#'                   which = 2, npoints = 5, only_extremes = TRUE)
#'
#' ## If model is a DImodels object, the don't need to specify prop
#' DI_mod <- DI(y = "response", prop = c("p1", "p2", "p3", "p4"),
#'              DImodel = "FULL", data = sim1)
#' model_diagnostics(DI_mod, which = 1)
#' }
model_diagnostics <- function(model, which = c(1,2,3,5), prop = NULL,
                              npoints = 3, cook.levels = c(0.5, 1),
                              pie_radius = 0.25, pie_colours = NULL,
                              only_extremes = FALSE, label_size = 4,
                              points_size = 3, nrow = 0, ncol = 0){
  # Sanity checks
  sanity_checks(model = model,
                numerics = list("which" = which, "npoints" = npoints,
                                "cook.levels" = cook.levels,
                                "pie_radius" = pie_radius,
                                "label_size" = label_size,
                                "points_size" = points_size,
                                "nrow" = nrow, "ncol" = ncol),
                booleans = list("only_extremes" = only_extremes),
                colours = pie_colours,
                unit_lengths = list("pie_radius" = pie_radius,
                                    "label_size" = label_size,
                                    "points_size" = points_size,
                                    "nrow" = nrow, "ncol" = ncol,
                                    "npoints" = npoints,
                                    "only_extremes" = only_extremes))

  # Data will all diagnostic metrics added
  plot_data <- fortify(model = model)

  if(is.null(prop)){
    if(inherits(model, "DI")){
      # Get original data used to fit the model
      original_data <- model$original_data

      # Get all species in the model
      #IDs <- eval(model$DIcall$ID)
      prop <- eval(model$DIcall$prop)
      plot_data <- cbind(plot_data, original_data[, prop])

      model_species <- original_data[, prop] %>% colnames()
      pies <- TRUE
    }  else {
      pies <- FALSE
    }
  } else {
    # Ensure proportions are present in the model and they sum to 1
    sanity_checks(data = plot_data, prop = prop)
    model_species <- plot_data %>% select(all_of(prop)) %>% colnames()
    pies <- TRUE
  }

  # If only_extremes in TRUE, the show pies only for extreme
  if(only_extremes){
    if(!pies){
      cli_warn(c("No values were specified in {.var pies} nor can they
                 be interepreted from the model object.",
                 "i" = "The extreme values will be shown a black points
                        and be identified with their observation numbers."))
      only_extremes <- FALSE
    }
  }

  # Sanity checks from plot.lm function
  if (!is.numeric(which) || any(which < 1) || any(which > 6)){
    cli::cli_abort(c("{.var which} must be a numeric vector with values
                     between 1 and 6.",
                     "i" = "The value{?s} specified {?was/were}
                           {as.character(which)}."))
  }

  #binomialLike <- family(model)$family == "binomial"

  # Decide which plots to show
  show <- rep(FALSE, 6)
  show[which] <- TRUE


  # Prepare the data for plotting
  plot_data$Obs <- seq_along(plot_data$.fitted)
  plot_data$Label <- names(residuals(model))
  plot_data$.qq <- qqnorm(plot_data$.stdresid, plot.it = F)$x

  # If weighted regression then omit observations with weights zero
  w <- weights(model)
  if (!is.null(w)) {
    plot_data <- plot_data %>%
      mutate(weights = w) %>%
      dplyr::filter(.data$weights!=0)
  }

  # Get at values to use for plots 5 and 6
  hii <- plot_data$.hat

  # Get slope and intercept for QQ plot
  qq_y <- quantile(plot_data$.stdresid, names = FALSE,
                   probs = c(0.25, 0.75), na.rm = TRUE)
  qq_x <- qnorm(c(0.25, 0.75))
  qq_slope <- diff(qq_y) / diff(qq_x)
  qq_intercept <- qq_y[[1L]] - qq_slope * qq_x[[1L]]

  # Handling how extreme points are shown
  if(npoints < 0 || npoints > nrow(plot_data)){
    cli::cli_abort(c("{.var npoints} should be an integer between 0 and
                      {nrow(plot_data)}",
                     "i" = "The specified value{?s} was
                            {as.character(npoints)}."))
  }

  # Store extreme points in data.frame
  if(any(show[1:3])){
    show.r <- plot_data %>%
                arrange(desc(abs(.data$.resid))) %>%
                slice(1:npoints)
  }
  if(any(show[4:6])){
    show.cook <- plot_data %>%
                   arrange(desc(abs(.data$.cooksd))) %>%
                   slice(1:npoints)
  }

  # Colours for the pie-glyph slices
  if(pies && is.null(pie_colours)){
    colours <- get_colours(model_species, FG = eval(model$DIcall$FG))
  } else if (pies && !is.null(pie_colours)){
    # sanity_checks(colours = pie_colours)
    if(length(pie_colours) != length(model_species)){
      cli::cli_abort(c("The number of {.var colours} specified should be
                       same as the number of columns specified in prop.",
                       "i" = "There are {length(model_species)} in prop
                              but only {length(colours)} colour{?s}
                              {?was/were} specified."))
    }
    colours <- pie_colours
  }

  # To adjust the positioning of hanging text for extreme observations
  pie_grob <- pieGrob(values = 1:4,
                      radius = ifelse(pies, pie_radius, points_size/15))

  # Progress bar to be shown if plots take long time
  p_bar <- cli_progress_bar(
    total = length(which),
    format = "Creating plot {pb_current} of {pb_total} |
             {pb_bar} {pb_percent} | ETA: {pb_eta}"
  )

  plots <- list()

  if(show[1L]){
    # Best fit line across points
    smoothing_line <- smoothing_fun(plot_data$.fitted,  plot_data$.resid)

    # Base plot with points and other aesthetics
    pl <- ggplot(plot_data, aes(.data$.fitted, .data$.resid))+
            theme_bw()+
            geom_point(size = 3)+
            geom_line(data = smoothing_line,
                      aes(x = .data$x, y = .data$y),
                      colour = "red", linetype = 1) +
            geom_hline(yintercept=0, col="black", linetype="dashed")

    # If possible then replace points with pie-chart glyphs showing proportions
    if(pies){
      # Show extreme observations with pie-chart glyphs
      if(only_extremes){
        pl <- add_pies(pl = pl, colours = colours,
                       data = show.r, slices = model_species,
                       colour = "black", radius = pie_radius)
      # Else show all points as pie-glyphs
      } else {
        pl <- add_pies(pl = pl, colours = colours,
                       slices = model_species, colour = "black",
                       radius = pie_radius)
      }
    }

    # Flag extreme observations
    pl <- add_label(pl = pl, data = show.r,
                    grob_obj = pie_grob,
                    label_size = label_size) +
      labs(x = "Fitted Values", y = "Residuals",
           title = "Residual vs Fitted")+
      scale_y_continuous(expand = expansion(mult = 0.1))

    plots[["ResivsFit"]] <- pl
    cli_progress_update(id = p_bar, set = sum(show[1L]))
  }

  if(show[2L]){
    pl <- ggplot(plot_data, aes(.data$.qq, .data$.stdresid))+
            theme_bw()+
            geom_point(size = 3, na.rm = TRUE)+
            geom_abline(slope = qq_slope, intercept = qq_intercept)

    if(pies){
      # Show extreme observations with pie-chart glyphs
      if(only_extremes){
        pl <- add_pies(pl = pl, colours = colours,
                       data = show.r, slices = model_species,
                       colour = "black", radius = pie_radius)
      # Else show all points as pie-glyphs
      } else {
        pl <- add_pies(pl = pl, colours = colours,
                       slices = model_species, colour = "black",
                       radius = pie_radius)
      }
    }

    pl <- add_label(pl = pl, data = show.r,
                    grob_obj = pie_grob,
                    label_size = label_size) +
      xlab("Theoretical Quantiles")+
      ylab("Std. Pearson Residuals")+
      scale_y_continuous(expand = expansion(mult = 0.1))+
      ggtitle("Normal Q-Q")

    plots[["QQplot"]] <- pl
    cli_progress_update(id = p_bar, set = sum(show[1L:2L]))
  }

  if(show[3L]){
    smoothing_line <- smoothing_fun(plot_data$.fitted,
                                    sqrt(abs(plot_data$.stdresid)))

    pl <- ggplot(plot_data, aes(.data$.fitted, sqrt(abs(.data$.stdresid))))+
            theme_bw()+
            geom_point(size = 3, na.rm=TRUE)+
            geom_line(data = smoothing_line,
                      aes(x = .data$x, y = .data$y),
                      colour = "red", linetype = 1)

    if(pies){
      # Show extreme observations with pie-chart glyphs
      if(only_extremes){
        pl <- add_pies(pl = pl, colours = colours,
                       data = show.r, slices = model_species,
                       colour = "black", radius = pie_radius)
      # Else show all points as pie-glyphs
      } else {
        pl <- add_pies(pl = pl, colours = colours,
                       slices = model_species, colour = "black",
                       radius = pie_radius)
      }
    }

    pl <- add_label(pl = pl, data = show.r,
                    grob_obj = pie_grob,
                    label_size = label_size) +
      xlab("Fitted Value")+
      ylab(expression(sqrt("|Std. Pearson residuals|")))+
      scale_y_continuous(expand = expansion(mult = 0.1))+
      ggtitle("Scale-Location")

    plots[["Scale-Location"]] <- pl
    cli_progress_update(id = p_bar, set = sum(show[1L:3L]))
  }

  if(show[4L]){
    pl <- ggplot(plot_data, aes(x = .data$Obs, y = .data$.cooksd))+
            theme_bw()+
            geom_point(size = 3)+
            geom_linerange(aes(ymax = .data$.cooksd, ymin = 0),
                           linewidth = 1)

    if(pies){
      # Show extreme observations with pie-chart glyphs
      if(only_extremes){
        pl <- add_pies(pl = pl, colours = colours,
                       data = show.cook, slices = model_species,
                       colour = "black", radius = pie_radius)
      # Else show all points as pie-glyphs
      } else {
        pl <- add_pies(pl = pl, colours = colours,
                       slices = model_species, colour = "black",
                       radius = pie_radius)
      }
    }

    pl <- add_label(pl = pl, data = show.cook,
                    grob_obj = pie_grob,
                    label_size = label_size) +
      xlab("Obs. Number")+
      ylab("Cook's distance")+
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
      ggtitle("Cook's distance")

    plots[["CooksD"]] <- pl
    cli_progress_update(id = p_bar, set = sum(show[1L:4L]))
  }

  if(show[5L]){

    r.hat <- range(hii, na.rm = TRUE)
    isConst.hat <- all(r.hat == 0) ||
                       diff(r.hat) < 1e-10 * mean(hii, na.rm = TRUE)
    y_thresh <- max(abs(plot_data$.stdresid)) + 0.5
    p <- model$rank
    xmax <- max(hii) + 0.01
    hh <- seq.int(min(r.hat[1L], r.hat[2L]/100), xmax, length.out = 101)
    cook_contours <- data.frame(x = vector(),
                                y = vector(),
                                group = vector(),
                                label = vector())

    if (length(cook.levels)) {
      for (i in seq_along(cook.levels)) {
        cl.h <- sqrt(cook.levels[i] * p * (1 - hh)/hh)
        x <- rep(hh, times = 2)
        y <- c(cl.h, -1 * cl.h)
        group <- rep(c(2*i - 1, 2*i), each = 101)
        label <- cook.levels[i]
        cook_contours <- bind_rows(cook_contours,
                                   data.frame(x, y, group, label))
      }
    }

    cook_contours <- cook_contours[abs(cook_contours$y) <= y_thresh, ]
    smoothing_line <- smoothing_fun(plot_data$.hat, plot_data$.stdresid)

    panel_plot <- ggplot(plot_data, aes(.data$.hat, .data$.stdresid))+
                    theme_bw()+
                    geom_point(size = 3, na.rm=TRUE)+
                    geom_line(data = smoothing_line,
                              aes(x = .data$x, y = .data$y),
                              colour = "red", linetype = 1)

    if(nrow(cook_contours) > 0){
      panel_plot <- panel_plot +
        geom_line(data = cook_contours,
                  aes(x = .data$x, y = .data$y, group = group),
                  colour = 'grey20', linetype = 2) +
        geom_text(data = cook_contours %>% filter(x == max(x)),
                  aes(x = .data$x + 0.01, y = .data$y, label = .data$label),
                  size = label_size, colour = 'grey20')
    }

    if(pies){
      # Show extreme observations with pie-chart glyphs
      if(only_extremes){
        panel_plot <- add_pies(pl = panel_plot, colours = colours,
                               data = show.cook, slices = model_species,
                               colour = "black", radius = pie_radius)
      # Else show all points as pie-glyphs
      } else {
        panel_plot <- add_pies(pl = panel_plot, colours = colours,
                               slices = model_species, colour = "black",
                               radius = pie_radius)
      }
    }

    panel_plot <- add_label(pl = panel_plot, data = show.cook,
                            grob_obj = pie_grob,
                            label_size = label_size,
                            lab_pos = c(0, 1.75, 0, 0)) +
      xlab("Leverage")+
      ylab("Std. Pearson Residuals")+
      ggtitle("Residual vs Leverage Plot")+
      scale_y_continuous(expand = expansion(mult = 0.1))+
      guides(radius = "none")+
      xlim(c(0, NA))

    plots[["Leverage"]] <- panel_plot
    cli_progress_update(id = p_bar, set = sum(show[1L:5L]))
  }

  if(show[6L]){
    p <- model$rank
    g <- dropInf(hii/(1 - hii), hii)
    bval <- pretty(sqrt(p * plot_data$.cooksd/g), 5)
    xmax <- max(plot_data$.hat)
    ymax <- max(plot_data$.cooksd)

    smoothing_line <- smoothing_fun(plot_data$.hat, plot_data$.cooksd)

    panel_plot <- ggplot(plot_data, aes(.data$.hat, .data$.cooksd)) +
                    theme_bw()+
                    geom_point(size = 3, na.rm=TRUE)+
                    geom_line(data = smoothing_line,
                              aes(x = .data$x, y = .data$y),
                              colour = "red", linetype = 1)

    abline_data <- data.frame("intercept" = vector(),
                              "slope" = vector())
    segment_data <- data.frame(x = vector(),
                               y = vector(),
                               xend = vector(),
                               yend = vector())

    for (i in seq_along(bval)) {
      bi2 <- bval[i]^2
      if (p * ymax > bi2 * xmax) {
        xi <- xmax + 1/3
        yi <- bi2 * xi/p
        abline_data <- bind_rows(abline_data,
                                 c("intercept" = 0, "slope" = bi2))

      } else {
        yi <- ymax + 0.01
        xi <- p * yi/bi2
        segment_data <- bind_rows(segment_data,
                                  c("x" = 0 ,"y" = 0,
                                    "xend" = xi, "yend" = yi))
      }
    }

    panel_plot <- panel_plot +
      geom_abline(data = abline_data,
                  aes(intercept = .data$intercept,
                      slope = .data$slope), linetype = 2)+
      geom_segment(data = segment_data,
                   aes(x = .data$x, y = .data$y,
                       xend = .data$xend, yend = .data$yend), linetype = 2)
    if(pies){
      # Show extreme observations with pie-chart glyphs
      if(only_extremes){
        panel_plot <- add_pies(pl = panel_plot, colours = colours,
                               data = show.cook, slices = model_species,
                               colour = "black", radius = pie_radius)
      # Else show all points as pie-glyphs
      } else {
        panel_plot <- add_pies(pl = panel_plot, colours = colours,
                               slices = model_species, colour = "black",
                               radius = pie_radius)
      }
    }

    panel_plot <- add_label(pl = panel_plot, data = show.cook,
                            grob_obj = pie_grob,
                            label_size = label_size,
                            lab_pos = c(0, 1.75, 0, 0)) +
      xlab(str2expression("Leverage~h[ii]"))+
      ylab("Cook's Distance")+
      xlim(c(0, NA))+
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
      ggtitle(expression("Cook's dist vs Leverage* " * h[ii]/(1 - h[ii])))

    plots[["CooksvsLev"]] <- panel_plot
    cli_progress_update(id = p_bar, set = sum(show[1L:6L]))
  }

  if(length(plots) > 1){
    plot <- new("ggmultiplot", plots = plots, nrow = nrow, ncol = ncol)
  } else {
    plot <- plots[[1]]
  }
  cli_process_done(id = p_bar)
  cli::cli_alert_success("Created all plots.")
  return(plot)
}

# Helper functions to add pies
add_pies <- function(pl, colours, legend_position = "top", ...){
  pl +
    geom_pie_glyph(...)+
    theme(legend.position = legend_position)+
    scale_fill_manual(values = colours)
}

# Helper function to labels at fixed positions
add_label <- function(pl, data, grob_obj, label_size = 4,
                      lab_pos = c(0, 0, 1.75, 0)){
  pl + geom_richtext(data = data,
                     aes(label = .data$Label),
                     fill = NA, label.colour = NA, size = label_size,
                     label.margin = unit(lab_pos, "grobheight",
                                         data = grob_obj))
}
