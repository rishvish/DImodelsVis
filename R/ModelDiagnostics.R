#' @title Regression diagnostics plots
#' @description This function returns regression diagnostics plots for a Diversity Interactions (DI) model with points replaced by pie-chart glyphs making it easier to track various communities in the plots
#'
#' @importFrom PieGlyph geom_pie_glyph
#' @importFrom dplyr bind_rows desc filter arrange slice
#' @importFrom stats family residuals qqnorm weights quantile qnorm
#' @importFrom ggplot2 fortify xlim ylim xlab ylab ggtitle geom_hline geom_text geom_abline geom_linerange geom_segment
#' @importFrom methods new
#' @importClassesFrom ggfortify ggmultiplot
#'
#' @param model A Diversity Interactions model object fit by using the \code{\link[DImodels:DI]{DI()}} function from the \code{\link[DImodels:DImodels-package]{DImodels}} package.
#' @param which a subset of the numbers 1 to 6, by default 1, 2, 3, and 5, referring to
#'                  1. - "Residuals vs Fitted", aka ‘Tukey-Anscombe’ plot
#'                  2. - "Normal Q-Q" plot, an enhanced qqnorm(resid(.))
#'                  3. - "Scale-Location"
#'                  4. - "Cook's distance"
#'                  5. - "Residuals vs Leverage"
#'                  6. - "Cook's dist vs Lev./(1-Lev.)"
#' @param npoints number of points to be labelled in each plot, starting with the most extreme.
#' @param cook.levels levels of Cook's distance at which to draw contours.
#' @param title Title for the plot
#' @param nrow 	Number of rows in which to arrange the final plot
#' @param ncol 	Number of columns in which to arrange the final plot
#'
#' @return A ggmultiplot class object
#' @export
#'
#' @examples
#' library(DImodels)
#'
#' ## Load data
#' data(sim1)
#'
#' ## Fit DI model
#' mod1 <- DI(y = "response", prop = 3:6, data = sim1, DImodel = "ADD")
#'
#' ## Get diagnostics plot
#' diagnostics <- model_diagnostics(mod1)
#' print(diagnostics)
#'
#' ## Get all plots
#' diagnostics <- model_diagnostics(mod1, which = 1:6)
#'
#' ## Access individual plots
#' print(diagnostics[[1]])
#' print(diagnostics[[6]])
#'
#' ## Change plot arrangement
#' model_diagnostics(mod1, which = c(1, 3), nrow = 2, ncol = 1)
#'
#' ## Modify contour levels for residual vs leverage plot
#' model_diagnostics(mod1, which = 5, cook.levels = c(0.02, 0.04, 0.06))
#'
#' ## Change number of extreme points to be flagged
#' model_diagnostics(mod1, which = 2, npoints = 5)
model_diagnostics <- function(model, which = c(1,2,3,5), npoints = 3,
                              cook.levels = c(0.5, 1),
                              title = "Diagnostic Plots",
                              nrow = 0, ncol = 0){
  # Sanity checks
  # Ensure model is a DImodels object
  if(!inherits(model, "DI")){
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

  # Sanity checks from plot.lm function
  if (!is.numeric(which) || any(which < 1) || any(which > 6)){
    stop("'which' must be a numeric vector with values between 1 and 6")
  }

  binomialLike <- family(model)$family == "binomial"

  # Decide which plots to show
  show <- rep(FALSE, 6)
  show[which] <- TRUE


  # Prepare the data for plotting
  plot_data <- fortify(model = model)
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
  if (!is.integer(npoints) && length(npoints) != 1){
    stop("'npoints' should be an integer describing the number of extreme points to be highlighed")
  } else if(npoints < 0 || npoints > nrow(plot_data)){
    stop(glue::glue("'npoints' should be between 0 and {nrow(plot_data)}"))
  }

  # Store extreme points in data.frame
  if(any(show[1:3])){
    show.r <- plot_data %>% arrange(desc(abs(.data$.resid))) %>% slice(1:npoints)
  }
  if(any(show[4:6])){
    show.cook <- plot_data %>% arrange(desc(abs(.data$.cooksd))) %>% slice(1:npoints)
  }

  plots <- list()

  if(show[1L]){
    smoothing_line <- smoothing_fun(plot_data$.fitted,  plot_data$.resid)

    plots[["ResivsFit"]] <- ggplot(plot_data, aes(.data$.fitted, .data$.resid))+
      geom_point()+
      geom_line(data = smoothing_line, aes(x = .data$x, y = .data$y),
                colour = "red", linetype = 1) +
      geom_hline(yintercept=0, col="black", linetype="dashed")+
      geom_pie_glyph(slices = model_species, colour = "black")+
      geom_text(data = show.r, aes(label = .data$Label, y = .data$.resid + 0.2))+
      labs(x = "Fitted Values", y = "Residuals",
           title = "Residual vs Fitted")+
      theme_bw()+
      theme(legend.position = "top")
  }

  if(show[2L]){
    plots[["QQplot"]] <- ggplot(plot_data, aes(.data$.qq, .data$.stdresid))+
      geom_point(na.rm = TRUE)+
      geom_abline(slope = qq_slope, intercept = qq_intercept)+
      geom_pie_glyph(slices = model_species, colour = "black")+
      geom_text(data = show.r, aes(label = .data$Label, y = .data$.stdresid + 0.2))+
      xlab("Theoretical Quantiles")+
      ylab("Std. Pearson Residuals")+
      ggtitle("Normal Q-Q")+
      theme_bw()+
      theme(legend.position = "top")
  }

  if(show[3L]){
    smoothing_line <- smoothing_fun(plot_data$.fitted,  sqrt(abs(plot_data$.stdresid)))

    plots[["Scale-Location"]] <- ggplot(plot_data, aes(.data$.fitted, sqrt(abs(.data$.stdresid))))+
      geom_point(na.rm=TRUE)+
      geom_line(data = smoothing_line, aes(x = .data$x, y = .data$y),
                colour = "red", linetype = 1) +
      geom_pie_glyph(slices = model_species, colour = "black")+
      geom_text(data = show.r, aes(label = .data$Label, y = sqrt(abs(.data$.stdresid)) + 0.2))+
      xlab("Fitted Value")+
      ylab(expression(sqrt("|Std. Pearson residuals|")))+
      ggtitle("Scale-Location")+
      theme_bw()+
      theme(legend.position = "top")
  }

  if(show[4L]){
    plots[["CooksD"]] <- ggplot(plot_data, aes(x = .data$Obs, y = .data$.cooksd))+
      geom_linerange(aes(ymax = .data$.cooksd, ymin = 0), linewidth = 1)+
      geom_pie_glyph(slices = model_species, colour = "black")+
      geom_text(data = show.cook, aes(label = .data$Label, y = .data$.cooksd + 0.0025))+
      xlab("Obs. Number")+
      ylab("Cook's distance")+
      ggtitle("Cook's distance")+
      theme_bw()+
      theme(legend.position = "top")
  }

  if(show[5L]){

    r.hat <- range(hii, na.rm = TRUE)
    isConst.hat <- all(r.hat == 0) || diff(r.hat) < 1e-10 * mean(hii, na.rm = TRUE)
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
        cook_contours <- bind_rows(cook_contours, data.frame(x, y, group, label))
      }
    }

    cook_contours <- cook_contours[abs(cook_contours$y) <= y_thresh, ]
    smoothing_line <- smoothing_fun(plot_data$.hat, plot_data$.stdresid)

    panel_plot <- ggplot(plot_data, aes(.data$.hat, .data$.stdresid))+
      geom_point(na.rm=TRUE)+
      geom_line(data = smoothing_line, aes(x = .data$x, y = .data$y),
                colour = "red", linetype = 1) +
      geom_pie_glyph(slices = model_species, colour = "black")+
      geom_text(data = show.cook, aes(label = .data$Label, x = .data$.hat - 0.01))+
      xlab("Leverage")+
      ylab("Std. Pearson Residuals")+
      ggtitle("Residual vs Leverage Plot")+
      #scale_size_continuous("Cook's Distance", range=c(1,5))+
      guides(radius = "none")+
      xlim(c(0, NA))+
      theme_bw()+
      theme(legend.position="top")

    if(nrow(cook_contours) > 0){
      panel_plot <- panel_plot +
        geom_line(data = cook_contours, aes(x = .data$x, y = .data$y, group = group),
                  colour = 'grey20', linetype = 2) +
        geom_text(data = cook_contours %>% filter(x == max(x)),
                  aes(x = .data$x + 0.01, y = .data$y, label = .data$label), colour = 'grey20')
    }

    plots[["Leverage"]] <- panel_plot
  }

  if(show[6L]){
    p <- model$rank
    g <- dropInf(hii/(1 - hii), hii)
    bval <- pretty(sqrt(p * plot_data$.cooksd/g), 5)
    xmax <- max(plot_data$.hat)
    ymax <- max(plot_data$.cooksd)

    smoothing_line <- smoothing_fun(plot_data$.hat, plot_data$.cooksd)

    panel_plot <- ggplot(plot_data, aes(.data$.hat, .data$.cooksd))+
      geom_point(na.rm=TRUE)+
      geom_line(data = smoothing_line, aes(x = .data$x, y = .data$y),
                colour = "red", linetype = 1) +
      geom_pie_glyph(slices = model_species, colour = "black")+
      geom_text(data = show.cook, aes(label = .data$Label, x = .data$.hat - 0.01))+
      xlab(str2expression("Leverage~h[ii]"))+
      ylab("Cook's Distance")+
      xlim(c(0, NA))+
      ggtitle(expression("Cook's dist vs Leverage* " * h[ii]/(1 - h[ii])))+
      theme_bw()+
      theme(legend.position = 'top')

    abline_data <- data.frame('intercept' = vector(),
                              'slope' = vector())
    segment_data <- data.frame(x = vector(),
                               y = vector(),
                               xend = vector(),
                               yend = vector())

    for (i in seq_along(bval)) {
      bi2 <- bval[i]^2
      if (p * ymax > bi2 * xmax) {
        xi <- xmax + 1/3
        yi <- bi2 * xi/p
        abline_data <- bind_rows(abline_data, c('intercept' = 0, 'slope' = bi2))

      } else {
        yi <- ymax + 0.01
        xi <- p * yi/bi2
        segment_data <- bind_rows(segment_data,
                                  c('x' = 0 ,'y' = 0, 'xend' = xi, 'yend' = yi))
      }
    }

    panel_plot <- panel_plot +
      geom_abline(data = abline_data,
                  aes(intercept = .data$intercept,
                      slope = .data$slope), linetype = 2)+
      geom_segment(data = segment_data,
                   aes(x = .data$x, y = .data$y,
                       xend = .data$xend, yend = .data$yend), linetype = 2)

    plots[["CooksvsLev"]] <- panel_plot
  }
  new("ggmultiplot", plots = plots, nrow = nrow, ncol = ncol)
}


