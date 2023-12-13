#' @title DImodelsVis: Model interpretation and visualisation for compositional data
#'
#' @description Statistical models fit to compositional data are often difficult to
#' interpret due to the sum to 1 constraint on data variables. `DImodelsVis` provides
#' novel visualisations tools to aid with the interpretation of models fit to
#' compositional data. All visualisations in the package are created using the
#' `ggplot2` plotting framework and can be extended like every other ggplot object.
#'
#' @details
#' \strong{Introduction to Diversity-Interactions models (this feels out of place now?)} \cr
#' Diversity-Interactions (DI) models (Kirwan et al 2009) are a set of tools for analysing and interpreting data from experiments that explore the effects of species diversity on community-level responses; for example, the effect of increasing community diversity on biomass production in a grassland ecosystem. Most analyses of diversity experiments quantify community diversity in terms of species richness, the number of species present. The DI method modifies this presence/absence approach in mixed communities by taking species relative abundance in the community into account. So, instead of ignoring differences in community responses across communities with the same species but with different species relative abundances, the DI approach aims to understand and explain these differences.
#'
#' @author
#' \strong{Maintainter:} Rishabh Vishwakarma \email{vishwakr@@tcd.ie} (\href{https://orcid.org/0000-0002-4847-3494}{ORCID})
#'
#' \strong{Authors:} \cr
#' \itemize{
#'     \item{Caroline Brophy}
#'     \item{Catherine Hurley}
#' }
#'
#' @references
#' \itemize{
#'     \item{Moral, R.A., Vishwakarma, R., Connolly, J., Byrne, L., Hurley, C., Finn, J.A. and Brophy, C., 2023. Going beyond richness: Modelling the BEF relationship using species identity, evenness, richness and species interactions via the DImodels R package. Methods in Ecology and Evolution, 14(9), pp.2250-2258.}
#'     \item{Brophy C, A Dooley, L Kirwan, JA Finn, J McDonnell, T Bell, MW Cadotte and J Connolly. (2017) Biodiversity and ecosystem function: Making sense of numerous species interactions in multi-species communities. Ecology 98, 1771-1778.}
#'     \item{Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A Lüscher, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355.}
#'     \item{Kirwan L, J Connolly, JA Finn, C Brophy, A Lüscher, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038.}
#' }
#'
#'
#' @examples
#' ## Load libraries
#' library(DImodels)
#' library(DImodelsVis)
#'
#' ## Load data
#' data(sim2)
#' sim2 <- sim2[sim2$block == 1, ]
#'
#' ## Fit model with compositional data
#' mod <- DI(y = "response", prop = 3:6,
#'           DImodel = "AV", data = sim2)
#'
#' ## Model diagnostics plots but points are replaced by
#' ## pie-glyphs showing the proportions of the compositional variables
#' ## See \code{\link{model_diagnostics}} for more information
#' model_diagnostics(model = mod)
#'
#' ## Visualise the predicted response variable as contributions
#' ## (\eqn{\beta} * predictor value) from the individual
#' ## terms in the model
#' ## See \code{\link{prediction_contributions}} for more information
#' prediction_contributions(model = mod)
#'
#' ## Visualise the change in average response over a diversity gradient
#' ## This plot shows the change in the response over a diversity gradient
#' ## We use richness (number of non-zero variables in a given observation)
#' ## as our gradient in this example. The black line shows the average response
#' ## at each level of richness while the position of the
#' ## `\code{\link[PieGlyph:PieGlyph-package]{pie-glyphs}}` show variations
#' ## about this average whilst also showing the relative abundances of each
#' ## variable in the composition.
#' ## See \code{\link{gradient_change}} for more information
#' plot_data <- get_equi_comms(nvars = 4, variables = c("p1", "p2", "p3", "p4"))
#' gradient_change(model = mod, data = plot_data)
#'
#' ## Visualise effects of increasing or decreasing a variable
#' ## within a set of compositional variables
#' ## This plot shows the effect of increasing the proportion of p1
#' ## in several different initial compositions of the variables
#' ## p1, p2, p3, and p4. Each curve shows the effect of increasing
#' ## the proportion of p1 whilst keeping the relative proportions of
#' ## the other three variables unchanged
#' ## See \code{\link{visualise_effects}} for more information
#' visualise_effects(model = mod, var_interest = "p1")
#'
#' ## Visualise slices of the n-dimensional simplex as ternary diagrams.
#' ## 2-d slices of the n-dimensional simplex are created by conditioning
#' ## certain compositional variables at a specific value `x` while the
#' ## remaining variables are allowed to vary within the range `0` to `1-x`.
#' ## In this example variable p1 is conditioned to have values `0`, `0.2`, and `0.5`
#' ## One ternary diagram is created for each case where p2, p3, and p4 are
#' ## allowed to vary from `0` upto `1`, `0.8`, and `0.5`, respectively.
#' ## This is equivalent to taking multiple slices of the n-dimensional simplex
#' ## and viewing multiple slices would enable us to get a picture the change
#' ## in the response across the n-dimensional simplex.
#' ## For example the response is maximised where p1 is 0.2
#' ## See \code{\link{conditional_ternary}} for more information
#' conditional_ternary(model = mod, tern_vars = c("p2", "p3", "p4"),
#'                     conditional = list("p1" = c(0, 0.2, 0.5)),
#'                     resolution = 1)
#'
#'
#' @keywords internal
"_PACKAGE"
NULL
