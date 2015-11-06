#' Plot membership indicator
#'
#' @inheritParams plot_avgsurf
#' @author Yuchen Wang, Jiaxun Chen
#'
#' @export
#' @examples
#' fit <- sppmix::est_mix_damcmc(pp = redwood, m = 3, truncate = FALSE,
#'                               L = 50000, LL = 100)
#' plot_ind(fit)

plot_ind <- function(fit, burnin = length(fit$allgens_List) / 10) {
  m <- dim(fit$genmus)[1]
  L <- length(fit$allgens_List)
  zs <- GetAvgLabelsDiscrete2Multinomial_sppmix(fit$genzs[(burnin + 1):L, ], m)
  plot_df <- tidyr::gather(data.frame(zs, point = 1:nrow(zs)),
                           comp, probability, -point)
  plot_df$component <- as.integer(gsub("X", "", plot_df$comp)) - 0.5

  ggplot2::qplot(point, component, data = plot_df, geom = "segment",
                 col = probability, xend = point, yend = component + 1,
                 size = I(5)) +
    ggplot2::scale_color_gradient(low = "white", high = "black") +
    ggplot2::coord_cartesian(ylim = c(.5, m + .5), xlim = c(0, nrow(zs) + 1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, size = 2)) +
    ggplot2::scale_y_discrete() +
    ggplot2::ggtitle("Plot of membership indicator")
}
