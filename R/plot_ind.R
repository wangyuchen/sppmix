#' Plot membership indicator
#'
#' @param dares list of results from DAMCMC
#'
#'
#' @export
#' @examples
#' mix1 <- rnormmix(3, .01, 4, square(1))
#' pp1=rsppmix(200, mix1, spatstat::square(1))
#' data=cbind(pp1$x,pp1$y)
#' post=DAMCMC2d_sppmix(data,c(0,1),c(0,1),3, 5000, 1000, 50, trunc = FALSE)
#' plot_ind(post)

plot_ind <- function(dares) {
  data <- dares$meanz
  ind_df <- data.frame(point = 1:nrow(data),
                       ind = apply(data, 1, which.max) - .5)
  ggplot2::qplot(point, ind, data = ind_df, geom = "segment",
                 xend = point, yend = ind + 1, size = I(1.5)) +
    ggplot2::coord_cartesian(ylim = c(.5, ncol(data) + .5)) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, size = 2)) +
    ggplot2::ggtitle("Plot of membership indicator")
}
