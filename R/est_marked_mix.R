est_marked_mix <- function(pp) {
  # test if this pp has mark
  if(missing(pp$marks)) {
    stop("This point pattern doesn't include marks.")
  }

}
