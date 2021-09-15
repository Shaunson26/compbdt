#' Check zeros in a vector
#'
#'check_zero(c(1,1), 'all')
#'check_zero(c(1,1), 'any')
#'check_zero(c(0,1), 'any')
#'check_zero(c(0,1), 'all')
#'check_zero(c(0,0), 'all')
#'
#' @export
check_zero <- function(x, type = c('any', 'all')){

  type <- type[1]

  if (type == 'any') {
    stopifnot('A value of x is zero' = !any(x == 0))
  } else if (type == 'all') {
    stopifnot("All x are zero" = !all(x == 0))
  }

}
