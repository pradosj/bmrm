
#' @name gradient
#' @rdname gradient
#' @aliases gradient.default
#' @aliases gradient<-
#' @aliases gradient<-.default
#' @title Return or set gradient attribute
#' @param x any R object
#' @param value new gradient value to set
#' @return attr(x,"gradient")
#' @export
gradient <- function(x,...) UseMethod("gradient")

#' @rdname gradient
#' @export
gradient.default <- function(x) attr(x, "gradient")

#' @rdname gradient
#' @export
"gradient<-" <- function(x,value,...) UseMethod("gradient<-")

#' @rdname gradient
#' @export
"gradient<-.default" <- function(x,value) {attr(x, "gradient") <- value;x}

