

#' @export
gradient <- function(x,...) UseMethod("gradient")

#' @export
gradient.default <- function(x) attr(x, "gradient")

#' @export
"gradient<-" <- function(x,...) UseMethod("gradient<-")

#' @export
"gradient<-.default" <- function(x,value) {attr(x, "gradient") <- value;x}

