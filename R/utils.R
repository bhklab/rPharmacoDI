#' Prepend zeros to an integer until its string representation has `string_length`
#'   characters.
#' 
#' @param x `integer(1)` An integer to prepend zeros to so the returned string
#'   length of `string_length`.
#' @param string_length `integer(1)` The desired length of the returned string.
#'
#' @export
zeroPad <- function(x, string_length) {
    stopifnot(
        is.numeric(x), 
        is.numeric(string_length),
        length(x) == 1,
        length(string_length) == 1
    )
    num_char <- nchar(as.character(x))
    padding <- string_length - num_char
    if (padding < 0) stop("The `x` argument has more characters than
        `string_length`!")
    if (padding > 0) return(
        paste0(
            paste0(rep('0', padding), collapse=''), 
            x
        ))
    else return(as.character(x))
}

#'
#'
#'
#' @export
zeroPadVector <- function(x, string_length) 
    vapply(x, FUN=zeroPad, string_length=string_length, FUN.VALUE=character(1))