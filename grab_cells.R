grab_cells <- function(object, feature = NULL) {
  if (is.numeric(feature)) {
    return(names(object@ident[object@ident == feature]))
  } 
}