VarType <- function (values,
                     variable,
                     typeCategorical,
                     typeContinuous
                     ) {

  if (is.numeric(values) &&
      all(values %in% c(0, 1)) &&
      !(!is.null(typeContinuous) &&
        variable %in% typeContinuous) &&
      !(!is.null(typeCategorical) &&
        variable %in% typeCategorical)) {
    return("binary")
  }
  else if ((is.numeric(values) ||
            (!is.null(typeContinuous) &&
             variable %in% typeContinuous)) &&
           !(!is.null(typeCategorical) &&
             variable %in% typeCategorical)) {
    return("continuous")
  }
  else {
    return("categorical")
  }
}
