ConditionalR2 <- function(model) {
  var_fixed <- var(as.numeric(fixef(model) %*% t(model@pp$X)))
  var_random <- sum(sapply(VarCorr(model), function(x) attr(x, "stddev")^2))
  var_residual <- attr(VarCorr(model), "sc")^2
  r2_conditional <- (var_fixed + var_random) / (var_fixed + var_random + var_residual)

  return(r2_conditional)
}
