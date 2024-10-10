ConditionalR2 <- function(model) {
  var_fixed <- var(as.numeric(lme4::fixef(model) %*% t(model@pp$X)))
  var_random <- sum(sapply(lme4::VarCorr(model), function(x) attr(x, "stddev")^2))
  var_residual <- attr(lme4::VarCorr(model), "sc")^2
  r2_conditional <- (var_fixed + var_random) / (var_fixed + var_random + var_residual)

  return(r2_conditional)
}
