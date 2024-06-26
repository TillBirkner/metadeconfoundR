Based on feedback on submission of metadeconfoundR v1.0.0, the following changes were made in v1.0.1:
  - removal of any print() statements, or rewrite into futile.logger statements, that can easily be suppressed
  - removal of all statements writing to file by default
  
results from "devtools::check(args = c('--as-cran'), build_args = c('--resave-data'))"
── R CMD check results ──────────────────────────── metadeconfoundR 1.0.1 ────
Duration: 58.7s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded
