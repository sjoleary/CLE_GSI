myexport <- function(...) {
  arg.list <- list(...)
  names <- all.names(match.call())[-1]
  for (i in seq_along(names)) assign(names[i],arg.list[[i]],.GlobalEnv)
}