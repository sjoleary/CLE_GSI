get.num.gens <- function(m,N){
  # Calculates half life to equilibrium and multiplies it by 50
  # Whitlock ref
    L <- (1-m)^2*(1-1/N)
		t.half <- log(1/2)/log(L)
		gens <- 50*t.half
    return(round(gens))
}