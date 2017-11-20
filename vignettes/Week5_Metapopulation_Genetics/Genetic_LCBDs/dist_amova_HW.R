dist_amova_HW <- function (x) 
{
  if (is(x, "locus")) 
    x <- data.frame(LOCUS = x)
  if (!is(x, "data.frame")) 
    stop(paste("The function dist_amova() requires a locus vector or a data frame of locus vectors.  You passed a '", 
               class(x), "' object.", sep = ""))
  N <- dim(x)[1]
  ret <- matrix(0, ncol = N, nrow = N)
  p <- ploidy(x)
  if (any(p$Ploidy != round(p$Ploidy))) {
    data <- FALSE
    if (any(p$Ploidy > 2)) 
      stop("As currently implemented, the 2gener amova distance is limited to diploid individuals.")
    loci <- p$Locus
    for (locus_name in loci) {
      locus <- x[[locus_name]]
      freqs <- frequencies(locus)
      y <- to_mv(locus)
      p <- ploidy(locus)
      for (i in seq(1:N)[p == 2]) {
        f <- freqs$Frequency * y[i, ]
        y[i, ] <- f/sum(f)/2
      }
      if (is(data, "logical")) 
        data <- y
      else data <- cbind(data, y)
    }
  }
  else data <- to_mv(x, drop.allele = FALSE)
  for (i in 1:N) {
    x <- data[i, ]
    for (j in 1:i) {
      if (i != j) {
        y <- data[j, ]
        ret[i, j] <- ret[j, i] <- sum(2 * t(x - y) %*% 
                                        (x - y))
      }
    }
  }
  res <- list(dist=ret, data=data)
  return(res)
}