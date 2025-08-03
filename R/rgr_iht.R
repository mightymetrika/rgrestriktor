rgr_iht <- function(model = NULL, data = NULL, rgcontraints = NULL, constraints = NULL){

  iht <- restriktor::iht(model, constraints = constraints)

  print(iht)


  rg_iht <- vector(mode = "list", length = length(rgcontraints))
  for (i in 1:length(rgcontraints)){
    rgc <- rgcontraints[i]

    rg_iht[i] <- restriktor::iht(model, constraints = rgc)

  }

  res <- list(iht = iht, rg_iht = rg_iht)
  class(res) <- "rgr_iht"

  return(res)

}
