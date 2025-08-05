rgr_iht_sim <- function(Its = 1000,
                        rgcontraints = NULL,
                        constraints = NULL,
                        dgp,
                        .formula){

  # Simulation set up
  res <- vector(mode = "list", length = Its) # store results
  .formula <- stats::as.formula(.formula) # ensure formula is a valid formula


  # Set up main simulation loop
  for (i in 1:Its){
    dat <- dgp() # get data
    mod <- stats::lm(.formula, data = dat) #fit model

    # use rgr_iht to obtain iht results
    rgr_iht_res <- rgr_iht(mod,
                   rgcontraints = rgcontraints,
                   constraints = constraints)

    # store results
    res[[i]] <- rgr_iht_res
  }

  return(res)


}
