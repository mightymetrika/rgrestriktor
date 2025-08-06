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
    anova_tab <- stats::anova(mod) # get anova results
    # get pairwise t-test results
    pttest <- rstatix::pairwise_t_test(data = dat,
                                       formula = .formula,
                                       p.adjust.method = "bonferroni") |>
      dplyr::mutate(modname = "pairwise_t_test") |>
      dplyr::select(modname, p, p.adj, group1, group2, n1, n2)

    # use rgr_iht to obtain iht results
    rgr_iht_res <- rgr_iht(mod,
                   rgcontraints = rgcontraints,
                   constraints = constraints)

    # Combine anova, pairwise t-test, and iht_df dataframes
    rgr_iht_res$iht_df <- dplyr::bind_rows(data.frame(modname = "anova",
                                                      statistic = anova_tab$`F value`[1],
                                                      p = anova_tab$`Pr(>F)`[1]),
                                           as.data.frame(pttest),
                                           rgr_iht_res$iht_df) |>
      dplyr::mutate(N = nrow(dat)) |>
      dplyr::select(N, dplyr::everything())

    # store results
    res[[i]] <- rgr_iht_res
  }

  return(res)


}
