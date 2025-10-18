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

    dat <- dgp()                                   # get data

    # (1) IHT model: no intercept so means are parameters
    mod_iht <- stats::lm(.formula, data = dat)

    # (2) Classical ANOVA & pairwise t-tests should use an intercepted formula
    outcome   <- all.vars(.formula)[1]
    rhs_terms <- attr(stats::terms(.formula), "term.labels")
    if (length(rhs_terms) != 1L) stop("Expecting a single grouping factor.")
    group_var <- rhs_terms[1]
    anv_formula <- stats::as.formula(paste(outcome, "~", group_var))
    mod_anv <- stats::lm(anv_formula, data = dat)
    anova_tab <- stats::anova(mod_anv)            # standard omnibus F for y ~ group

    # pairwise t-tests on the intercepted formula
    pttest <- rstatix::pairwise_t_test(
      data = dat, formula = anv_formula, p.adjust.method = "bonferroni"
    ) |>
      dplyr::mutate(modname = "pairwise_t_test") |>
      dplyr::select(modname, p, p.adj, group1, group2, n1, n2)

    # run IHT on the no-intercept model
    rgr_iht_res <- rgr_iht(mod_iht, rgcontraints = rgcontraints, constraints = constraints)

    # combine
    rgr_iht_res$iht_df <- dplyr::bind_rows(
      data.frame(modname = "anova", statistic = anova_tab$`F value`[1],
                 p = anova_tab$`Pr(>F)`[1]),
      as.data.frame(pttest),
      rgr_iht_res$iht_df
    ) |>
      dplyr::mutate(N = nrow(dat)) |>
      dplyr::select(N, dplyr::everything())

    # store results
    res[[i]] <- rgr_iht_res
  }

  return(res)


}
