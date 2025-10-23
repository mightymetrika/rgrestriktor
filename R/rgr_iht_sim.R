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

    # ---- NEW: enrich pairwise rows with direction & adjacency info ----
    outcome   <- all.vars(.formula)[1]
    rhs_terms <- attr(terms(.formula), "term.labels")
    group_var <- rhs_terms[1]

    # group means (for sign of difference)
    means_df  <- dat |>
      dplyr::group_by(!!rlang::sym(group_var)) |>
      dplyr::summarise(.mean = mean(.data[[outcome]]), .groups = "drop") |>
      dplyr::rename(group = !!group_var) |>
      dplyr::mutate(group = as.character(group))

    # Try to infer intended order. If groups look numeric, use that numeric order.
    levs <- means_df$group
    if (all(suppressWarnings(!is.na(as.numeric(levs))))) {
      ord <- levs[order(as.numeric(levs))]
    } else {
      # fallback: alphabetical order of labels
      ord <- sort(levs)
    }
    # adjacency map: consecutive elements in 'ord' are adjacent
    adj_pairs <- tibble::tibble(
      group1 = head(ord, -1),
      group2 = tail(ord, -1)
    )

    pttest <- pttest |>
      dplyr::mutate(group1 = as.character(group1),
                    group2 = as.character(group2)) |>
      dplyr::left_join(dplyr::rename(means_df, group1 = group, mean1 = .mean), by = "group1") |>
      dplyr::left_join(dplyr::rename(means_df, group2 = group, mean2 = .mean), by = "group2") |>
      dplyr::mutate(
        diff   = mean2 - mean1,
        dir_ok = diff > 0
      ) |>
      dplyr::left_join(dplyr::mutate(adj_pairs, adj = TRUE),
                       by = c("group1","group2")) |>
      dplyr::mutate(adj = dplyr::coalesce(adj, FALSE))

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
