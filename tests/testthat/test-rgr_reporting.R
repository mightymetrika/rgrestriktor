test_that("multiplication works", {

  iterations <- 1000
  k <- 5; levs <- as.character(seq_len(k))

  monotone <- list(
    name = "Monotone delta=1, n=5",
    rgconstraints = list(' group1 < group2 < group3 ',
                         ' group3 < group4 < group5 '),
    constraints  = ' group1 < group2 < group3 < group4 < group5 ',
    dgp = function() {
      data.frame(
        y = unlist(Map(function(mu) rnorm(5, mu, 1), 0:(k-1))),
        group = factor(rep(levs, each = 5), levels = levs)
      )
    },
    formula = "y ~ -1 + group"
  )

  middle_inversion <- list(
    name = "Middle inversion delta=1, n=5",
    rgconstraints = monotone$rgconstraints,
    constraints   = monotone$constraints,
    dgp = function() {
      # means: 0,3,2,1,4
      mu <- c(0,3,2,1,4)
      data.frame(
        y = unlist(Map(function(m) rnorm(5, m, 1), mu)),
        group = factor(rep(levs, each = 5), levels = levs)
      )
    },
    formula = "y ~ -1 + group"
  )

  summ <- rgr_run_scenarios(
    scenarios = list(monotone, middle_inversion),
    Its = iterations,
    include_strict = TRUE,
    include_anova_gate = TRUE,
    gate_args = list(
      require_omnibus = TRUE,
      segments = list(c("1","2","3"), c("3","4","5")),
      adj_method = "holm",
      alternative = c("greater","two.sided")
    )
  )


  abbr <- c("Monotone delta=1, n=5" = "MON",
            "Middle inversion delta=1, n=5" = "MID")

  cap <- "MON = Monotone delta=1, n=5;  MID = Middle inversion delta=1, n=5"

  # Table: group by scenario, order MON then MID
  rgp <- rgr_gt_power(
    summ,
    group_by = "scenario",
    scenario_labels = abbr,
    order_scenarios = c("MON","MID"),
    notes = cap,
    show_scenario = FALSE,
    separate = "scenario"
  )

  # Plot: facet by both family and scenario (MON left, MID right)
  rpp <- rgr_plot_power(
    summ, Its = iterations, caption = cap,
    facet_by = "both",
    scenario_labels = abbr,
    order_scenarios = c("MON","MID"),
    strip_right = FALSE,
    separate = "scenario"
  )


  expect_equal(2 * 2, 4)
})
