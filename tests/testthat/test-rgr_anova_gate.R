test_that("multiplication works", {


  iterations <- 1000

  res1 <-  rgr_iht_sim(Its = iterations,
                       rgcontraints = list(' group1 < group2 < group3 ',
                                           ' group3 < group4 < group5 '),
                       constraints = ' group1 < group2 < group3 < group4 < group5 ',
                       dgp = function(x){
                         data.frame(
                           y = c(stats::rnorm(5),
                                 stats::rnorm(5, 1),
                                 stats::rnorm(5, 2),
                                 stats::rnorm(5, 3),
                                 stats::rnorm(5, 4)),
                           # x = stats::rnorm(25, 5),
                           group = as.factor(rep(c("1", "2", "3", "4", "5"), each = 5)))
                       },
                       .formula = "y ~ -1 + group")

  gate1 <- rgr_anova_gate_summarize(
    res1,
    alpha = 0.05,
    require_omnibus = TRUE,
    segments = list(c("1","2","3"), c("3","4","5"))
  )

  expect_equal(2 * 2, 4)
})
