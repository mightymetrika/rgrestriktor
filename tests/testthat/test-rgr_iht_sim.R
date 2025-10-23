test_that("rhr_iht_sim works", {

  iterations <- 2

  res <-  rgr_iht_sim(Its = iterations,
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
  expect_equal(length(res), iterations)
})

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

rgr_summarize(
  res1,
  include_strict = TRUE,
  include_anova_gate = TRUE,
  gate_args = list(
    require_omnibus = TRUE,
    segments = list(c("1","2","3"), c("3","4","5")),
    adj_method = "holm",
    alternative = c("greater","two.sided")  # <-- both variants appended
  )
)

res2 <-  rgr_iht_sim(Its = iterations,
                    rgcontraints = list(' group1 < group2 < group3 ',
                                        ' group3 < group4 < group5 '),
                    constraints = ' group1 < group2 < group3 < group4 < group5 ',
                    dgp = function(x){
                      data.frame(
                        y = c(stats::rnorm(5),
                              stats::rnorm(5, 3),
                              stats::rnorm(5, 2),
                              stats::rnorm(5, 1),
                              stats::rnorm(5, 4)),
                        # x = stats::rnorm(25, 5),
                        group = as.factor(rep(c("1", "2", "3", "4", "5"), each = 5)))
                    },
                    .formula = "y ~ -1 + group")

rgr_summarize(
  res2,
  include_strict = TRUE,
  include_anova_gate = TRUE,
  gate_args = list(
    require_omnibus = TRUE,
    segments = list(c("1","2","3"), c("3","4","5")),
    adj_method = "holm",
    alternative = c("greater","two.sided")  # <-- both variants appended
  )
)
