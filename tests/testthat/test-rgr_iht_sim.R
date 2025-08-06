test_that("rhr_iht_sim works", {

  iterations <- 2

  res <-  rgr_iht_sim(Its = iterations,
                    rgcontraints = list(' group1 < group2 < group3 ',
                                        ' group3 < group4 < group5 '),
                    constraints = ' group1 < group2 < group3 < group4 < group5 ',
                    dgp = function(x){
                      data.frame(
                        y = stats::rnorm(25),
                        # x = stats::rnorm(25, 5),
                        group = as.factor(rep(c("1", "2", "3", "4", "5"), each = 5)))
                      },
                    .formula = "y ~ -1 + group")
  expect_equal(length(res), iterations)
})

