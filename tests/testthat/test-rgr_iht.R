test_that("test rgr_iht works", {

  #########################
  ## Example modified from ?restriktor::iht ##
  #########################
  # data generation process
  n <- 10
  means <- c(1,2,1,3, 4)
  nm <- length(means)
  group <- as.factor(rep(1:nm, each = n))
  y <- rnorm(n * nm, rep(means, each = n))
  DATA2 <- data.frame(y, group)

  # fit unrestricted linear model
  fit2.lm <- lm(y ~ -1 + group, data = DATA2)
  coef(fit2.lm)

  ## increasing means with five groups
  myConstraints2 <- ' group1 < group2 < group3 < group4 < group5 '

  myrgConstraints2 <- list(' group1 < group2 < group3 < group4 < group5 ',
                           ' group3 < group4 < group5 ')

  # fit rgr_iht model
  res <- rgr_iht(fit2.lm,
                 rgcontraints = myrgConstraints2,
                 constraints = myConstraints2) |>
    suppressWarnings()
  expect_equal(length(res), 2)
})
