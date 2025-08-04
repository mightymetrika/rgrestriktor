test_that("test rgr_iht works", {

  #########################
  ## Example modified from ?restriktor::iht ##
  #########################
  # data generation process
  n <- 10
  means <- c(1,2,1,3,4)
  nm <- length(means)
  group <- as.factor(rep(1:nm, each = n))
  y <- rnorm(n * nm, rep(means, each = n))
  DATA2 <- data.frame(y, group)

  # fit unrestricted linear model
  fit2.lm <- lm(y ~ -1 + group, data = DATA2)
  coef(fit2.lm)

  ## Set constraints
  myConstraints2 <- ' group1 < group2 < group3 < group4 < group5 '
  myrgConstraints2 <- list(' group1 < group2 < group3 ',
                           ' group3 < group4 < group5 ')

  # fit rgr_iht model
  res <- rgr_iht(fit2.lm,
                 rgcontraints = myrgConstraints2,
                 constraints = myConstraints2) |>
    suppressWarnings()

  expect_equal(length(res), 3)
})

# res$iht$global$pvalue
#
# global_Ts <- res$iht$global$Ts
# global_pval <- res$iht$global$pvalue[[1]]
# global_R2org <- res$iht$global$R2.org
# global_R2red <- res$iht$global$R2.reduced
# global_reuH0 <- res$iht$global$b.eqrestr
# global_reuHA <- res$iht$global$b.unrestr
#
# global_df <- data.frame(global_Ts = global_Ts,
#                         global_pval = global_pval,
#                         global_R2org = global_R2red)
#
# names(global_reuH0) <- paste0("global_reuH0_", names(global_reuH0))
# global_df <- cbind(global_df, t(as.data.frame(global_reuH0)))
# names(global_reuHA) <- paste0("global_reuHA_", names(global_reuHA))
# global_df <- cbind(global_df, t(as.data.frame(global_reuHA)))
#
#
# for (i in length(1:length(global_reuH0))){
#   col_nm <- names(global_reuH0)[i]
#   global_df <- cbind(global_df, data.frame(assign("name", col_nm) = global_reuH0[i]))
# }
