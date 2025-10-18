rgr_iht <- function(model = NULL, data = NULL, rgcontraints = NULL, constraints = NULL){

  # Set up data storage lists
  iht_df <- vector("list", sum(1, length(rgcontraints))) # tidy iht results
  rg_iht <- vector(mode = "list", length = length(rgcontraints)) # rg_iht models

  # Run iht with full constraints
  iht <- restriktor::iht(model, constraints = constraints)
  iht_df[[1]] <- extract_iht(iht, modname = "iht", constr = constraints)

  # Run iht with reduced graph constraints
  for (i in 1:length(rgcontraints)){
    rg_iht[[i]] <- restriktor::iht(model, constraints = rgcontraints[[i]])
    iht_df[[1 + i]] <- extract_iht(rg_iht[[i]], modname = paste0("rg_iht", i), constr = rgcontraints[[i]])

  }

  # Combine tidy results inot one dataframe
  iht_df <- do.call(dplyr::bind_rows, iht_df)

  # Return results
  return(list(iht_df = iht_df, iht = iht, rg_iht = rg_iht))
}


extract_iht <- function(ihtres, modname = "", constr = ""){

  # Global
  global_Ts <- ihtres$global$Ts
  global_pval <- ihtres$global$pvalue[[1]]
  global_R2org <- ihtres$global$R2.org
  global_R2red <- ihtres$global$R2.reduced
  global_reuH0 <- ihtres$global$b.eqrestr
  global_reuHA <- ihtres$global$b.unrestr

  global_df <- data.frame(global_Ts = global_Ts,
                          global_pval = global_pval,
                          global_R2org = global_R2org,
                          global_R2red = global_R2red)

  names(global_reuH0) <- paste0("global_reuH0_", names(global_reuH0))
  global_df <- cbind(global_df, t(as.data.frame(global_reuH0)))
  names(global_reuHA) <- paste0("global_reuHA_", names(global_reuHA))
  global_df <- cbind(global_df, t(as.data.frame(global_reuHA)))

  # Type A
  A_Ts <- ihtres$A$Ts
  A_pval <- ihtres$A$pvalue[[1]]
  A_R2org <- ihtres$A$R2.org
  A_R2red <- ihtres$A$R2.reduced
  A_reuH0 <- ihtres$A$b.eqrestr
  A_reuHA <- ihtres$A$b.unrestr

  A_df <- data.frame(A_Ts = A_Ts, A_pval = A_pval,
                     A_R2org = A_R2org, A_R2red = A_R2red)

  names(A_reuH0) <- paste0("A_reuH0_", names(A_reuH0))
  A_df <- cbind(A_df, t(as.data.frame(A_reuH0)))
  names(A_reuHA) <- paste0("A_reuHA_", names(A_reuHA))
  A_df <- cbind(A_df, t(as.data.frame(A_reuHA)))

  # Type B
  B_Ts <- ihtres$B$Ts
  B_pval <- ihtres$B$pvalue[[1]]
  B_R2org <- ihtres$B$R2.org
  B_R2red <- ihtres$B$R2.reduced
  B_reuH0 <- ihtres$B$b.restr
  B_ure <- ihtres$B$b.unrestr

  B_df <- data.frame(B_Ts = B_Ts, B_pval = B_pval,
                     B_R2org = B_R2org, B_R2red = B_R2red)

  names(B_reuH0) <- paste0("B_reuH0_", names(B_reuH0))
  B_df <- cbind(B_df, t(as.data.frame(B_reuH0)))
  names(B_ure) <- paste0("B_ure_", names(B_ure))
  B_df <- cbind(B_df, t(as.data.frame(B_ure)))

  # Type C
  C_Ts <- ihtres$C$Ts
  C_pval <- ihtres$C$pvalue[[1]]
  C_R2org <- ihtres$C$R2.org
  C_R2red <- ihtres$C$R2.reduced
  C_ure <- ihtres$C$b.unrestr

  C_df <- data.frame(C_Ts = C_Ts, C_pval = C_pval,
                     C_R2org = C_R2org, C_R2red = C_R2red)

  names(C_ure) <- paste0("C_ure_", names(C_ure))
  C_df <- cbind(C_df, t(as.data.frame(C_ure)))

  # Get dataframe information
  df_info <- data.frame(modname = modname, constr = constr)

  # Combine dataframes
  ihtres_df <- cbind(df_info, global_df, A_df, B_df, C_df)

  return(ihtres_df)
}
