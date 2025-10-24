# Decisions for informative hypothesis testing (B then A), with optional STRICT (C then B)
# and optional ANOVA adjacent-pairs gatekeeper appended to the main summary.

#' Per-iteration decisions for IHT objects (B then A)
#' (unchanged from your current version, included here for completeness)
#' @export
rgr_decisions <- function(iht_df, alpha = 0.05) {
  make_decision_rows <- function(df) {
    if (nrow(df) == 0) {
      return(tibble::tibble(
        modname   = character(0),
        support   = logical(0),
        A_p       = numeric(0),
        B_p       = numeric(0),
        C_p       = numeric(0),
        global_p  = numeric(0),
        fit_ok    = logical(0)
      ))
    }
    A_p <- df$A_pval
    B_p <- df$B_pval
    C_p <- if ("C_pval" %in% names(df)) df$C_pval else rep(NA_real_, length(A_p))
    fit_ok_col <- if ("fit_ok" %in% names(df)) df$fit_ok else rep(NA, length(A_p))
    tibble::tibble(
      modname  = df$modname,
      support  = (B_p >= alpha) & (A_p < alpha),
      A_p      = A_p,
      B_p      = B_p,
      C_p      = C_p,
      global_p = df$global_pval,
      fit_ok   = fit_ok_col
    )
  }

  full_out <- iht_df |>
    dplyr::filter(.data$modname == "iht") |>
    dplyr::slice(1) |>
    make_decision_rows()

  segs_out <- iht_df |>
    dplyr::filter(grepl("^rg_iht", .data$modname)) |>
    make_decision_rows()

  dplyr::bind_rows(full_out, segs_out)
}
# rgr_decisions <- function(iht_df, alpha = 0.05) {
#   make_decision_rows <- function(df) {
#     if (nrow(df) == 0) {
#       return(tibble::tibble(
#         modname   = character(0),
#         support   = logical(0),
#         A_p       = numeric(0),
#         B_p       = numeric(0),
#         C_p       = numeric(0),
#         global_p  = numeric(0)
#       ))
#     }
#     A_p <- df$A_pval
#     B_p <- df$B_pval
#     C_p <- if ("C_pval" %in% names(df)) df$C_pval else rep(NA_real_, length(A_p))
#     tibble::tibble(
#       modname  = df$modname,
#       support  = (B_p >= alpha) & (A_p < alpha),
#       A_p      = A_p,
#       B_p      = B_p,
#       C_p      = C_p,
#       global_p = df$global_pval
#     )
#   }
#
#   full_out <- iht_df |>
#     dplyr::filter(.data$modname == "iht") |>
#     dplyr::slice(1) |>
#     make_decision_rows()
#
#   segs_out <- iht_df |>
#     dplyr::filter(grepl("^rg_iht", .data$modname)) |>
#     make_decision_rows()
#
#   dplyr::bind_rows(full_out, segs_out)
# }

#' Summaries across simulation iterations with optional STRICT and ANOVA gate rows.
#' Adds a per-model fit_rate column and returns families in a fixed, meaningful order.
#'
#' @param sim_list         Output of rgr_iht_sim(): list of per-iteration results with $iht_df
#' @param alpha            Significance level (default 0.05).
#' @param include_strict   If TRUE, append STRICT rows using (B_p >= alpha) & (C_p < alpha).
#' @param include_anova_gate If TRUE, append ANOVA AP-gate rows (see gate_args).
#' @param gate_args        Named list passed to rgr_anova_gate_summarize().
#' @return Tibble with columns: family, modname, power, fit_rate
#' @export
rgr_summarize <- function(sim_list,
                          alpha = 0.05,
                          include_strict = FALSE,
                          include_anova_gate = FALSE,
                          gate_args = list()) {

  # Per-iteration IHT decisions (+ fit_ok), plus ANOVA flag for omnibus
  decs <- purrr::imap_dfr(sim_list, function(one, i) {
    d <- rgr_decisions(one$iht_df, alpha = alpha)
    anova_p <- dplyr::filter(one$iht_df, .data$modname == "anova")$p[1]
    dplyr::mutate(d, iter = i, anova_sig = !is.na(anova_p) && (anova_p < alpha))
  })

  # Power and fit-rate for direct IHT and RG-IHT models
  pow_models <- decs |>
    dplyr::group_by(.data$modname) |>
    dplyr::summarise(
      power    = mean(.data$support, na.rm = TRUE),
      fit_rate = mean(.data$fit_ok, na.rm = TRUE),
      .groups  = "drop"
    )

  # ANOVA omnibus (no fit_rate; set NA)
  pow_anova <- decs |>
    dplyr::distinct(.data$iter, .data$anova_sig) |>
    dplyr::summarise(power = mean(.data$anova_sig, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(modname = "anova", fit_rate = NA_real_) |>
    dplyr::select(.data$modname, .data$power, .data$fit_rate)

  # IUT/any across RG segments; define fit_rate for these composites
  seg_only <- decs |>
    dplyr::filter(grepl("^rg_iht", .data$modname))

  add_iut_block <- function(df, tag = "") {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    by_iter <- df |>
      dplyr::group_by(.data$iter) |>
      dplyr::summarise(
        all_IUT     = all(.data$support, na.rm = TRUE),
        any_unadj   = any(.data$support, na.rm = TRUE),
        all_fit_ok  = all(.data$fit_ok == TRUE, na.rm = FALSE),
        any_fit_ok  = any(.data$fit_ok == TRUE, na.rm = FALSE),
        .groups = "drop"
      )
    tibble::tibble(
      modname  = c(paste0("rg_all_IUT", tag), paste0("rg_any_unadj", tag)),
      power    = c(mean(by_iter$all_IUT, na.rm = TRUE),
                   mean(by_iter$any_unadj, na.rm = TRUE)),
      fit_rate = c(mean(by_iter$all_fit_ok, na.rm = TRUE),
                   mean(by_iter$any_fit_ok, na.rm = TRUE))
    )
  }

  out_BA <- dplyr::bind_rows(pow_models, pow_anova, add_iut_block(seg_only, tag = ""))

  # STRICT readout: C & B (carry fit_ok through unchanged)
  out <- out_BA
  if (isTRUE(include_strict)) {
    decs_strict <- decs |>
      dplyr::mutate(support = (.data$B_p >= alpha) & (.data$C_p < alpha),
                    modname  = paste0(.data$modname, "_C&B"))

    pow_models_strict <- decs_strict |>
      dplyr::group_by(.data$modname) |>
      dplyr::summarise(
        power    = mean(.data$support, na.rm = TRUE),
        fit_rate = mean(.data$fit_ok, na.rm = TRUE),
        .groups  = "drop"
      )

    seg_only_strict <- decs_strict |>
      dplyr::filter(grepl("^rg_iht", .data$modname))

    out <- dplyr::bind_rows(
      out,
      pow_models_strict,
      add_iut_block(seg_only_strict, tag = "_C&B")
    )
  }

  # Optional ANOVA AP-gate comparator (no fit_rate; set NA)
  if (isTRUE(include_anova_gate)) {
    gate_tbl <- do.call(
      rgr_anova_gate_summarize,
      c(list(sim_list = sim_list, alpha = alpha), gate_args)
    ) |>
      dplyr::mutate(fit_rate = NA_real_)
    out <- dplyr::bind_rows(out, gate_tbl)
  }

  # ---- Family labeling and stable family order ----
  # Recognize both suffixed (_1s/_2s) and unsuffixed AP-gate names.
  add_family <- function(name) {
    if (identical(name, "iht")) return("IHT (B->A)")
    if (grepl("^rg_iht", name) || name %in% c("rg_all_IUT", "rg_any_unadj")) return("RG-IHT (B->A)")
    if (identical(name, "iht_C&B")) return("IHT strict (C&B)")
    if (grepl("^rg_.*_C&B$", name) || name %in% c("rg_all_IUT_C&B", "rg_any_unadj_C&B")) return("RG-IHT strict (C&B)")
    if (identical(name, "anova")) return("ANOVA omnibus")
    # AP-gate families
    if (grepl("(_1s)$", name)) return("ANOVA AP-gate (1-sided)")
    if (grepl("(_2s)$", name)) return("ANOVA AP-gate (2-sided)")
    if (grepl("^(rg_anova_AP_gate[0-9]|AP_gate_IUT|anova_AP_gate_full)$", name)) return("ANOVA AP-gate (1-sided)")
    return("Other")
  }

  family_levels <- c(
    "IHT (B->A)",
    "RG-IHT (B->A)",
    "IHT strict (C&B)",
    "RG-IHT strict (C&B)",
    "ANOVA omnibus",
    "ANOVA AP-gate (1-sided)",
    "ANOVA AP-gate (2-sided)"
  )

  out <- out |>
    dplyr::mutate(
      family = vapply(.data$modname, add_family, character(1L)),
      family = factor(.data$family, levels = family_levels)
    ) |>
    dplyr::select(.data$family, .data$modname, .data$power, .data$fit_rate)

  out
}


# rgr_summarize <- function(sim_list,
#                           alpha = 0.05,
#                           include_strict = FALSE,
#                           include_anova_gate = FALSE,
#                           gate_args = list()) {
#   # Gather per-iteration decisions + ANOVA flag
#   decs <- purrr::imap_dfr(sim_list, function(one, i) {
#     d <- rgr_decisions(one$iht_df, alpha = alpha)
#     anova_p <- dplyr::filter(one$iht_df, .data$modname == "anova")$p[1]
#     dplyr::mutate(d, iter = i, anova_sig = !is.na(anova_p) && (anova_p < alpha))
#   })
#
#   # Primary (B→A)
#   pow_models <- decs |>
#     dplyr::group_by(.data$modname) |>
#     dplyr::summarise(power = mean(.data$support, na.rm = TRUE), .groups = "drop")
#
#   pow_anova <- decs |>
#     dplyr::distinct(.data$iter, .data$anova_sig) |>
#     dplyr::summarise(power = mean(.data$anova_sig, na.rm = TRUE), .groups = "drop") |>
#     dplyr::mutate(modname = "anova") |>
#     dplyr::select(.data$modname, .data$power)
#
#   seg_only <- decs |>
#     dplyr::filter(grepl("^rg_iht", .data$modname))
#
#   add_iut_block <- function(df, tag = "") {
#     if (is.null(df) || nrow(df) == 0) return(NULL)
#     by_iter <- df |>
#       dplyr::group_by(.data$iter) |>
#       dplyr::summarise(
#         all_IUT   = all(.data$support, na.rm = TRUE),
#         any_unadj = any(.data$support, na.rm = TRUE),
#         .groups = "drop"
#       )
#     tibble::tibble(
#       modname = c(paste0("rg_all_IUT", tag), paste0("rg_any_unadj", tag)),
#       power   = c(mean(by_iter$all_IUT, na.rm = TRUE),
#                   mean(by_iter$any_unadj, na.rm = TRUE))
#     )
#   }
#
#   out <- dplyr::bind_rows(pow_models, pow_anova, add_iut_block(seg_only, tag = ""))
#
#   # Optional STRICT (C & B)
#   if (isTRUE(include_strict)) {
#     decs_strict <- decs |>
#       dplyr::mutate(support = (B_p >= alpha) & (C_p < alpha),
#                     modname  = paste0(.data$modname, "_C&B"))
#
#     pow_models_strict <- decs_strict |>
#       dplyr::group_by(.data$modname) |>
#       dplyr::summarise(power = mean(.data$support, na.rm = TRUE), .groups = "drop")
#
#     seg_only_strict <- decs_strict |>
#       dplyr::filter(grepl("^rg_iht", .data$modname))
#
#     out <- dplyr::bind_rows(
#       out,
#       pow_models_strict,
#       add_iut_block(seg_only_strict, tag = "_C&B")
#     )
#   }
#
#   # # Optional ANOVA AP-gate comparator (adjacent pairs, order-directed)
#   # if (isTRUE(include_anova_gate)) {
#   #   gate_tbl <- do.call(
#   #     rgr_anova_gate_summarize,
#   #     c(list(sim_list = sim_list, alpha = alpha), gate_args)
#   #   )
#   #   out <- dplyr::bind_rows(out, gate_tbl)
#   # }
#   # Optional ANOVA AP-gate comparator (adjacent pairs) — allow both 1s & 2s in one call
#   if (isTRUE(include_anova_gate)) {
#     gate_defaults <- list(
#       require_omnibus = TRUE,
#       segments = NULL,
#       adj_method = "holm",
#       alternative = "greater"   # can be length 1 or c("greater","two.sided")
#     )
#     args <- utils::modifyList(gate_defaults, gate_args)
#
#     alts <- args$alternative
#     if (is.null(alts)) alts <- "greater"
#     alts <- as.character(alts)
#
#     for (alt in alts) {
#       args$alternative <- alt
#       gate_tbl <- do.call(
#         rgr_anova_gate_summarize,
#         c(list(sim_list = sim_list, alpha = alpha), args)
#       )
#       suff <- if (identical(alt, "greater")) "_1s" else "_2s"
#       gate_tbl$modname <- paste0(gate_tbl$modname, suff)
#       out <- dplyr::bind_rows(out, gate_tbl)
#     }
#   }
#
#   # --- Add a slim family label for presentation ---
#   # out <- out |>
#   #   dplyr::mutate(
#   #     family = dplyr::case_when(
#   #       .data$modname == "anova" ~ "ANOVA omnibus",
#   #       grepl("_C&B$", .data$modname) ~ "IHT strict (C&B)",
#   #       grepl("^anova_AP_gate_full", .data$modname) |
#   #         grepl("^rg_anova_AP_gate", .data$modname) |
#   #         grepl("^AP_gate_IUT", .data$modname) ~ NA_character_,  # fill below
#   #       TRUE ~ "IHT (B→A)"
#   #     )
#   #   ) |>
#   #   dplyr::mutate(
#   #     family = dplyr::if_else(
#   #       is.na(.data$family) & grepl("_1s$", .data$modname),
#   #       "ANOVA AP-gate (1-sided)",
#   #       dplyr::if_else(
#   #         is.na(.data$family) & grepl("_2s$", .data$modname),
#   #         "ANOVA AP-gate (2-sided)",
#   #         dplyr::if_else(is.na(.data$family), "ANOVA AP-gate", .data$family)
#   #       )
#   #     )
#   #   ) |>
#   #   dplyr::relocate(family, .before = .data$modname)
#   #
#   # # (optional) sort by family for nicer tables
#   # fam_levels <- c("IHT (B→A)", "IHT strict (C&B)", "ANOVA omnibus",
#   #                 "ANOVA AP-gate (1-sided)", "ANOVA AP-gate (2-sided)")
#   # out <- out |>
#   #   dplyr::mutate(family = factor(.data$family, levels = fam_levels)) |>
#   #   dplyr::arrange(.data$family, .data$modname)
#
#   out <- out |>
#     dplyr::mutate(
#       is_strict   = grepl("_C&B$", .data$modname),
#       is_rg       = grepl("^rg_iht", .data$modname) |
#         .data$modname %in% c("rg_all_IUT", "rg_any_unadj",
#                              "rg_all_IUT_C&B", "rg_any_unadj_C&B"),
#       is_iht_full = .data$modname %in% c("iht", "iht_C&B"),
#
#       # Family labels
#       family = dplyr::case_when(
#         .data$modname == "anova" ~ "ANOVA omnibus",
#
#         # ANOVA AP-gate (suffixes from our earlier patch)
#         grepl("_1s$", .data$modname) &
#           (grepl("^anova_AP_gate_full", .data$modname) |
#              grepl("^rg_anova_AP_gate", .data$modname) |
#              grepl("^AP_gate_IUT", .data$modname)) ~ "ANOVA AP-gate (1-sided)",
#         grepl("_2s$", .data$modname) &
#           (grepl("^anova_AP_gate_full", .data$modname) |
#              grepl("^rg_anova_AP_gate", .data$modname) |
#              grepl("^AP_gate_IUT", .data$modname)) ~ "ANOVA AP-gate (2-sided)",
#
#         # IHT vs RG-IHT, B→A vs strict C&B
#         .data$is_iht_full & !.data$is_strict ~ "IHT (B→A)",
#         .data$is_rg       & !.data$is_strict ~ "RG-IHT (B→A)",
#         .data$is_iht_full &  .data$is_strict ~ "IHT strict (C&B)",
#         .data$is_rg       &  .data$is_strict ~ "RG-IHT strict (C&B)",
#
#         TRUE ~ "Other"
#       ),
#
#       # Optional: a friendly sort order within family
#       mod_order = dplyr::case_when(
#         # IHT
#         .data$modname == "iht" ~ 1,
#         # RG-IHT B→A
#         .data$modname == "rg_all_IUT" ~ 1,
#         grepl("^rg_iht", .data$modname) & !.data$is_strict ~ 2,
#         .data$modname == "rg_any_unadj" ~ 3,
#         # IHT strict
#         .data$modname == "iht_C&B" ~ 1,
#         # RG-IHT strict
#         .data$modname == "rg_all_IUT_C&B" ~ 1,
#         grepl("^rg_iht", .data$modname) & .data$is_strict ~ 2,
#         .data$modname == "rg_any_unadj_C&B" ~ 3,
#         # ANOVA omnibus
#         .data$modname == "anova" ~ 1,
#         # ANOVA AP-gate (put full before segments, then IUT)
#         grepl("^anova_AP_gate_full", .data$modname) ~ 1,
#         grepl("^rg_anova_AP_gate", .data$modname) ~ 2,
#         grepl("^AP_gate_IUT", .data$modname) ~ 3,
#         TRUE ~ 99
#       )
#     ) |>
#     dplyr::relocate(.data$family, .before = .data$modname)
#
#   # Order families as: IHT → RG-IHT → IHT strict → RG-IHT strict → ANOVA omnibus → AP-gate (1s/2s)
#   fam_levels <- c("IHT (B→A)", "RG-IHT (B→A)",
#                   "IHT strict (C&B)", "RG-IHT strict (C&B)",
#                   "ANOVA omnibus",
#                   "ANOVA AP-gate (1-sided)", "ANOVA AP-gate (2-sided)")
#   out <- out |>
#     dplyr::mutate(family = factor(.data$family, levels = fam_levels)) |>
#     dplyr::arrange(.data$family, .data$mod_order, .data$modname) |>
#     dplyr::select(-.data$mod_order, -is_strict, -is_rg, -is_iht_full)
#
#   out
#
# }
