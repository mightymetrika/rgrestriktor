# Decisions for informative hypothesis testing (B then A), with optional STRICT readout (C then B).
# Default behaviour is unchanged: we report B→A results plus IUT across rg_iht* segments,
# and ANOVA omnibus power as a comparator. If include_strict = TRUE, we append parallel
# "_C&B" rows that use C (strict) gated by B (compatibility).

#' Per-iteration decisions for IHT objects (B then A)
#'
#' @param iht_df Data frame produced by your extractor with rows for "iht" and "rg_iht*"
#'               and columns A_pval, B_pval, C_pval, global_pval, modname.
#' @param alpha  Significance level (default 0.05).
#' @return Tibble with columns: modname, support (B→A), A_p, B_p, C_p, global_p.
#'         `support` is TRUE iff (B_p >= alpha) & (A_p < alpha).
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
        global_p  = numeric(0)
      ))
    }
    A_p <- df$A_pval
    B_p <- df$B_pval
    C_p <- if ("C_pval" %in% names(df)) df$C_pval else rep(NA_real_, length(A_p))
    tibble::tibble(
      modname  = df$modname,
      # Canonical "evidence" rule: B not sig AND A sig
      support  = (B_p >= alpha) & (A_p < alpha),
      A_p      = A_p,
      B_p      = B_p,
      C_p      = C_p,
      global_p = df$global_pval
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

#' Summaries across simulation iterations
#'
#' @param sim_list       Output of rgr_iht_sim(): list of per-iteration results with $iht_df
#' @param alpha          Significance level (default 0.05).
#' @param include_strict If TRUE, also append rows using STRICT rule (C & B):
#'                       support_strict = (B_p >= alpha) & (C_p < alpha).
#'                       Default FALSE -> behaviour unchanged.
#' @return Tibble with rows:
#'   - "iht", "rg_iht*": power for B→A
#'   - "anova": omnibus power (comparator)
#'   - "rg_all_IUT": IUT across all segments (B→A)
#'   - "rg_any_unadj": at least one segment passes (exploratory)
#'   - (optional) corresponding "_C&B" rows when include_strict = TRUE
#' @export
rgr_summarize <- function(sim_list, alpha = 0.05, include_strict = FALSE) {
  # Gather per-iteration decisions + ANOVA flag
  decs <- purrr::imap_dfr(sim_list, function(one, i) {
    d <- rgr_decisions(one$iht_df, alpha = alpha)
    anova_p <- dplyr::filter(one$iht_df, .data$modname == "anova")$p[1]
    dplyr::mutate(d, iter = i, anova_sig = !is.na(anova_p) && (anova_p < alpha))
  })

  # ---- Primary (your original) B→A power ----
  pow_models <- decs |>
    dplyr::group_by(.data$modname) |>
    dplyr::summarise(power = mean(.data$support, na.rm = TRUE), .groups = "drop")

  pow_anova <- decs |>
    dplyr::distinct(.data$iter, .data$anova_sig) |>
    dplyr::summarise(power = mean(.data$anova_sig, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(modname = "anova") |>
    dplyr::select(.data$modname, .data$power)

  seg_only <- decs |>
    dplyr::filter(grepl("^rg_iht", .data$modname))

  add_iut_block <- function(df, tag = "") {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    by_iter <- df |>
      dplyr::group_by(.data$iter) |>
      dplyr::summarise(
        all_IUT   = all(.data$support, na.rm = TRUE),
        any_unadj = any(.data$support, na.rm = TRUE),
        .groups = "drop"
      )
    tibble::tibble(
      modname = c(paste0("rg_all_IUT", tag), paste0("rg_any_unadj", tag)),
      power   = c(mean(by_iter$all_IUT, na.rm = TRUE),
                  mean(by_iter$any_unadj, na.rm = TRUE))
    )
  }

  out <- dplyr::bind_rows(pow_models, pow_anova, add_iut_block(seg_only, tag = ""))

  # ---- Optional STRICT readout: C & B (OFF by default) ----
  if (isTRUE(include_strict)) {
    decs_strict <- decs |>
      dplyr::mutate(support = (B_p >= alpha) & (C_p < alpha),
                    modname  = paste0(.data$modname, "_C&B"))

    pow_models_strict <- decs_strict |>
      dplyr::group_by(.data$modname) |>
      dplyr::summarise(power = mean(.data$support, na.rm = TRUE), .groups = "drop")

    seg_only_strict <- decs_strict |>
      dplyr::filter(grepl("^rg_iht", .data$modname))

    out <- dplyr::bind_rows(
      out,
      pow_models_strict,
      add_iut_block(seg_only_strict, tag = "_C&B")
    )
  }

  out
}
