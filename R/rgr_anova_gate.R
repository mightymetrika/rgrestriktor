##' Adjacent-Pairs Gatekeeper summaries for ANOVA post-hoc
##'
##' Tests only the (k-1) adjacent differences along a prespecified order,
##' mirroring the minimal set of constraints behind a monotone order.
##' A run "supports" the ordered claim if EVERY adjacent pair is in the
##' correct direction (mean2>mean1) AND has adjusted p < alpha.
##' Optionally require the omnibus F-test to be significant as a gate.
##'
##' You can also evaluate restricted-graph segments (e.g., c("1","2","3"))
##' and compute an IUT across those segments.
##'
##' @param sim_list  Output of rgr_iht_sim(): list of per-iteration lists with $iht_df
##' @param alpha     Significance level (default 0.05)
##' @param require_omnibus  If TRUE, require ANOVA p<alpha as a gate (default TRUE)
##' @param order     Optional character vector giving the intended order of levels.
##'                  If NULL, will infer a numeric order from labels if possible,
##'                  otherwise sort labels alphabetically (per iteration).
##' @param segments  Optional list of character vectors, e.g.,
##'                  list(c("1","2","3"), c("3","4","5")) for restricted graphs.
##'                  If provided, segment-level AP-gate power and IUT are reported.
##' @param use_adjusted If TRUE (default), use p.adj from pairwise_t_test;
##'                     if FALSE, use raw p (not recommended).
##' @return tibble with power rows:
##'   - AP_gate_full            (full chain)
##'   - anova_AP_gate_full      (if require_omnibus=TRUE)
##'   - rg_anova_AP_gate{j}     (each segment, if segments provided)
##'   - AP_gate_IUT             (IUT across all segments, if segments provided)
##' @export
# rgr_anova_gate_summarize <- function(sim_list,
#                                      alpha = 0.05,
#                                      require_omnibus = TRUE,
#                                      order = NULL,
#                                      segments = NULL,
#                                      use_adjusted = TRUE) {
#
#   # Decide 'support' for one iteration
#   decide_iter <- function(iht_df) {
#     # Pull omnibus p
#     anova_p <- dplyr::filter(iht_df, .data$modname == "anova")$p[1]
#     pass_omni <- (!require_omnibus) || (!is.na(anova_p) && anova_p < alpha)
#
#     pairs <- dplyr::filter(iht_df, .data$modname == "pairwise_t_test")
#     need <- c("group1","group2","dir_ok","adj")
#     if (!all(need %in% names(pairs))) {
#       # If enrichment wasn't added, fail safely (no support)
#       return(dplyr::tibble(modname = character(0), support = logical(0)))
#     }
#
#     # Choose p-value column
#     pcol <- if (use_adjusted && "p.adj" %in% names(pairs)) "p.adj" else "p"
#
#     # Order of levels for "full chain"
#     if (is.null(order)) {
#       labs <- sort(unique(c(as.character(pairs$group1), as.character(pairs$group2))))
#       if (all(suppressWarnings(!is.na(as.numeric(labs))))) {
#         ord <- labs[order(as.numeric(labs))]
#       } else {
#         ord <- labs
#       }
#     } else {
#       ord <- as.character(order)
#     }
#
#     # Adjacent pairs we expect in the full chain
#     full_adj <- tibble::tibble(
#       group1 = head(ord, -1),
#       group2 = tail(ord, -1)
#     )
#
#     # Join to get only required adjacent comparisons
#     full_rows <- full_adj |>
#       dplyr::left_join(pairs, by = c("group1","group2"))
#
#     full_ok <- nrow(full_rows) == (length(ord) - 1) &&
#       all(full_rows$dir_ok, na.rm = FALSE) &&
#       all(!is.na(full_rows[[pcol]]) & (full_rows[[pcol]] < alpha))
#
#     out <- dplyr::tibble(
#       modname = if (require_omnibus) "anova_AP_gate_full" else "AP_gate_full",
#       support = pass_omni && full_ok
#     )
#
#     # Segment-wise AP-gate (restricted graph), if requested
#     if (!is.null(segments) && length(segments) > 0) {
#       seg_tbl <- purrr::imap_dfr(segments, function(seg, j) {
#         seg <- as.character(seg)
#         if (length(seg) < 2) {
#           return(dplyr::tibble(modname = paste0("rg_anova_AP_gate", j), support = FALSE))
#         }
#         seg_adj <- tibble::tibble(group1 = head(seg, -1), group2 = tail(seg, -1))
#         seg_rows <- seg_adj |>
#           dplyr::left_join(pairs, by = c("group1","group2"))
#
#         ok <- nrow(seg_rows) == (length(seg) - 1) &&
#           all(seg_rows$dir_ok, na.rm = FALSE) &&
#           all(!is.na(seg_rows[[pcol]]) & (seg_rows[[pcol]] < alpha))
#
#         dplyr::tibble(modname = paste0("rg_anova_AP_gate", j), support = pass_omni && ok)
#       })
#
#       out <- dplyr::bind_rows(out, seg_tbl)
#
#       # IUT across segments
#       out <- dplyr::bind_rows(
#         out,
#         dplyr::tibble(
#           modname = if (require_omnibus) "AP_gate_IUT" else "AP_gate_IUT",
#           support = all(dplyr::filter(out, grepl("^rg_anova_AP_gate", .data$modname))$support, na.rm = FALSE)
#         )
#       )
#     }
#
#     out
#   }
#
#   # iterate → aggregate power
#   decs <- purrr::imap_dfr(sim_list, function(one, i) {
#     decide_iter(one$iht_df) |>
#       dplyr::mutate(iter = i)
#   })
#
#   decs |>
#     dplyr::group_by(.data$modname) |>
#     dplyr::summarise(power = mean(.data$support, na.rm = TRUE), .groups = "drop")
# }




#' Adjacent-Pairs Gatekeeper summaries for ANOVA post-hoc (order-directed)
#'
#' Compares ANOVA post-hoc to IHT by testing only the (k-1) adjacent
#' differences along a prespecified order. A run "supports" the ordered claim
#' if EVERY adjacent pair is (a) in the correct direction (mean2 > mean1) and
#' (b) significant after multiplicity correction applied *within the adjacent set*.
#' Optionally, require the omnibus F-test to be significant as a gate.
#'
#' We reconstruct one-sided p-values from stored two-sided pairwise p's using
#' the observed direction (dir_ok). For symmetric tests (t), p_one = p_two/2
#' when dir_ok==TRUE, and p_one = 1 - p_two/2 otherwise.
#'
#' You can also evaluate restricted-graph segments (e.g., c("1","2","3"))
#' and compute an IUT across those segments.
#'
#' @param sim_list  Output of rgr_iht_sim(): list of per-iteration lists with $iht_df
#' @param alpha     Significance level (default 0.05)
#' @param require_omnibus  If TRUE, require ANOVA p<alpha as a gate (default TRUE)
#' @param order     Optional character vector giving the intended order of levels.
#'                  If NULL, infer a numeric order when labels are numeric; else
#'                  use sorted labels (per iteration).
#' @param segments  Optional list of character vectors, e.g.,
#'                  list(c("1","2","3"), c("3","4","5")) for restricted graphs.
#' @param adj_method P-adjust method within adjacent sets: "holm" (default),
#'                   "bonferroni", "BH", or "none".
#' @param alternative  "greater" (default, order-directed) or "two.sided".
#' @return tibble with power rows:
#'   - AP_gate_full (or anova_AP_gate_full if require_omnibus=TRUE)
#'   - rg_anova_AP_gate{j}     (each segment, if segments provided)
#'   - AP_gate_IUT             (IUT across segments, if provided)
#' @export
# rgr_anova_gate_summarize <- function(sim_list,
#                                      alpha = 0.05,
#                                      require_omnibus = TRUE,
#                                      order = NULL,
#                                      segments = NULL,
#                                      adj_method = c("holm","bonferroni","BH","none"),
#                                      alternative = c("greater","two.sided")) {
#   adj_method  <- match.arg(adj_method)
#   alternative <- match.arg(alternative)
#
#   # Helper: infer full order from labels if not supplied
#   infer_order <- function(pairs) {
#     labs <- sort(unique(c(as.character(pairs$group1), as.character(pairs$group2))))
#     if (all(suppressWarnings(!is.na(as.numeric(labs))))) {
#       labs[order(as.numeric(labs))]
#     } else {
#       labs
#     }
#   }
#
#   # Build one-sided p-values from two-sided and direction
#   to_one_sided <- function(p_two, dir_ok) {
#     if (alternative == "two.sided") return(p_two)
#     # "greater" (mean2 > mean1)
#     ifelse(is.na(p_two) | is.na(dir_ok),
#            NA_real_,
#            ifelse(dir_ok, p_two/2, 1 - p_two/2))
#   }
#
#   # Adjust p-values within a given set
#   padj_within <- function(p) {
#     if (adj_method == "none") p else stats::p.adjust(p, method = adj_method)
#   }
#
#   # Decide support for one iteration
#   decide_iter <- function(iht_df) {
#     # Omnibus gate
#     anova_p <- dplyr::filter(iht_df, .data$modname == "anova")$p[1]
#     pass_omni <- (!require_omnibus) || (!is.na(anova_p) && anova_p < alpha)
#
#     pairs <- dplyr::filter(iht_df, .data$modname == "pairwise_t_test")
#     need  <- c("group1","group2","dir_ok","adj","p")
#     if (!all(need %in% names(pairs))) {
#       return(dplyr::tibble(modname = character(0), support = logical(0)))
#     }
#
#     # Use supplied order or infer per iteration
#     ord <- if (is.null(order)) infer_order(pairs) else as.character(order)
#
#     # FULL-CHAIN adjacent pairs we expect
#     full_adj <- tibble::tibble(
#       group1 = head(ord, -1),
#       group2 = tail(ord, -1)
#     )
#
#     # Join to actual pairwise rows (must exist)
#     full_rows <- full_adj |>
#       dplyr::left_join(pairs, by = c("group1","group2"))
#
#     # One-sided (or two-sided) raw p-values for required pairs
#     p_one <- to_one_sided(full_rows$p, full_rows$dir_ok)
#     # Adjust WITHIN the adjacent set
#     p_adj <- padj_within(p_one)
#
#     full_ok <- nrow(full_rows) == (length(ord) - 1) &&
#       all(full_rows$dir_ok, na.rm = FALSE) &&
#       all(!is.na(p_adj) & (p_adj < alpha))
#
#     out <- dplyr::tibble(
#       modname = if (require_omnibus) "anova_AP_gate_full" else "AP_gate_full",
#       support = pass_omni && full_ok
#     )
#
#     # Segment-wise (restricted graph) AP-gate
#     if (!is.null(segments) && length(segments) > 0) {
#       seg_tbl <- purrr::imap_dfr(segments, function(seg, j) {
#         seg <- as.character(seg)
#         if (length(seg) < 2) {
#           return(dplyr::tibble(modname = paste0("rg_anova_AP_gate", j), support = FALSE))
#         }
#         seg_adj <- tibble::tibble(group1 = head(seg, -1), group2 = tail(seg, -1))
#         seg_rows <- seg_adj |>
#           dplyr::left_join(pairs, by = c("group1","group2"))
#         p1 <- to_one_sided(seg_rows$p, seg_rows$dir_ok)
#         pA <- padj_within(p1)
#         ok <- nrow(seg_rows) == (length(seg) - 1) &&
#           all(seg_rows$dir_ok, na.rm = FALSE) &&
#           all(!is.na(pA) & (pA < alpha))
#         dplyr::tibble(modname = paste0("rg_anova_AP_gate", j), support = pass_omni && ok)
#       })
#
#       out <- dplyr::bind_rows(out, seg_tbl)
#
#       out <- dplyr::bind_rows(
#         out,
#         dplyr::tibble(
#           modname = "AP_gate_IUT",
#           support = all(dplyr::filter(out, grepl("^rg_anova_AP_gate", .data$modname))$support, na.rm = FALSE)
#         )
#       )
#     }
#
#     out
#   }
#
#   # Iterate and aggregate
#   decs <- purrr::imap_dfr(sim_list, function(one, i) {
#     decide_iter(one$iht_df) |>
#       dplyr::mutate(iter = i)
#   })
#
#   decs |>
#     dplyr::group_by(.data$modname) |>
#     dplyr::summarise(power = mean(.data$support, na.rm = TRUE), .groups = "drop")
# }

rgr_anova_gate_summarize <- function(sim_list,
                                     alpha = 0.05,
                                     require_omnibus = TRUE,
                                     order = NULL,
                                     segments = NULL,
                                     adj_method = c("holm","bonferroni","BH","none"),
                                     alternative = c("greater","two.sided")) {
  adj_method  <- match.arg(adj_method)
  alternative <- match.arg(alternative)

  infer_order <- function(pairs) {
    labs <- sort(unique(c(as.character(pairs$group1), as.character(pairs$group2))))
    if (all(suppressWarnings(!is.na(as.numeric(labs))))) labs[order(as.numeric(labs))] else labs
  }

  padj_within <- function(p) if (adj_method == "none") p else stats::p.adjust(p, method = adj_method)

  # core helper: one-sided p from t, df, aligned so that positive t favors mean2>mean1
  # falls back to halving two-sided if t/df/diff are missing
  to_one_sided <- function(p_two, dir_ok, tstat, df, diff) {
    if (alternative == "two.sided") return(p_two)
    if (all(!is.na(tstat)) && all(!is.na(df)) && all(!is.na(diff))) {
      sign_diff <- ifelse(diff > 0, 1, ifelse(diff < 0, -1, 0))
      t_aligned <- tstat * sign_diff   # >0 when sample supports mean2>mean1
      return(stats::pt(t_aligned, df = df, lower.tail = FALSE))
    } else {
      # symmetric t fallback
      return(ifelse(is.na(p_two) | is.na(dir_ok),
                    NA_real_,
                    ifelse(dir_ok, p_two/2, 1 - p_two/2)))
    }
  }

  decide_iter <- function(iht_df) {
    anova_p <- dplyr::filter(iht_df, .data$modname == "anova")$p[1]
    pass_omni <- (!require_omnibus) || (!is.na(anova_p) && anova_p < alpha)

    pairs <- dplyr::filter(iht_df, .data$modname == "pairwise_t_test")
    need  <- c("group1","group2","dir_ok","p")  # statistic/df/diff are optional but preferred
    if (!all(need %in% names(pairs))) {
      return(dplyr::tibble(modname = character(0), support = logical(0)))
    }

    ord <- if (is.null(order)) infer_order(pairs) else as.character(order)

    full_adj <- tibble::tibble(group1 = head(ord, -1), group2 = tail(ord, -1))
    full_rows <- full_adj |>
      dplyr::left_join(pairs, by = c("group1","group2"))

    p_one <- to_one_sided(
      p_two = full_rows$p,
      dir_ok = full_rows$dir_ok,
      tstat = full_rows$statistic %||% NA_real_,
      df    = full_rows$df        %||% NA_real_,
      diff  = full_rows$diff      %||% NA_real_
    )
    p_adj <- padj_within(p_one)

    full_ok <- nrow(full_rows) == (length(ord) - 1) &&
      all(full_rows$dir_ok, na.rm = FALSE) &&
      all(!is.na(p_adj) & (p_adj < alpha))

    out <- dplyr::tibble(
      modname = if (require_omnibus) "anova_AP_gate_full" else "AP_gate_full",
      support = pass_omni && full_ok
    )

    if (!is.null(segments) && length(segments) > 0) {
      seg_tbl <- purrr::imap_dfr(segments, function(seg, j) {
        seg <- as.character(seg)
        if (length(seg) < 2) {
          return(dplyr::tibble(modname = paste0("rg_anova_AP_gate", j), support = FALSE))
        }
        seg_adj <- tibble::tibble(group1 = head(seg, -1), group2 = tail(seg, -1))
        seg_rows <- seg_adj |>
          dplyr::left_join(pairs, by = c("group1","group2"))
        p1 <- to_one_sided(
          p_two = seg_rows$p,
          dir_ok = seg_rows$dir_ok,
          tstat = seg_rows$statistic %||% NA_real_,
          df    = seg_rows$df        %||% NA_real_,
          diff  = seg_rows$diff      %||% NA_real_
        )
        pA <- padj_within(p1)
        ok <- nrow(seg_rows) == (length(seg) - 1) &&
          all(seg_rows$dir_ok, na.rm = FALSE) &&
          all(!is.na(pA) & (pA < alpha))
        dplyr::tibble(modname = paste0("rg_anova_AP_gate", j), support = pass_omni && ok)
      })

      out <- dplyr::bind_rows(out, seg_tbl)

      out <- dplyr::bind_rows(
        out,
        dplyr::tibble(
          modname = "AP_gate_IUT",
          support = all(dplyr::filter(out, grepl("^rg_anova_AP_gate", .data$modname))$support, na.rm = FALSE)
        )
      )
    }

    out
  }

  decs <- purrr::imap_dfr(sim_list, function(one, i) {
    decide_iter(one$iht_df) |>
      dplyr::mutate(iter = i)
  })

  decs |>
    dplyr::group_by(.data$modname) |>
    dplyr::summarise(power = mean(.data$support, na.rm = TRUE), .groups = "drop")
}

