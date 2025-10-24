#' Run a set of simulation scenarios and return one tidy summary table
#'
#' @param scenarios list of scenarios; each is a list with:
#'   - name:        character label for the scenario
#'   - rgconstraints: list of constraint strings for segments (can be NULL)
#'   - constraints: single full-chain constraint string (can be NULL)
#'   - dgp:         function() -> data.frame with columns y and group
#'   - formula:     model formula string for IHT, e.g., "y ~ -1 + group"
#' @param Its iterations per scenario
#' @param alpha significance level for IHT & ANOVA comparators
#' @param include_strict include strict (C&B) rows in the summary
#' @param include_anova_gate include ANOVA adjacent-pairs gatekeeper rows
#' @param gate_args named list passed to rgr_anova_gate_summarize()
#' @return tibble with columns: scenario, family, modname, power
#' @export
rgr_run_scenarios <- function(
    scenarios,
    Its = 1000,
    alpha = 0.05,
    include_strict = TRUE,
    include_anova_gate = TRUE,
    gate_args = list(
      require_omnibus = TRUE,
      segments = NULL,
      adj_method = "holm",
      alternative = c("greater","two.sided")
    )
) {
  stopifnot(is.list(scenarios), length(scenarios) > 0)

  out <- purrr::map_dfr(scenarios, function(sc) {
    stopifnot(is.list(sc), !is.null(sc$name), !is.null(sc$dgp), !is.null(sc$formula))

    res <- rgr_iht_sim(
      Its = Its,
      rgconstraints = sc$rgconstraints %||% NULL,  # <- fix spelling here
      constraints  = sc$constraints  %||% NULL,
      dgp          = sc$dgp,
      .formula     = sc$formula
    )

    summ <- rgr_summarize(
      res,
      alpha = alpha,
      include_strict = include_strict,
      include_anova_gate = include_anova_gate,
      gate_args = gate_args
    )

    dplyr::mutate(summ, scenario = sc$name, .before = 1)
  })

  # order scenarios by appearance
  out$scenario <- factor(out$scenario, levels = unique(out$scenario))
  out
}

#' Quick plot of power, grouped by family; optional caption and scenario facets
#'
#' @param tbl tibble from rgr_summarize() or rgr_run_scenarios()
#' @param Its iterations used (if provided, 95% & 68% Wilson CIs are drawn)
#' @param drop_families families to drop from plotting
#' @param caption optional plot caption
#' @param facet_by one of "auto","family","scenario","both"
#'        - "auto": if >1 scenario present, facet by both; else by family.
#' @param scenario_labels optional named vector full->label for facet strips
#' @param order_scenarios optional vector with desired scenario facet order (use labels if provided)
#' @param order_families optional vector with desired family facet order
#' @export
rgr_plot_power <- function(tbl,
                           Its = NULL,
                           drop_families = c("RG-IHT strict (C&B)", "IHT strict (C&B)"),
                           caption = NULL,
                           facet_by = c("auto","family","scenario","both","none"),
                           scenario_labels = NULL,
                           order_scenarios = NULL,
                           order_families = NULL,
                           strip_right = TRUE,
                           separate = c("none","scenario","family","both"),
                           ...) {

  facet_by <- match.arg(facet_by)
  separate <- match.arg(separate)

  # -- internal helpers (no new deps) -----------------------------------------
  .reorder_within <- function(x, by, within, fun = mean, desc = TRUE, sep = "___") {
    # x: character/factor; by: numeric; within: data.frame of facet vars
    w <- interaction(within, drop = TRUE)
    x_chr <- as.character(x)
    x_w <- paste(x_chr, w, sep = sep)
    o <- if (desc) stats::reorder(x_w, -by, FUN = fun) else stats::reorder(x_w, by, FUN = fun)
    o
  }
  .scale_x_reordered <- function(..., sep = "___") {
    ggplot2::scale_x_discrete(labels = function(x) sub(paste0(sep, ".*$"), "", x), ...)
  }

  # --- "separate" mode: split first, then plot with *sensible* facets ----------
  if (separate != "none") {
    df <- dplyr::as_tibble(tbl)

    split_vars <- switch(
      separate,
      scenario = "scenario",
      family   = "family",
      both     = c("scenario","family")
    )
    split_vars <- split_vars[split_vars %in% names(df)]
    stopifnot(length(split_vars) > 0L)

    # preserve user’s facet intent as much as possible
    local_facet <- facet_by
    if (separate == "scenario" && facet_by %in% c("auto","both")) local_facet <- "family"
    if (separate == "family"   && facet_by %in% c("auto","both")) local_facet <- "scenario"
    if (separate == "both") local_facet <- "none"

    grp_keys <- unique(df[split_vars])
    out <- vector("list", nrow(grp_keys))
    nm  <- apply(grp_keys, 1, function(r) paste0(names(r), "=", unname(unlist(r)), collapse = " | "))

    for (i in seq_len(nrow(grp_keys))) {
      idx <- rep(TRUE, nrow(df))
      for (v in split_vars) idx <- idx & df[[v]] == grp_keys[[v]][i]
      out[[i]] <- rgr_plot_power(
        df[idx, , drop = FALSE],
        Its = Its,
        caption = caption,
        facet_by = local_facet,
        scenario_labels = scenario_labels,
        order_scenarios = order_scenarios,
        order_families  = order_families,
        strip_right = strip_right,
        separate = "none",
        ...
      )
    }
    names(out) <- nm
    class(out) <- c("rgr_plot_list","list")
    return(out)
  }

  # -------------------- standard path (no splitting) --------------------------
  df <- dplyr::as_tibble(tbl)
  stopifnot("family" %in% names(df), "scenario" %in% names(df), "modname" %in% names(df), "power" %in% names(df))

  if (!is.null(drop_families)) df <- dplyr::filter(df, !.data$family %in% drop_families)

  # scenario labels + ordering (preserve original order unless user overrides)
  orig_levels <- if (is.factor(df$scenario)) levels(df$scenario) else unique(df$scenario)
  df$scenario_display <- as.character(df$scenario)
  if (!is.null(scenario_labels)) {
    lab <- unname(scenario_labels[df$scenario_display])
    df$scenario_display <- ifelse(is.na(lab), df$scenario_display, lab)
  }
  default_ord <- df |>
    dplyr::distinct(scenario, scenario_display) |>
    dplyr::arrange(factor(.data$scenario, levels = orig_levels)) |>
    dplyr::pull(scenario_display)
  scen_levels <- if (is.null(order_scenarios)) default_ord else order_scenarios
  df$scenario_display <- factor(df$scenario_display, levels = scen_levels)

  # family ordering (optional)
  if (!is.null(order_families)) df$family <- factor(df$family, levels = order_families)

  # 95%/68% Wilson CIs (same as before)
  if (!is.null(Its)) {
    wilson <- function(p, n, z) {
      denom  <- 1 + z^2 / n
      center <- (p + z^2/(2*n)) / denom
      half   <- z * sqrt((p*(1-p) + z^2/(4*n)) / n) / denom
      list(l = pmax(0, center - half), u = pmin(1, center + half))
    }
    ci95 <- wilson(df$power, Its, z = 1.96)
    ci68 <- wilson(df$power, Its, z = 1.00)
    df$lo95 <- ci95$l; df$hi95 <- ci95$u
    df$lo68 <- ci68$l; df$hi68 <- ci68$u
  }

  # Decide what "within" means for reordering the bars
  within_vars <- switch(
    facet_by,
    family   = c("family"),
    scenario = c("scenario_display"),
    both     = c("family","scenario_display"),
    auto     = if (length(levels(df$scenario_display)) > 1) c("family","scenario_display") else c("family"),
    none     = character(0)
  )
  # compute facet-aware x with descending power
  if (length(within_vars)) {
    df$modname_ord <- .reorder_within(df$modname, df$power, df[within_vars], fun = mean, desc = TRUE)
  } else {
    # no facets: one global order
    df$modname_ord <- stats::reorder(df$modname, -df$power, FUN = mean)
  }

  base <- ggplot2::ggplot(df, ggplot2::aes(x = modname_ord, y = power, fill = family)) +
    ggplot2::geom_col(width = 0.72) +
    { if (!is.null(Its)) ggplot2::geom_errorbar(ggplot2::aes(ymin = lo95, ymax = hi95), width = 0, alpha = 0.45) } +
    { if (!is.null(Its)) ggplot2::geom_errorbar(ggplot2::aes(ymin = lo68, ymax = hi68), width = 0, alpha = 0.8, linewidth = 1) } +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100*x), "%"), limits = c(0,1)) +
    .scale_x_reordered() +
    ggplot2::labs(x = NULL, y = "Power", fill = NULL, caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor   = ggplot2::element_blank(),
                   legend.position    = "none",
                   strip.text         = ggplot2::element_text(face = "bold"))

  # choose faceting
  if (facet_by == "auto") {
    facet_by <- if (length(levels(df$scenario_display)) > 1) "both" else "family"
  }
  p <- switch(
    facet_by,
    family   = base + ggplot2::facet_wrap(~ family, ncol = 2, scales = "free_y"),
    scenario = base + ggplot2::facet_wrap(~ scenario_display, ncol = 2, scales = "free_y"),
    both     = base + ggplot2::facet_grid(rows = ggplot2::vars(family),
                                          cols  = ggplot2::vars(scenario_display),
                                          scales = "free_y"),
    none     = base
  )

  if (!strip_right) {
    p <- p + ggplot2::theme(
      strip.text.y       = ggplot2::element_blank(),
      strip.text.y.right = ggplot2::element_blank(),
      strip.text.y.left  = ggplot2::element_blank(),
      strip.background.y = ggplot2::element_blank()
    )
  }

  p
}
# rgr_plot_power <- function(tbl,
#                            Its = NULL,
#                            drop_families = c("RG-IHT strict (C&B)", "IHT strict (C&B)"),
#                            caption = NULL,
#                            facet_by = c("auto","family","scenario","both","none"),
#                            scenario_labels = NULL,
#                            order_scenarios = NULL,
#                            order_families = NULL,
#                            strip_right = TRUE,
#                            separate = c("none","scenario","family","both"),
#                            ...) {
#
#   facet_by  <- match.arg(facet_by)
#   separate  <- match.arg(separate)
#
#   # --- NEW: "separate" mode (returns a named list of ggplots) ---
#   if (separate != "none") {
#     split_vars <- switch(
#       separate,
#       scenario = "scenario",
#       family   = "family",
#       both     = c("scenario","family")
#     )
#     df <- dplyr::as_tibble(tbl)            # ensure df exists here too
#     split_vars <- split_vars[split_vars %in% names(df)]
#     stopifnot(length(split_vars) > 0L)
#
#     grp_keys <- unique(df[split_vars])
#     out <- vector("list", nrow(grp_keys))
#     nm  <- apply(grp_keys, 1, function(r)
#       paste0(names(r), "=", unname(unlist(r)), collapse = " | ")
#     )
#
#     for (i in seq_len(nrow(grp_keys))) {
#       idx <- rep(TRUE, nrow(df))
#       for (v in split_vars) idx <- idx & df[[v]] == grp_keys[[v]][i]
#       out[[i]] <- rgr_plot_power(
#         df[idx, , drop = FALSE],
#         Its = Its,
#         caption = caption,
#         facet_by = "none",
#         scenario_labels = scenario_labels,
#         order_scenarios = order_scenarios,
#         strip_right = strip_right,
#         separate = "none",
#         ...
#       )
#     }
#     names(out) <- nm
#     class(out) <- c("rgr_plot_list","list")
#     return(out)
#   }
#
#
#
#   df <- dplyr::as_tibble(tbl)
#   stopifnot("family" %in% names(df), "scenario" %in% names(df))
#
#   if (!is.null(drop_families)) df <- dplyr::filter(df, !.data$family %in% drop_families)
#
#   # scenario display + ordering
#   orig_levels <- if (is.factor(df$scenario)) levels(df$scenario) else unique(df$scenario)
#   df$scenario_display <- as.character(df$scenario)
#   if (!is.null(scenario_labels)) {
#     lab <- unname(scenario_labels[df$scenario_display])
#     df$scenario_display <- ifelse(is.na(lab), df$scenario_display, lab)
#   }
#   default_ord <- df |>
#     dplyr::distinct(scenario, scenario_display) |>
#     dplyr::arrange(factor(.data$scenario, levels = orig_levels)) |>
#     dplyr::pull(scenario_display)
#   scen_levels <- if (is.null(order_scenarios)) default_ord else order_scenarios
#   df$scenario_display <- factor(df$scenario_display, levels = scen_levels)
#
#   # family ordering
#   if (!is.null(order_families)) df$family <- factor(df$family, levels = order_families)
#
#   # CIs
#   wilson <- function(p, n, z = 1.96) {
#     if (is.null(n)) return(list(l = NA_real_, u = NA_real_))
#     denom <- 1 + z^2 / n
#     center <- (p + z^2/(2*n)) / denom
#     half   <- z * sqrt((p*(1-p) + z^2/(4*n)) / n) / denom
#     list(l = pmax(0, center - half), u = pmin(1, center + half))
#   }
#   if (!is.null(Its)) {
#     ci95 <- wilson(df$power, Its, z = 1.96)
#     ci68 <- wilson(df$power, Its, z = 1.00)
#     df$lo95 <- ci95$l; df$hi95 <- ci95$u
#     df$lo68 <- ci68$l; df$hi68 <- ci68$u
#   }
#
#   df <- df |>
#     dplyr::group_by(.data$family) |>
#     dplyr::mutate(modname_ord = reorder(.data$modname, .data$power, FUN = mean)) |>
#     dplyr::ungroup()
#
#   # choose faceting
#   if (facet_by == "auto") {
#     facet_by <- if (length(levels(df$scenario_display)) > 1) "both" else "family"
#   }
#
#   base <- ggplot2::ggplot(df, ggplot2::aes(x = modname_ord, y = power, fill = family)) +
#     ggplot2::geom_col(width = 0.72) +
#     { if (!is.null(Its)) ggplot2::geom_errorbar(ggplot2::aes(ymin = lo95, ymax = hi95), width = 0, alpha = 0.45) } +
#     { if (!is.null(Its)) ggplot2::geom_errorbar(ggplot2::aes(ymin = lo68, ymax = hi68), width = 0, alpha = 0.8, linewidth = 1) } +
#     ggplot2::coord_flip() +
#     ggplot2::scale_y_continuous(labels = function(x) paste0(round(100*x), "%"), limits = c(0,1)) +
#     ggplot2::labs(x = NULL, y = "Power", fill = NULL, caption = caption) +
#     ggplot2::theme_minimal(base_size = 12) +
#     ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
#                    panel.grid.minor = ggplot2::element_blank(),
#                    legend.position = "none",
#                    strip.text = ggplot2::element_text(face = "bold"))
#
#   p <- switch(facet_by,
#               family   = base + ggplot2::facet_wrap(~ family, ncol = 2, scales = "free_y"),
#               scenario = base + ggplot2::facet_wrap(~ scenario_display, ncol = 2, scales = "free_y"),
#               both     = base + ggplot2::facet_grid(rows = ggplot2::vars(family),
#                                                     cols  = ggplot2::vars(scenario_display),
#                                                     scales = "free_y"),
#               none     = base
#   )
#
#   if (!strip_right) {
#     p <- p + ggplot2::theme(
#       strip.text.y        = ggplot2::element_blank(),
#       strip.text.y.right  = ggplot2::element_blank(),
#       strip.text.y.left   = ggplot2::element_blank(),
#       strip.background.y  = ggplot2::element_blank()
#     )
#   }
#
#   return(p)
# }


#' Pretty gt table for power (grouped by *scenario* or *family*)
#'
#' @param tbl tibble from rgr_summarize() or rgr_run_scenarios()
#' @param digits number of digits for power
#' @param group_by which column to use as row-group header: "scenario" (default) or "family"
#' @param scenario_labels optional named character vector mapping
#'   full scenario names -> display labels (e.g., c("Middle inversion delta=1, n=5"="MID")).
#'   Only affects display; the original factor order is preserved unless you override it below.
#' @param order_scenarios optional character vector giving the desired order of
#'   row-groups when group_by == "scenario" (use display labels if you passed scenario_labels).
#'   If NULL, the function preserves the original factor order from `rgr_run_scenarios()`.
#' @param order_families optional character vector giving the desired order of families
#'   when group_by == "family". If NULL, whatever is in `tbl$family` is kept.
#' @param notes optional character string appended as a source note (free-form).
#' @param show_scenario logical; if FALSE, hide the Scenario column in the body.
#' @return gt table object
#' @export
rgr_gt_power <- function(tbl,
                         digits = 3,
                         group_by = c("scenario", "family"),
                         scenario_labels = NULL,
                         order_scenarios = NULL,
                         order_families  = NULL,
                         notes = NULL,
                         show_scenario = TRUE,
                         separate = c("none","scenario","family"),
                         ...) {
  requireNamespace("gt", quietly = TRUE)
  group_by <- match.arg(group_by)

  df <- dplyr::as_tibble(tbl)
  stopifnot(all(c("scenario","family","modname","power") %in% names(df)))

  # Preserve original scenario order from rgr_run_scenarios()
  # (it's a factor in that function)
  orig_levels <- if (is.factor(df$scenario)) levels(df$scenario) else unique(df$scenario)

  # Build a display column for scenario WITHOUT destroying the original ordering
  df$scenario_display <- as.character(df$scenario)
  if (!is.null(scenario_labels)) {
    # map full -> label; keep unmapped as-is
    mapped <- unname(scenario_labels[df$scenario_display])
    df$scenario_display <- ifelse(is.na(mapped), df$scenario_display, mapped)
  }

  # Decide group orders
  if (group_by == "scenario") {
    # default order follows original scenario appearance (after mapping)
    default_ord <- df |>
      dplyr::distinct(scenario, scenario_display) |>
      dplyr::arrange(factor(.data$scenario, levels = orig_levels)) |>
      dplyr::pull(scenario_display)
    scen_levels <- if (is.null(order_scenarios)) default_ord else order_scenarios
    df$scenario_display <- factor(df$scenario_display, levels = scen_levels)
    ord <- df |>
      dplyr::arrange(.data$scenario_display, .data$family, dplyr::desc(.data$power))
  } else {
    if (!is.null(order_families)) {
      df$family <- factor(df$family, levels = order_families)
    }
    ord <- df |>
      dplyr::arrange(.data$family, factor(.data$scenario, levels = orig_levels), dplyr::desc(.data$power))
  }

  separate <- match.arg(separate)
  if (separate != "none") {
    split_var <- separate
    stopifnot(split_var %in% names(df))
    sp  <- split(df, df[[split_var]], drop = TRUE)
    out <- lapply(sp, function(dfi) {
      rgr_gt_power(
        dfi,
        group_by = group_by,
        scenario_labels = scenario_labels,
        order_scenarios = order_scenarios,
        order_families  = order_families,
        notes = notes,
        show_scenario = show_scenario,
        separate = "none",
        ...
      )
    })
    class(out) <- c("rgr_gt_list","list")
    return(out)
  }



  # Build table. Use groupname_col = chosen header; keep both Scenario/Family columns visible by default.
  gt_tbl <- gt::gt(
    ord,
    groupname_col = if (group_by == "scenario") "scenario_display" else "family"
  ) |>
    gt::fmt_number(columns = "power", decimals = digits) |>
    gt::tab_options(
      table.border.top.width = gt::px(0),
      table.border.bottom.width = gt::px(0),
      row.striping.include_table_body = TRUE
    ) |>
    gt::cols_label(
      scenario          = "Scenario",
      scenario_display  = "Scenario",
      family            = "Family",
      modname           = "Model",
      power             = "Power"
    )

  # Hide one of the scenario columns if desired
  if (group_by == "scenario") {
    gt_tbl <- gt_tbl |> gt::cols_hide(c("scenario"))  # keep the display one
    if (!show_scenario) gt_tbl <- gt_tbl |> gt::cols_hide("scenario_display")
  } else {
    if (!show_scenario) gt_tbl <- gt_tbl |> gt::cols_hide("scenario")
    gt_tbl <- gt_tbl |> gt::cols_hide("scenario_display")
  }

  # Free-form note from user (no automation)
  if (!is.null(notes)) {
    gt_tbl <- gt_tbl |> gt::tab_source_note(notes)
  }

  gt_tbl
}
