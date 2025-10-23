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
      rgcontraints = sc$rgconstraints %||% NULL,
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


##' Quick plot of power, grouped by family; optional caption for notes/abbrevs
##'
##' @param tbl  tibble from rgr_summarize() or rgr_run_scenarios()
##' @param Its  iterations used (if provided, 95% & 68% Wilson CIs are drawn)
##' @param drop_families families to drop from plotting
##' @param caption optional string placed as plot caption (e.g., scenario abbrev legend)
##' @return ggplot
##' @export
# rgr_plot_power <- function(tbl,
#                            Its = NULL,
#                            drop_families = c("RG-IHT strict (C&B)", "IHT strict (C&B)"),
#                            caption = NULL) {
#   requireNamespace("ggplot2", quietly = TRUE)
#
#   df <- dplyr::as_tibble(tbl)
#   if (!("family" %in% names(df))) {
#     stop("Expected a 'family' column; call rgr_summarize(...), which adds it.")
#   }
#
#   if (!is.null(drop_families)) {
#     df <- dplyr::filter(df, !.data$family %in% drop_families)
#   }
#
#   # Wilson CI helper
#   wilson <- function(p, n, z = 1.96) {
#     if (is.null(n)) return(list(l = NA_real_, u = NA_real_))
#     denom <- 1 + z^2 / n
#     center <- (p + z^2/(2*n)) / denom
#     half   <- z * sqrt((p*(1-p) + z^2/(4*n)) / n) / denom
#     list(l = pmax(0, center - half), u = pmin(1, center + half))
#   }
#
#   if (!is.null(Its)) {
#     ci95 <- wilson(df$power, Its, z = 1.96)
#     ci68 <- wilson(df$power, Its, z = 1.00)
#     df$lo95 <- ci95$l; df$hi95 <- ci95$u
#     df$lo68 <- ci68$l; df$hi68 <- ci68$u
#   }
#
#   # Reorder within family by power
#   df <- df |>
#     dplyr::group_by(.data$family) |>
#     dplyr::mutate(modname_ord = reorder(.data$modname, .data$power, FUN = mean)) |>
#     dplyr::ungroup()
#
#   p <- ggplot2::ggplot(df, ggplot2::aes(x = modname_ord, y = power, fill = family)) +
#     ggplot2::geom_col(width = 0.72) +
#     { if (!is.null(Its)) ggplot2::geom_errorbar(ggplot2::aes(ymin = lo95, ymax = hi95), width = 0, alpha = 0.45) } +
#     { if (!is.null(Its)) ggplot2::geom_errorbar(ggplot2::aes(ymin = lo68, ymax = hi68), width = 0, alpha = 0.8, linewidth = 1) } +
#     ggplot2::coord_flip() +
#     ggplot2::facet_wrap(~ family, ncol = 2, scales = "free_y") +
#     ggplot2::scale_y_continuous(labels = function(x) paste0(round(100*x), "%"), limits = c(0,1)) +
#     ggplot2::labs(x = NULL, y = "Power", fill = NULL, caption = caption) +
#     ggplot2::theme_minimal(base_size = 12) +
#     ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
#                    panel.grid.minor = ggplot2::element_blank(),
#                    legend.position = "none",
#                    strip.text = ggplot2::element_text(face = "bold"))
#
#   p
# }


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
                           facet_by = c("auto","family","scenario","both"),
                           scenario_labels = NULL,
                           order_scenarios = NULL,
                           order_families = NULL) {
  requireNamespace("ggplot2", quietly = TRUE)
  facet_by <- match.arg(facet_by)

  df <- dplyr::as_tibble(tbl)
  stopifnot("family" %in% names(df), "scenario" %in% names(df))

  if (!is.null(drop_families)) df <- dplyr::filter(df, !.data$family %in% drop_families)

  # scenario display + ordering
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

  # family ordering
  if (!is.null(order_families)) df$family <- factor(df$family, levels = order_families)

  # CIs
  wilson <- function(p, n, z = 1.96) {
    if (is.null(n)) return(list(l = NA_real_, u = NA_real_))
    denom <- 1 + z^2 / n
    center <- (p + z^2/(2*n)) / denom
    half   <- z * sqrt((p*(1-p) + z^2/(4*n)) / n) / denom
    list(l = pmax(0, center - half), u = pmin(1, center + half))
  }
  if (!is.null(Its)) {
    ci95 <- wilson(df$power, Its, z = 1.96)
    ci68 <- wilson(df$power, Its, z = 1.00)
    df$lo95 <- ci95$l; df$hi95 <- ci95$u
    df$lo68 <- ci68$l; df$hi68 <- ci68$u
  }

  df <- df |>
    dplyr::group_by(.data$family) |>
    dplyr::mutate(modname_ord = reorder(.data$modname, .data$power, FUN = mean)) |>
    dplyr::ungroup()

  # choose faceting
  if (facet_by == "auto") {
    facet_by <- if (length(levels(df$scenario_display)) > 1) "both" else "family"
  }

  base <- ggplot2::ggplot(df, ggplot2::aes(x = modname_ord, y = power, fill = family)) +
    ggplot2::geom_col(width = 0.72) +
    { if (!is.null(Its)) ggplot2::geom_errorbar(ggplot2::aes(ymin = lo95, ymax = hi95), width = 0, alpha = 0.45) } +
    { if (!is.null(Its)) ggplot2::geom_errorbar(ggplot2::aes(ymin = lo68, ymax = hi68), width = 0, alpha = 0.8, linewidth = 1) } +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100*x), "%"), limits = c(0,1)) +
    ggplot2::labs(x = NULL, y = "Power", fill = NULL, caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none",
                   strip.text = ggplot2::element_text(face = "bold"))

  p <- switch(facet_by,
              family   = base + ggplot2::facet_wrap(~ family, ncol = 2, scales = "free_y"),
              scenario = base + ggplot2::facet_wrap(~ scenario_display, ncol = 2, scales = "free_y"),
              both     = base + ggplot2::facet_grid(rows = ggplot2::vars(family),
                                                    cols  = ggplot2::vars(scenario_display), scales = "free_y")
  )
  p
}



##' Pretty gt table for power (grouped by *scenario* or *family*)
##'
##' @param tbl tibble from rgr_summarize() or rgr_run_scenarios()
##' @param digits number of digits for power
##' @param group_by which column to use as row-group header: "scenario" (default) or "family"
##' @param scenario_abbrev optional named character vector mapping
##'   full scenario names -> abbreviations (e.g., c("Middle inversion Δ=1, n=5"="MID"))
##'   If supplied, the table will display the *abbreviation*; a note is added automatically.
##' @param notes optional character string to append as a source note
##'              (if omitted and scenario_abbrev is supplied, an auto note is built)
##' @return gt table object
##' @export
# rgr_gt_power <- function(tbl,
#                          digits = 3,
#                          group_by = c("scenario", "family"),
#                          scenario_abbrev = NULL,
#                          notes = NULL) {
#   requireNamespace("gt", quietly = TRUE)
#
#   group_by <- match.arg(group_by)
#   df <- dplyr::as_tibble(tbl)
#
#   if (!all(c("scenario","family","modname","power") %in% names(df))) {
#     stop("Expected columns: scenario, family, modname, power. Did you pass rgr_summarize()/rgr_run_scenarios() output?")
#   }
#
#   # Optional: replace scenario names with abbreviations for display
#   if (!is.null(scenario_abbrev)) {
#     df$scenario <- as.character(df$scenario)
#     look <- names(scenario_abbrev)
#     repl <- unname(scenario_abbrev[ df$scenario ])
#     df$scenario <- ifelse(is.na(repl), df$scenario, repl)
#   }
#
#   # Arrange for readability
#   if (group_by == "scenario") {
#     ord <- df |>
#       dplyr::arrange(.data$scenario, .data$family, dplyr::desc(.data$power))
#   } else {
#     ord <- df |>
#       dplyr::arrange(.data$family, .data$scenario, dplyr::desc(.data$power))
#   }
#
#   # Build the table, grouping by the chosen column
#   gt_tbl <- gt::gt(ord, groupname_col = group_by) |>
#     gt::fmt_number(columns = "power", decimals = digits) |>
#     gt::tab_options(
#       table.border.top.width = gt::px(0),
#       table.border.bottom.width = gt::px(0),
#       row.striping.include_table_body = TRUE
#     ) |>
#     gt::cols_label(
#       scenario = "Scenario",
#       family   = "Family",
#       modname  = "Model",
#       power    = "Power"
#     )
#
#   # Auto-note for scenario abbreviations, unless user supplies their own 'notes'
#   if (is.null(notes) && !is.null(scenario_abbrev)) {
#     auto <- paste(paste0(scenario_abbrev, " = ", names(scenario_abbrev)), collapse = "  ·  ")
#     gt_tbl <- gt_tbl |> gt::tab_source_note(auto)
#   } else if (!is.null(notes)) {
#     gt_tbl <- gt_tbl |> gt::tab_source_note(notes)
#   }
#
#   gt_tbl
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
                         show_scenario = TRUE) {
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
