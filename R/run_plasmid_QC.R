#' Run single-library plasmid QC analysis
#'
#' Reads `distribution_all_clones.csv` from the results folder produced by
#' [run_extraction()] and runs five analysis sections:
#'
#' 1. **Clone-level statistics** — total clones, matched %, above DIs %.
#' 2. **TXT report** — detailed read-count breakdown saved to disk.
#' 3. **Clone frequency histogram** — log-scale plot distinguishing noise
#'    from real clones.
#' 4. **shRNA representation 3-panel** — % UCIs / raw UCI count /
#'    % reads per shRNA.
#' 5. **Birthday-problem analysis** — for each target number of distinct
#'    UCIs per shRNA (m = 1..10): duplication rate, cells to nucleofect,
#'    and probability of recovery.
#'
#' @param results_dir      Path to the `results/` folder from [run_extraction()].
#'   Must contain `distribution_all_clones.csv`. Default `"results/"`.
#' @param DIs              DIs threshold — same value used in [run_extraction()].
#'   Clones with fewer reads are treated as artefacts. Default `300`.
#' @param transfect_clones Number of clones recovered in your **test**
#'   nucleofection experiment. Used to estimate efficiency. Default `100`.
#' @param transfect_cells  Number of cells nucleofected in that test.
#'   Default `2000000`.
#' @param nucleofect_cells Number of cells you **plan** to nucleofect in the
#'   actual experiment. Determines expected clones per shRNA (K).
#'   Default `2000000`.
#'
#' @details
#' **K is never a user input.** It is derived as:
#' \deqn{K = \frac{\text{nucleofect\_cells} \times \text{efficiency}}{N_{\text{shRNAs}}}}
#' where efficiency = `transfect_clones / transfect_cells`.
#'
#' **Birthday problem model.** Each recovered clone independently received
#' one plasmid at random from the pool. If a shRNA has *n* distinct UCIs in
#' the library and *K* clones are recovered for it, some will share a UCI
#' (duplicates). The expected number of distinct UCIs is:
#' \deqn{E[\text{distinct}] = n \left(1 - \left(1 - \frac{1}{n}\right)^K\right)}
#'
#' @return A named list (returned invisibly) containing:
#' \describe{
#'   \item{`stats`}{Clone-level summary statistics (list).}
#'   \item{`shrna`}{Per-shRNA data frame with UCI counts and birthday metrics.}
#'   \item{`birthday`}{Birthday analysis results including `results` table
#'     and per-shRNA `shrna` data frame.}
#' }
#'
#' @section Output files written to `results_dir`:
#' \describe{
#'   \item{`DIs_threshold_stats_DIs<N>.txt`}{Clone count breakdown.}
#'   \item{`clone_frequency_dist_DIs<N>.pdf/.jpg`}{Log-scale histogram.}
#'   \item{`shrna_representation_3panel_DIs<N>.pdf/.jpg`}{3-panel figure.}
#'   \item{`UCI_birthday_analysis_DIs<N>.txt`}{Birthday analysis report.}
#'   \item{`calc_cells_and_dup_DIs<N>.pdf/.jpg`}{Cells and dup rate per m.}
#'   \item{`calc_pct_dup_per_shrna_DIs<N>.pdf/.jpg`}{% dup per shRNA.}
#'   \item{`calc_recovery_prob_DIs<N>.pdf/.jpg`}{Recovery probability.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- run_plasmid_QC(
#'   results_dir      = "CATCHER1_new/results/",
#'   DIs              = 300,
#'   transfect_clones = 100,
#'   transfect_cells  = 2000000,
#'   nucleofect_cells = 2000000
#' )
#'
#' # Inspect birthday results programmatically
#' result$birthday$results
#' result$shrna
#' }
#'
#' @seealso [run_extraction()], [run_combined_QC()], [run_plasmid_plots()]
#' @export
run_plasmid_QC <- function(results_dir      = "results/",
                            DIs              = 300,
                            transfect_clones = 100,
                            transfect_cells  = 2000000,
                            nucleofect_cells = 2000000) {

  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  eff <- transfect_clones / transfect_cells

  message("=============================================================")
  message("  catcheR::run_plasmid_QC()")
  message("=============================================================")
  message(sprintf("  Results dir        : %s", results_dir))
  message(sprintf("  DIs                : %d", DIs))
  message(sprintf("  Efficiency         : %d per %s cells",
                  transfect_clones, format(transfect_cells, big.mark = ",")))
  message(sprintf("  Planned nucleofection: %s cells",
                  format(nucleofect_cells, big.mark = ",")))

  # Load
  clones <- .load_clones(results_dir)

  # --- Section 1: stats -------------------------------------------------------
  message("\n-- Section 1: clone-level statistics --")
  s1 <- .clone_stats(clones, DIs)
  message(sprintf("  Total: %s | Above DIs: %s (%.1f%%) | Matched: %s (%.1f%%)",
                  format(s1$total,         big.mark = ","),
                  format(s1$above_DIs,     big.mark = ","),
                  s1$above_DIs / s1$total * 100,
                  format(s1$matched,       big.mark = ","),
                  s1$matched   / s1$total * 100))

  # --- Section 2: TXT report --------------------------------------------------
  message("\n-- Section 2: TXT report --")
  notes_map <- c("1" = "almost certainly artefacts", "2-10" = "likely noise",
                 "11-50" = "", "51-100" = "", "101-500" = "", "501-1,000" = "")
  txt_out <- file.path(results_dir,
                       paste0("DIs_threshold_stats_DIs", DIs, ".txt"))
  sink(txt_out)
  cat("=============================================================\n")
  cat("  DIs THRESHOLD STATISTICS\n")
  cat("=============================================================\n")
  cat(sprintf("  DIs threshold  : %s reads\n", format(DIs, big.mark = ",")))
  cat("-------------------------------------------------------------\n\n")
  cat("OVERVIEW\n")
  cat(sprintf("  Total clones         : %s\n",  format(s1$total,         big.mark=",")))
  cat(sprintf("  Matched barcodes     : %s  (%.1f%%)\n",
              format(s1$matched,   big.mark=","), s1$matched   / s1$total * 100))
  cat(sprintf("  Unmatched            : %s  (%.1f%%)\n",
              format(s1$unmatched, big.mark=","), s1$unmatched / s1$total * 100))
  cat(sprintf("  Below DIs (< %s)    : %s  (%.1f%%)  -- filtered out\n",
              format(DIs, big.mark=","),
              format(s1$below_DIs, big.mark=","), s1$below_DIs / s1$total * 100))
  cat(sprintf("  Above DIs (>= %s)   : %s  (%.1f%%)  -- kept\n",
              format(DIs, big.mark=","),
              format(s1$above_DIs, big.mark=","), s1$above_DIs / s1$total * 100))
  cat(sprintf("    of which matched   : %s  (%.1f%%)\n",
              format(s1$above_matched, big.mark=","),
              s1$above_matched / s1$above_DIs * 100))
  cat("\n-------------------------------------------------------------\n")
  cat("READ COUNT BREAKDOWN\n\n")
  cat(sprintf("  %-20s  %10s  %7s  %s\n",
              "Read count range", "Clones", "%", "Notes"))
  cat(sprintf("  %s\n", strrep("-", 60)))
  for (i in seq_len(nrow(s1$bc))) {
    lbl  <- as.character(s1$bc$bucket[i])
    note <- ifelse(lbl %in% names(notes_map), notes_map[[lbl]], "")
    cat(sprintf("  %-20s  %10s  %6.1f%%  %s\n",
                lbl, format(s1$bc$n_clones[i], big.mark = ","),
                s1$bc$pct[i], note))
  }
  cat("\n=============================================================\n")
  sink()
  message("  Saved: DIs_threshold_stats_DIs", DIs, ".txt")

  # --- Section 3: histogram ---------------------------------------------------
  message("\n-- Section 3: clone frequency histogram --")
  cl3 <- s1$clones %>%
    dplyr::mutate(match_status = dplyr::if_else(
      name == "UNMATCHED", "Unmatched", "Matched"))
  p3 <- ggplot2::ggplot(cl3, ggplot2::aes(x = Freq, fill = match_status)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count)),
                            bins = 60, colour = NA, alpha = 0.9) +
    ggplot2::geom_vline(xintercept = DIs, colour = "#D85A30",
                        linewidth = 0.8, linetype = "dashed") +
    ggplot2::annotate("text", x = DIs * 1.15, y = Inf,
                      label = paste0("DIs = ", format(DIs, big.mark = ",")),
                      hjust = 0, vjust = 1.5, size = 3.5, colour = "#D85A30") +
    ggplot2::scale_x_log10(labels = scales::label_comma(),
                            breaks = 10^(0:6)) +
    ggplot2::scale_y_log10(labels = scales::label_comma()) +
    ggplot2::scale_fill_manual(
      values = c("Matched" = "#378ADD", "Unmatched" = "#B4B2A9"), name = NULL) +
    ggplot2::labs(x = "Reads per clone", y = "Number of clones (log scale)",
                  title = paste0("Clone frequency distribution  (DIs=", DIs, ")")) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "top",
                   plot.caption = ggplot2::element_text(colour = "grey50"))
  .save_plot(p3, results_dir, paste0("clone_frequency_dist_DIs", DIs))

  # --- Section 4: 3-panel -----------------------------------------------------
  message("\n-- Section 4: shRNA representation 3-panel --")
  above <- s1$clones %>% dplyr::filter(Freq >= DIs, name != "UNMATCHED")
  shrna <- .shrna_agg(above)
  N     <- nrow(shrna)

  sw_uci   <- stats::shapiro.test(shrna$pct_UCI)
  sw_reads <- stats::shapiro.test(shrna$pct_reads)
  sw_n     <- stats::shapiro.test(shrna$n_UCI)
  wx_uci   <- stats::wilcox.test(shrna$pct_UCI,   mu = 100 / N)
  wx_reads <- stats::wilcox.test(shrna$pct_reads, mu = 100 / N)

  fill_col <- "#BA7517"
  bt <- .base_theme()

  make_hist_panel <- function(vals, total_n, x_lab, title_letter,
                               norm_p, med_p, bw = 0.5) {
    df <- data.frame(x = vals)
    ggplot2::ggplot(df, ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = ggplot2::after_stat(count) / total_n * 100),
        binwidth = bw, fill = fill_col, colour = NA, boundary = 0) +
      ggplot2::geom_vline(xintercept = stats::median(vals),
                          linetype = "dashed", linewidth = 0.6) +
      ggplot2::annotate("text", x = max(vals) * 0.6, y = Inf,
                        vjust = 1.5, hjust = 0, size = 3, colour = "grey30",
                        lineheight = 1.3,
                        label = paste0("Normality\n", norm_p,
                                       "\n\nMedian\n", med_p)) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.06))) +
      ggplot2::labs(x = x_lab, y = "% of shRNAs", title = title_letter) + bt
  }

  p_A <- make_hist_panel(shrna$pct_UCI, N, "% of UCIs per shRNA",
                          "A", .fmt_p(sw_uci$p.value),   .fmt_p(wx_uci$p.value))
  p_B <- ggplot2::ggplot(shrna, ggplot2::aes(x = n_UCI)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(count) / N * 100),
      binwidth = 2, fill = fill_col, colour = NA, boundary = 0.5) +
    ggplot2::geom_vline(xintercept = stats::median(shrna$n_UCI),
                        linetype = "dashed", linewidth = 0.6) +
    ggplot2::annotate("text", x = max(shrna$n_UCI) * 0.6, y = Inf,
                      vjust = 1.5, hjust = 0, size = 3, colour = "grey30",
                      lineheight = 1.3,
                      label = paste0("Normality\n", .fmt_p(sw_n$p.value),
                                     sprintf("\n\nMedian: %g UCIs",
                                             stats::median(shrna$n_UCI)))) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max(shrna$n_UCI) + 2, by = 4),
      expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.06))) +
    ggplot2::labs(x = "Number of UCIs per shRNA", y = "% of shRNAs",
                  title = "B") + bt
  p_C <- make_hist_panel(shrna$pct_reads, N, "% of reads per shRNA",
                          "C", .fmt_p(sw_reads$p.value), .fmt_p(wx_reads$p.value))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- (p_A | p_B | p_C) +
      patchwork::plot_annotation(
        title    = paste0("shRNA library representation  (DIs=", DIs, ")"),
        subtitle = sprintf("%d shRNAs | %d UCIs | mean=%.1f | median=%.0f",
                           N, sum(shrna$n_UCI), mean(shrna$n_UCI),
                           stats::median(shrna$n_UCI)),
        theme = ggplot2::theme(
          plot.title    = ggplot2::element_text(size = 13, face = "plain"),
          plot.subtitle = ggplot2::element_text(size = 10, colour = "grey50")))
    .save_plot(combined, results_dir,
               paste0("shrna_representation_3panel_DIs", DIs), w = 18, h = 5)
  } else {
    message("  Install patchwork for the combined 3-panel figure:")
    message("  install.packages('patchwork')")
    for (pp in list(list(p = p_A, n = paste0("shrna_rep_A_DIs", DIs)),
                    list(p = p_B, n = paste0("shrna_rep_B_DIs", DIs)),
                    list(p = p_C, n = paste0("shrna_rep_C_DIs", DIs))))
      .save_plot(pp$p, results_dir, pp$n, w = 7, h = 5)
  }

  # --- Section 5: birthday ----------------------------------------------------
  message("\n-- Section 5: birthday problem analysis --")
  b <- .birthday_results(shrna, N, DIs, eff,
                          nucleofect_cells, transfect_clones, transfect_cells)

  # Write birthday TXT
  txt_b <- file.path(results_dir,
                     sprintf("UCI_birthday_analysis_DIs%d.txt", DIs))
  sink(txt_b)
  cat("=============================================================\n")
  cat("  UCI UNIQUENESS ANALYSIS - BIRTHDAY PROBLEM\n")
  cat("=============================================================\n")
  cat(sprintf("  Efficiency        : %d per %s cells\n",
              transfect_clones, format(transfect_cells, big.mark = ",")))
  cat(sprintf("  Planned nucleofection : %s cells\n",
              format(nucleofect_cells, big.mark = ",")))
  cat(sprintf("  Expected K/shRNA  : %.2f clones\n", b$K_planned))
  cat(sprintf("  Library           : %d shRNAs, mean=%.1f UCIs\n\n",
              N, mean(shrna$n_UCI)))
  cat("RESULTS PER TARGET m\n\n")
  cat(sprintf("  %-4s  %-6s  %-10s  %-8s  %-12s  %-13s  %s\n",
              "m", "K", "cells(M)", "% dup", "P(>=m UCI)", "% shRNAs OK",
              "verdict"))
  cat(sprintf("  %s\n", strrep("-", 72)))
  for (i in seq_len(nrow(b$results))) {
    r <- b$results[i, ]
    cat(sprintf("  %-4d  %-6d  %-10.1f  %-7.1f%%  %-11.1f%%  %-12.1f%%  %s\n",
                r$m_target, r$K, r$cells_M, r$pct_dup_mean,
                r$p_recover_mean, r$pct_shrnas_ok, r$verdict))
  }
  cat("\n-------------------------------------------------------------\n")
  cat("PER-shRNA TABLE\n\n")
  cat(sprintf("  %-14s %-10s %6s %8s %10s  %s\n",
              "shRNA", "gene", "n_UCI", "% dup", "E[distinct]", "status"))
  cat(sprintf("  %s\n", strrep("-", 62)))
  for (i in seq_len(nrow(b$shrna %>% dplyr::arrange(n_UCI)))) {
    r <- (b$shrna %>% dplyr::arrange(n_UCI))[i, ]
    cat(sprintf("  %-14s %-10s %6d %7.1f%% %10.2f  %s\n",
                r$name, r$gene, r$n_UCI,
                r$pct_dup_planned, r$exp_distinct_plan, r$dup_status))
  }
  cat("\n=============================================================\n")
  sink()
  message("  Saved: UCI_birthday_analysis_DIs", DIs, ".txt")

  # Birthday plots
  rp <- b$results %>%
    dplyr::mutate(
      m_label  = paste0("m=", m_target),
      bar_fill = dplyr::case_when(
        pct_dup_mean <= 20 & pct_shrnas_ok >= 80 ~ "OK (<20% dup, lib OK)",
        pct_dup_mean <= 20                        ~ "Library insufficient",
        TRUE                                       ~ ">20% duplicates"
      )
    )

  p_cells <- ggplot2::ggplot(
    rp, ggplot2::aes(x = factor(m_label, levels = paste0("m=", 1:10)),
                     y = cells_M, fill = bar_fill)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%\ndup\n%d%%\nrecov.",
                                    pct_dup_mean, round(p_recover_mean))),
      vjust = -0.15, size = 2.7, lineheight = 1.1, colour = "grey20") +
    ggplot2::geom_hline(yintercept = nucleofect_cells / 1e6,
                        linetype = "dotted", colour = "black", linewidth = 0.9) +
    ggplot2::annotate("text", x = 0.6, y = nucleofect_cells / 1e6,
                      hjust = 0, vjust = -0.5, size = 3, colour = "black",
                      label = sprintf("Your plan: %.0fM cells (K=%.1f/shRNA)",
                                      nucleofect_cells / 1e6, b$K_planned)) +
    ggplot2::scale_fill_manual(
      values = c("OK (<20% dup, lib OK)" = "#1D9E75",
                 "Library insufficient"  = "#BA7517",
                 ">20% duplicates"       = "#D85A30"),
      name = NULL) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(suffix = "M"),
      expand = ggplot2::expansion(mult = c(0, 0.20))) +
    ggplot2::labs(x = "Target distinct UCIs per shRNA (m)",
                  y = "Cells to nucleofect (millions)",
                  title = "Cells needed and duplication rate per target m") +
    .base_theme() + ggplot2::theme(legend.position = "top")
  .save_plot(p_cells, results_dir,
             sprintf("calc_cells_and_dup_DIs%d", DIs), w = 12, h = 7)

  so <- b$shrna %>% dplyr::arrange(n_UCI) %>%
    dplyr::mutate(name = factor(name, levels = name))
  p_dup <- ggplot2::ggplot(
    so, ggplot2::aes(x = name, y = pct_dup_planned, fill = dup_status)) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::geom_hline(yintercept = 20, linetype = "dashed",
                        colour = "black", linewidth = 0.7) +
    ggplot2::scale_fill_manual(
      values = c("GOOD  (<=20% dup)" = "#1D9E75",
                 "FAIR  (<=40% dup)" = "#BA7517",
                 "POOR  (>40% dup)"  = "#D85A30"), name = NULL) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(suffix = "%"),
      limits = c(0, max(b$shrna$pct_dup_planned) * 1.15),
      expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(x = "shRNA (ordered by n UCIs)",
                  y = "% duplicate clones",
                  title = "% duplicate clones per shRNA at planned nucleofection") +
    .base_theme() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 90, hjust = 1,
                                               vjust = 0.5, size = 6),
      legend.position = "top")
  .save_plot(p_dup, results_dir,
             sprintf("calc_pct_dup_per_shrna_DIs%d", DIs), w = 16, h = 6)

  rec_df <- dplyr::bind_rows(lapply(2:5, function(m) {
    b$shrna %>%
      dplyr::mutate(
        m_target = m, m_label = paste0("m=", m),
        p_rec    = sapply(n_UCI, function(n) .p_at_least_m(n, m, m) * 100),
        lib_ok   = n_UCI >= m)
  }))
  p_rec <- ggplot2::ggplot(
    rec_df %>% dplyr::mutate(
      name = factor(name,
                    levels = b$shrna %>% dplyr::arrange(n_UCI) %>%
                      dplyr::pull(name))),
    ggplot2::aes(x = name, y = p_rec, fill = lib_ok)) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::geom_hline(yintercept = 80, linetype = "dashed",
                        colour = "black", linewidth = 0.6) +
    ggplot2::facet_wrap(~m_label, ncol = 2) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#1D9E75", "FALSE" = "#D85A30"),
      labels = c("TRUE" = "n >= m (feasible)", "FALSE" = "n < m (impossible)"),
      name = NULL) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(suffix = "%"),
      limits = c(0, 105), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(x = "shRNA (ordered by UCI count)",
                  y = "P(recover >= m distinct UCIs) %",
                  title = "Probability of recovering >= m distinct UCIs per shRNA") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                           vjust = 0.5, size = 4),
      legend.position = "top",
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA))
  .save_plot(p_rec, results_dir,
             sprintf("calc_recovery_prob_DIs%d", DIs), w = 14, h = 10)

  message("\n=============================================================")
  message("  Done. All outputs written to: ", results_dir)
  message("=============================================================")

  invisible(list(stats = s1, shrna = shrna, birthday = b))
}
