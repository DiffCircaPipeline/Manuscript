# the original dodr algorithm forces to ouptut q-value, but here we need p-values
compareRhythms_dodr2 <- function(expr, exp_design, period=24, rhythm_fdr = 0.05,
                                compare_fdr = 0.05, amp_cutoff = 0.5,
                                just_classify = TRUE) {

  exp_design <- base::cbind(exp_design,
                            col_number = base::seq(base::nrow(exp_design)))

  group_id <- base::levels(exp_design$group)
  assertthat::are_equal(length(group_id), 2)

  exp_design_A <- exp_design[exp_design$group == group_id[1], ]

  exp_design_A <- exp_design_A[base::order(exp_design_A$time), ]

  exp_design_B <- exp_design[exp_design$group == group_id[2], ]

  exp_design_B <- exp_design_B[base::order(exp_design_B$time), ]

  deltat_A <- min(diff(base::unique(exp_design_A$time)))

  time_A <- base::seq(min(exp_design_A$time),
                      max(exp_design_A$time),
                      by = deltat_A)

  unique_times_A <- base::table(exp_design_A$time)

  measure_sequence_A <- base::vapply(time_A,
                                     function(t) {
                                       ifelse(any(names(unique_times_A) == t),
                                              unique_times_A[names(unique_times_A) == t],
                                              0)
                                     },
                                     integer(1))


  deltat_B <- min(diff(base::unique(exp_design_B$time)))

  time_B <- base::seq(min(exp_design_B$time),
                      max(exp_design_B$time),
                      by = deltat_B)

  unique_times_B <- base::table(exp_design_B$time)

  measure_sequence_B <- base::vapply(time_B,
                                     function(t) {
                                       ifelse(any(names(unique_times_B) == t),
                                              unique_times_B[names(unique_times_B) == t],
                                              0)
                                     },
                                     integer(1))


  expr_A <- expr[, exp_design_A$col_number]
  expr_B <- expr[, exp_design_B$col_number]

  assertthat::assert_that((sum(measure_sequence_A) >= 12) && (sum(measure_sequence_B) >= 12),
                          msg = "Not enough samples to run RAIN rhythmicity analysis confidently.")

  rain_A <- rain::rain(t(expr_A), deltat_A, period,
                       measure.sequence = measure_sequence_A)

  rain_B <- rain::rain(t(expr_B), deltat_B, period,
                       measure.sequence = measure_sequence_B)

  # should change this part
  # rain_results <- data.frame(stats::p.adjust(rain_A[, "pVal"], method = "BH"),
  #                            stats::p.adjust(rain_B[, "pVal"], method = "BH"))
  # colnames(rain_results) <- c("adj_p_val_A", "adj_p_val_B")

  rain_results <- data.frame(rain_A[, "pVal"],
                             rain_B[, "pVal"])
  colnames(rain_results) <- c("p_val_A", "p_val_B")

  compute_circ_params <- function(y, t, period) {

    inphase <- cos(2 * pi * t / period)
    outphase <- sin(2 * pi * t / period)

    X <- stats::model.matrix(~inphase + outphase)

    fit <- stats::lm.fit(X, t(y))

    amps <- 2 * sqrt(base::colSums(fit$coefficients[-1, ]^2))

    phases <- (atan2(fit$coefficients[3, ], fit$coefficients[2, ]) %% (2*pi))

    return(base::cbind(amps = amps, phases = phases))

  }
  circ_params_A <- compute_circ_params(expr_A, exp_design_A$time, period = period)

  circ_params_B <- compute_circ_params(expr_B, exp_design_B$time, period = period)

  rhythmic_in_A <- (rain_results$p_val_A < rhythm_fdr) &
    (circ_params_A[, "amps"] > amp_cutoff)

  rhythmic_in_B <- (rain_results$p_val_B < rhythm_fdr) &
    (circ_params_B[, "amps"] > amp_cutoff)

  rhythmic_in_either <- rhythmic_in_A | rhythmic_in_B

  assertthat::assert_that(sum(rhythmic_in_either) > 0,
                          msg = "Sorry no rhythmic genes in either dataset for the thresholds provided.")

  dodr_results <- DODR::robustDODR(t(expr_A[rhythmic_in_either, ]),
                                   t(expr_B[rhythmic_in_either, ]),
                                   times1 = exp_design_A$time,
                                   times2 = exp_design_B$time,
                                   norm = TRUE,
                                   period = period)
  # dodr_results$adj_p_val <- stats::p.adjust(dodr_results$p.value, method = "BH")

  # results <- data.frame(id = rownames(expr_A)[rhythmic_in_either],
  #                       rhythmic_in_A = rhythmic_in_A[rhythmic_in_either],
  #                       rhythmic_in_B = rhythmic_in_B[rhythmic_in_either],
  #                       diff_rhythmic = dodr_results$adj_p_val < compare_fdr,
  #                       stringsAsFactors = FALSE)
  results <- data.frame(id = rownames(expr_A)[rhythmic_in_either],
                        rhythmic_in_A = rhythmic_in_A[rhythmic_in_either],
                        rhythmic_in_B = rhythmic_in_B[rhythmic_in_either],
                        diff_rhythmic = dodr_results$p.value < compare_fdr,
                        stringsAsFactors = FALSE)

  results$category <- base::mapply(categorize,
                                   results$rhythmic_in_A,
                                   results$rhythmic_in_B,
                                   results$diff_rhythmic)
  results <- results[, c(1, 5, 2, 3, 4)]

  if (!just_classify) {
    expand_results <- data.frame(
      A_amp = circ_params_A[rhythmic_in_either, "amps"],
      A_phase = circ_params_A[rhythmic_in_either, "phases"],
      B_amp = circ_params_B[rhythmic_in_either, "amps"],
      B_phase = circ_params_B[rhythmic_in_either, "phases"],
      p_val_A = rain_results$p_val_A[rhythmic_in_either],
      p_val_B = rain_results$p_val_B[rhythmic_in_either],
      p_val_dodr = dodr_results$p.value
    )
    results <- base::cbind(results, expand_results)
  }

  rownames(results) <- NULL
  colnames(results) <- gsub("A", group_id[1], colnames(results))
  colnames(results) <- gsub("B", group_id[2], colnames(results))

  return(results)
}

compareRhythms2 <- function(data, exp_design, lengths=NULL,
                           method = "mod_sel", period=24, rhythm_fdr = 0.05,
                           compare_fdr = 0.05, amp_cutoff = 0.5,
                           criterion = "bic", schwarz_wt_cutoff = 0.6,
                           just_classify = TRUE, robust = TRUE, outliers = FALSE
) {

  assertthat::assert_that(
    is.matrix(data),
    assertthat::not_empty(data),
    is.data.frame(exp_design),
    assertthat::not_empty(exp_design),
    assertthat::are_equal(ncol(data), nrow(exp_design)),
    assertthat::has_name(exp_design, c("time", "group")),
    assertthat::noNA(data),
    method %in% c("mod_sel", "limma", "dodr", "voom", "deseq2", "edger"),
    is.factor(exp_design$group),
    is.numeric(exp_design$time),
    length(levels(exp_design$group)) == 2,
    assertthat::is.number(period),
    period > 0,
    assertthat::is.number(rhythm_fdr),
    rhythm_fdr <= 1.0 & rhythm_fdr >= 0,
    assertthat::is.number(compare_fdr),
    compare_fdr <= 1.0 & compare_fdr >= 0,
    assertthat::is.number(amp_cutoff),
    amp_cutoff >= 0,
    assertthat::is.string(criterion),
    criterion %in% c("bic", "aic"),
    assertthat::is.number(schwarz_wt_cutoff),
    schwarz_wt_cutoff >=0 & schwarz_wt_cutoff <=1.0,
    assertthat::is.flag(just_classify),
    assertthat::is.flag(robust),
    assertthat::is.flag(outliers)
  )
  if (method %in% c("deseq", "edger") && !is.null(lengths)) {
    assertthat::assert_that(all(lengths>0), msg = "All transcript lengths are not positive")
  }

  if (assertthat::has_name(exp_design, "batch")) {
    assertthat::assert_that(is.factor(exp_design$batch))
  }

  switch (method,
          mod_sel = compareRhythms_model_select(data = data,
                                                exp_design = exp_design,
                                                period = period,
                                                amp_cutoff = amp_cutoff,
                                                criterion = criterion,
                                                schwarz_wt_cutoff = schwarz_wt_cutoff,
                                                just_classify = just_classify),
          dodr = compareRhythms_dodr2(expr = data,
                                     exp_design = exp_design,
                                     period = period,
                                     rhythm_fdr = rhythm_fdr,
                                     compare_fdr = compare_fdr,
                                     amp_cutoff = amp_cutoff,
                                     just_classify = just_classify),
          limma = compareRhythms_limma(eset = data,
                                       exp_design = exp_design,
                                       period = period,
                                       rhythm_fdr = rhythm_fdr,
                                       amp_cutoff = amp_cutoff,
                                       compare_fdr = compare_fdr,
                                       just_classify = just_classify,
                                       rna_seq = FALSE,
                                       robust = robust),
          voom = compareRhythms_voom(counts = data,
                                     exp_design = exp_design,
                                     period = period,
                                     rhythm_fdr = rhythm_fdr,
                                     amp_cutoff = amp_cutoff,
                                     compare_fdr = compare_fdr,
                                     just_classify = just_classify,
                                     robust = robust,
                                     outliers = outliers),
          deseq2 = compareRhythms_deseq2(counts = data, exp_design = exp_design,
                                         lengths = lengths, period = period, rhythm_fdr = rhythm_fdr,
                                         amp_cutoff = amp_cutoff, compare_fdr = compare_fdr,
                                         just_classify = just_classify),
          edger = compareRhythms_edgeR(counts = data, exp_design = exp_design,
                                       lengths = lengths, period = period, rhythm_fdr = rhythm_fdr,
                                       amp_cutoff = amp_cutoff, compare_fdr = compare_fdr,
                                       just_classify = just_classify)
  )
}

compute_circ_params <- function(y, t, period) {

  inphase <- cos(2 * pi * t / period)
  outphase <- sin(2 * pi * t / period)

  X <- stats::model.matrix(~inphase + outphase)

  fit <- stats::lm.fit(X, t(y))

  amps <- 2 * sqrt(base::colSums(fit$coefficients[-1, ]^2))

  phases <- (atan2(fit$coefficients[3, ], fit$coefficients[2, ]) %% (2*pi))

  return(base::cbind(amps = amps, phases = phases))

}

compute_model_params <- function(y, group_id, d=NULL, type="fit") {

  if (type == "fit") {
    fit <- limma::lmFit(y, d)
    coeffs <- fit$coefficients
  }
  else if (type == "coef") {
    coeffs <- y
  }

  if (any(base::grepl(paste0(group_id[1], "_"),
                      colnames(coeffs)))) {
    rhy_params <- coeffs[, base::paste(group_id[1],
                                       c("inphase", "outphase"),
                                       sep = "_")]
    amps_A <- 2 * sqrt(base::rowSums(rhy_params^2))
    phases_A <- base::atan2(rhy_params[, 2], rhy_params[, 1]) %% (2*pi)
  } else {
    amps_A <- 0
    phases_A <- 0
  }

  if (any(base::grepl(paste0(group_id[2], "_"),
                      colnames(coeffs)))) {
    rhy_params <- coeffs[, base::paste(group_id[2],
                                       c("inphase", "outphase"),
                                       sep = "_")]
    amps_B <- 2 * sqrt(base::rowSums(rhy_params^2))
    phases_B <- base::atan2(rhy_params[, 2], rhy_params[, 1])  %% (2*pi)
  } else {
    amps_B <- 0
    phases_B <- 0
  }

  if (all(base::is.element(c("inphase", "outphase"),
                           colnames(coeffs)))) {
    rhy_params <- coeffs[, c("inphase", "outphase")]
    amps <- 2 * sqrt(base::rowSums(rhy_params^2))
    phases <- base::atan2(rhy_params[, 2], rhy_params[, 1])  %% (2*pi)
    model_params <- base::cbind(amps, phases, amps, phases)
  } else {
    model_params <- base::cbind(amps_A, phases_A, amps_B, phases_B)
  }

  colnames(model_params) <- base::paste(rep(group_id, each = 2),
                                        c("amp", "phase"), sep = "_")

  return(model_params)
}

categorize <- function(a, b, dr) {
  if (a && !b && dr) {
    category <- "loss"
  } else if (!a && b && dr) {
    category <- "gain"
  } else if (a && b && dr) {
    category <- "change"
  } else {
    category <- "same"
  }
  return(category)
}
