#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(MuMIn)
  library(car)
})

# -------------------------
# Paths
# -------------------------
in_csv <- "/users/kcl/code/ES/dfm_for_nlme_fixed900ms.csv"
false_count_src_csv <- "/users/kcl/code/ES/dfm_for_lme4_final_with_drowsy_fixed900ms.csv"
sig_coef_csv <- "/users/kcl/code/ES/3_pc_lmer_coefficients_FDR.csv"
output_dir <- getwd()
out_suffix <- "_ver4"

out_file <- function(stem) {
  file.path(output_dir, paste0(stem, out_suffix, ".csv"))
}

cat("Results will be saved to:", output_dir, "\n")
cat("False-count source CSV:", false_count_src_csv, "\n")
cat("[INFO] Input drowsy source is 900 ms window (cue-0.4s to cue+0.5s).\n")

build_false_response_map <- function(csv_path) {
  if (!file.exists(csv_path)) {
    cat("[WARN] false_response_count source not found:", csv_path, "\n")
    return(NULL)
  }

  x <- read.csv(csv_path, stringsAsFactors = FALSE)
  req <- c("subject_id", "trial", "false_response_count", "trial_ear_success_frame_count")
  miss <- setdiff(req, names(x))
  if (length(miss) > 0) {
    cat("[WARN] false_response_count source missing columns:", paste(miss, collapse = ", "), "\n")
    return(NULL)
  }

  x$subject_id_chr <- as.character(x$subject_id)
  x$trial_num <- suppressWarnings(as.integer(as.character(x$trial)))
  x$false_response_count <- suppressWarnings(as.numeric(as.character(x$false_response_count)))
  x$trial_ear_success_frame_count <- suppressWarnings(
    as.numeric(as.character(x$trial_ear_success_frame_count))
  )

  out <- x %>%
    select(subject_id_chr, trial_num, false_response_count, trial_ear_success_frame_count) %>%
    filter(!is.na(subject_id_chr), !is.na(trial_num)) %>%
    group_by(subject_id_chr, trial_num) %>%
    summarise(
      false_response_count = if (all(!is.finite(false_response_count))) NA_real_
      else max(false_response_count[is.finite(false_response_count)]),
      trial_ear_success_frame_count = if (all(!is.finite(trial_ear_success_frame_count))) NA_real_
      else max(trial_ear_success_frame_count[is.finite(trial_ear_success_frame_count)]),
      .groups = "drop"
    )

  as.data.frame(out)
}

# -------------------------
# 1) Load & Basic Preprocessing
# -------------------------
dfm <- read.csv(in_csv, stringsAsFactors = FALSE)
# -------------------------
# [ADD] Filter subjects (exclude-list)
# -------------------------
excluded_subjects <- c(
  1013, 1025, 1032, 1038, 1053, 1061, 1063, 1077
)
# excluded_subjects <- c(
#   1012, 1021, 1025, 1032, 1038, 1053, 1054, 1061, 1063, 1077
# )
# subject_id가 숫자/문자 혼재할 수 있어서 안전하게 숫자로 비교
if (!"subject_id" %in% names(dfm)) stop("[ERROR] subject_id column not found in dfm.")

dfm$subject_id_num <- suppressWarnings(as.integer(as.character(dfm$subject_id)))
before_n <- nrow(dfm)

dfm <- dfm[!is.na(dfm$subject_id_num) & !(dfm$subject_id_num %in% excluded_subjects), ]
dfm$subject_id_num <- NULL

cat(
  "[INFO] Subject exclusion filter applied (n_excluded_ids=",
  length(excluded_subjects),
  "). Rows:",
  before_n,
  "->",
  nrow(dfm),
  "\n",
  sep = ""
)

# 혹시 필터 후 subject_id가 전부 NA가 되면(파싱 실패) 바로 알려주기
if (nrow(dfm) == 0) {
  stop("[ERROR] After filtering, dfm has 0 rows. Check subject_id format in CSV.")
}

if (!"subject_id" %in% names(dfm)) {
  stop("[ERROR] 'subject_id' column not found.")
}

if (!"log_rt" %in% names(dfm)) {
  if ("response_time" %in% names(dfm)) {
    cat("[INFO] 'log_rt' not found. Creating from 'response_time'...\n")
    dfm$log_rt <- ifelse(dfm$response_time > 0, log(dfm$response_time), NA_real_)
  } else {
    stop("[ERROR] Neither 'log_rt' nor 'response_time' columns found.")
  }
}

# Join keys for trial-level covariates
dfm$subject_id_chr <- as.character(dfm$subject_id)
dfm$trial_num <- suppressWarnings(as.integer(as.character(dfm$trial)))

# Keep only trials with false_response_count == 0
false_map <- build_false_response_map(false_count_src_csv)
if (!is.null(false_map)) {
  before_false_filter <- nrow(dfm)
  dfm <- dfm %>% left_join(false_map, by = c("subject_id_chr", "trial_num"))
  # Resolve success-frame column robustly (join can create .x/.y suffixes).
  if ("trial_ear_success_frame_count" %in% names(dfm)) {
    succ_vals <- suppressWarnings(as.numeric(as.character(dfm$trial_ear_success_frame_count)))
  } else {
    succ_x <- if ("trial_ear_success_frame_count.x" %in% names(dfm)) {
      suppressWarnings(as.numeric(as.character(dfm$trial_ear_success_frame_count.x)))
    } else {
      rep(NA_real_, nrow(dfm))
    }
    succ_y <- if ("trial_ear_success_frame_count.y" %in% names(dfm)) {
      suppressWarnings(as.numeric(as.character(dfm$trial_ear_success_frame_count.y)))
    } else {
      rep(NA_real_, nrow(dfm))
    }
    succ_vals <- succ_x
    fill <- (!is.finite(succ_vals)) & is.finite(succ_y)
    succ_vals[fill] <- succ_y[fill]
  }
  dfm$trial_ear_success_frame_count_qc <- succ_vals

  keep_false0 <- !is.na(dfm$false_response_count) & (dfm$false_response_count == 0)
  keep_ear_success <- is.finite(dfm$trial_ear_success_frame_count_qc) & (dfm$trial_ear_success_frame_count_qc > 0)
  keep_mask <- keep_false0 & keep_ear_success
  removed_false <- before_false_filter - sum(keep_mask)
  missing_false <- sum(is.na(dfm$false_response_count))
  missing_success <- sum(!is.finite(dfm$trial_ear_success_frame_count_qc))
  removed_zero_success <- sum(keep_false0 & !keep_ear_success, na.rm = TRUE)
  dfm <- dfm[keep_mask, , drop = FALSE]
  cat(
    "[INFO] Applied trial QC filter (false_response_count==0 & success_frames>0). Rows:",
    before_false_filter, "->", nrow(dfm),
    "(removed=", removed_false,
    ", removed_zero_success=", removed_zero_success,
    ", missing_false_map=", missing_false,
    ", missing_success_map=", missing_success, ")\n"
  )
} else {
  cat("[WARN] false_response_count filter skipped (source unavailable).\n")
}

drowsy_interact_targets <- c(
  # "trial_ear_norm_baseline",    # disabled: remove PC x EAR_baseline interaction
  "trial_ear_norm_baseline_sd"
)
missing_drowsy_targets <- setdiff(drowsy_interact_targets, names(dfm))
if (length(missing_drowsy_targets) > 0) {
  stop(
    "[ERROR] Required drowsy interaction columns not found in dfm: ",
    paste(missing_drowsy_targets, collapse = ", ")
  )
}

# -------------------------
# 2) Type Casting & Scaling
# -------------------------
dfm$subject_id <- factor(dfm$subject_id)

if ("sound_norm" %in% names(dfm)) {
  dfm$sound_norm <- factor(dfm$sound_norm)
  bgm_id_col <- "sound_norm"
} else {
  cat("[WARN] 'sound_norm' column not found. Random effect for BGM will be skipped.\n")
  bgm_id_col <- NULL
}

# block order covariate (one block per sound, ordered by first trial within subject)
if (!is.null(bgm_id_col) && ("trial_num" %in% names(dfm))) {
  dfm$sound_key <- as.character(dfm[[bgm_id_col]])
  block_map <- dfm %>%
    filter(!is.na(subject_id_chr), !is.na(trial_num), !is.na(sound_key)) %>%
    group_by(subject_id_chr, sound_key) %>%
    summarise(block_start_trial = min(trial_num), .groups = "drop") %>%
    arrange(subject_id_chr, block_start_trial, sound_key) %>%
    group_by(subject_id_chr) %>%
    mutate(block_order = row_number()) %>%
    ungroup() %>%
    select(subject_id_chr, sound_key, block_order)

  dfm <- dfm %>% left_join(block_map, by = c("subject_id_chr", "sound_key"))
  dfm$sound_key <- NULL
  cat("[INFO] block_order ready (non-NA):", sum(!is.na(dfm$block_order)), "/", nrow(dfm), "\n")
} else {
  cat("[WARN] block_order could not be created (missing sound_norm or trial).\n")
}

fac_cols <- c("cueType", "arrowDirection", "congruency", "targetDirection", "drive_exp", "drowsiness_label")
for (cc in fac_cols) {
  if (cc %in% names(dfm)) dfm[[cc]] <- factor(dfm[[cc]])
}

if ("eye_status" %in% names(dfm)) {
  dfm$eye_status <- factor(as.integer(dfm$eye_status))
}

snd_pc_cols <- grep("^snd_pc[0-9]+$", names(dfm), value = TRUE)
psycho_vars <- c("Sharpness", "Roughness", "Tonality")
mental_vars <- c("BDI_II", "BAI", "PSS")
drowsy_cont_vars <- intersect(c("eye_status_ratio"), names(dfm))

cont_to_scale <- intersect(
  c(
    "age", "sex", "satisfaction_filled", snd_pc_cols, psycho_vars,
    drowsy_cont_vars, mental_vars, "block_order"
  ),
  names(dfm)
)

for (cc in cont_to_scale) {
  dfm[[paste0(cc, "_z")]] <- as.numeric(scale(dfm[[cc]]))
}

# -------------------------
# 3) Build Base Terms (same backbone as main pc model)
# -------------------------
task_terms <- intersect(c("cueType", "arrowDirection", "congruency", "targetDirection"), names(dfm))

demo_terms <- c()
if ("age_z" %in% names(dfm)) demo_terms <- c(demo_terms, "age_z")
if ("sex_z" %in% names(dfm)) demo_terms <- c(demo_terms, "sex_z")
if ("drive_exp" %in% names(dfm)) demo_terms <- c(demo_terms, "drive_exp")
for (mv in c("BDI_II_z", "BAI_z", "PSS_z")) {
  if (mv %in% names(dfm)) demo_terms <- c(demo_terms, mv)
}

satis_terms <- if ("satisfaction_filled_z" %in% names(dfm)) "satisfaction_filled_z" else character(0)
snd_terms <- intersect(paste0(snd_pc_cols, "_z"), names(dfm))

base_fixed_terms <- unique(c(task_terms, demo_terms, satis_terms, snd_terms, drowsy_interact_targets))

# Additional confounds
confound_terms <- c()
if ("block_order_z" %in% names(dfm)) {
  confound_terms <- c(confound_terms, "block_order_z")
} else if ("block_order" %in% names(dfm)) {
  confound_terms <- c(confound_terms, "block_order")
}
base_fixed_terms <- unique(c(base_fixed_terms, confound_terms))

# # -------------------------
# # 4) Load FDR-significant terms from main pc model and map to formula vars
# # -------------------------
# if (!file.exists(sig_coef_csv)) {
#   stop("[ERROR] Significant-term file not found: ", sig_coef_csv)
# }

# sig_df <- read.csv(sig_coef_csv, stringsAsFactors = FALSE)
# required_cols <- c("Term", "sig_fdr_0.05")
# if (!all(required_cols %in% names(sig_df))) {
#   stop("[ERROR] Required columns missing in sig file. Need: Term, sig_fdr_0.05")
# }

# sig_flag <- sig_df[["sig_fdr_0.05"]]
# if (!is.logical(sig_flag)) sig_flag <- as.logical(sig_flag)

# sig_terms <- unique(sig_df$Term[which(sig_flag %in% TRUE)])
# sig_terms <- sig_terms[!is.na(sig_terms) & sig_terms != "(Intercept)"]

# if (length(sig_terms) == 0) {
#   stop("[ERROR] No sig_fdr_0.05 == TRUE terms found in 3_pc_lmer_coefficients_FDR.csv")
# }

# factor_candidates <- intersect(c("cueType", "arrowDirection", "congruency", "targetDirection", "drive_exp"), base_fixed_terms)

# map_sig_term <- function(term) {
#   if (term %in% base_fixed_terms) return(term)
#   for (fc in factor_candidates) {
#     if (startsWith(term, fc)) return(fc)
#   }
#   NA_character_
# }

# mapped_terms <- vapply(sig_terms, map_sig_term, character(1))
# unmapped_terms <- sig_terms[is.na(mapped_terms)]
# if (length(unmapped_terms) > 0) {
#   cat("[WARN] Unmapped significant terms (skipped):", paste(unmapped_terms, collapse = ", "), "\n")
# }

# interaction_vars <- unique(mapped_terms[!is.na(mapped_terms)])
# interaction_vars <- setdiff(interaction_vars, drowsy_interact_targets)
# interaction_vars <- interaction_vars[interaction_vars %in% base_fixed_terms]

# if (length(interaction_vars) == 0) {
#   stop("[ERROR] No valid variables mapped for drowsy interactions.")
# }

# cat("[INFO] FDR-significant terms from main model:\n")
# cat(" ", paste(sig_terms, collapse = ", "), "\n")
# cat("[INFO] Drowsy targets:\n")
# cat(" ", paste(drowsy_interact_targets, collapse = ", "), "\n")
# cat("[INFO] Variables used for interaction with drowsy targets:\n")
# cat(" ", paste(interaction_vars, collapse = ", "), "\n")

# -------------------------
# 4) [REPLACEMENT] Manually Define Interaction Variables
# -------------------------
# Basic 모델 유의성과 상관없이, 가설에 따라 PC 1~7축을 상호작용 분석에 포함합니다.
interaction_vars <- intersect(paste0(snd_pc_cols, "_z"), names(dfm))

if (length(interaction_vars) == 0) {
  stop("[ERROR] No PC columns found for interaction.")
}

cat("[INFO] Manually selected variables for interaction with drowsy targets (PC1-PC7):\n")
cat(" ", paste(interaction_vars, collapse = ", "), "\n")

# -------------------------
# 5) NA Removal based on actual model vars
# -------------------------
vars_needed <- unique(c("log_rt", "subject_id", bgm_id_col, base_fixed_terms, interaction_vars))
vars_needed <- vars_needed[!is.na(vars_needed) & vars_needed %in% names(dfm)]

dfm <- dfm[complete.cases(dfm[, vars_needed]), ]
dfm <- droplevels(dfm)
cat("After final NA removal:", nrow(dfm), "rows\n")

base_fixed_terms <- setdiff(base_fixed_terms, bgm_id_col)
# -------------------------
# 6) Formula Construction (only sig vars x [avg_ear_norm, avg_ear_norm_sd] interactions)
# -------------------------
interaction_terms <- unlist(
  lapply(
    drowsy_interact_targets, # avg_ear_norm 등
    function(dv) paste0(interaction_vars, " * ", dv)
  ),
  use.names = FALSE
)

# random_part <- "+ (1 | subject_id)"
# if (!is.null(bgm_id_col)) {
#   random_part <- paste0(random_part, " + (1 | ", bgm_id_col, ")")
# }
if (!is.null(bgm_id_col) && (bgm_id_col %in% names(dfm))) {
  random_part <- paste0(
    "+ (1 | subject_id) ",
    "+ (1 | subject_id:", bgm_id_col, ")"
  )
} else {
  random_part <- "+ (1 | subject_id)"
}

# all_fixed <- unique(c(task_terms, demo_terms, satis_terms, interaction_terms))
all_fixed <- unique(c(task_terms, satis_terms, confound_terms, interaction_terms))
fixed_formula <- paste(all_fixed, collapse = " + ")
fF_str <- paste("log_rt ~", fixed_formula, random_part)

cat("\n=== REVISED INTERACTION MODEL FORMULA ===\n", fF_str, "\n")

# -------------------------
# 7) Fit Model (lmer)
# -------------------------
mF <- lme4::lmer(as.formula(fF_str), data = dfm, REML = TRUE, 
                 control = lmerControl(optimizer = "bobyqa"))

# 4. P-value 추출: 원래 쓰시던 방식 그대로 (vcov 안 쓰고 직접 계산)
fe <- fixef(mF)
V_unsc <- as.matrix(mF@pp$unsc()) # C언어단에서 터지는 vcov() 대신 안전한 내부 행렬 사용
se <- sigma(mF) * sqrt(diag(V_unsc))
tvals <- fe / se
pvals <- 2 * pnorm(abs(tvals), lower.tail = FALSE) # 자유도 무한대 가정 (데이터가 많아서 OK)

coef_df <- data.frame(
  Term = names(fe),
  Estimate = unname(fe),
  StdError = se,
  tValue = tvals,
  pValue = pvals
)

cat("\n[OK] Interaction model fitted using lmer.\n")

# -------------------------
# VIF from fixed-effects lm proxy
# -------------------------
cat("\n[INFO] Calculating VIF from fixed-effects lm proxy...\n")
fixed_only_str <- paste("log_rt ~", fixed_formula) 
mF_lm <- lm(as.formula(fixed_only_str), data = dfm)
vif_res <- tryCatch(car::vif(mF_lm), error = function(e) NULL)

if (!is.null(vif_res)) {
  if (is.matrix(vif_res)) {
    vif_df <- data.frame(Term = rownames(vif_res), vif_res)
  } else {
    vif_df <- data.frame(Term = names(vif_res), VIF = vif_res)
  }
  write.csv(vif_df, out_file("0_pc_interaction_lmer_model_VIF"), row.names = FALSE)
  print(vif_df)
}
rm(mF_lm)

# -------------------------
# [추가] STEP 0) Model Diagnostics Plot (ggplot2 직접 구현)
# -------------------------
cat("\n[INFO] Generating diagnostic plots using ggplot2 (failsafe)...\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# 잔차와 적합값 추출 (vcov()를 피하는 안전한 방법)
diag_df <- data.frame(
  fitted = fitted(mF),
  resid  = residuals(mF)
)

# 1. 등분산성: Residuals vs Fitted
p_resid <- ggplot(diag_df, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.3, color = "#2C3E50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#E74C3C", linewidth = 1) +
  geom_smooth(method = "loess", se = FALSE, color = "#2980B9") +
  theme_minimal(base_size = 14) +
  labs(title = "Residuals vs Fitted", 
       subtitle = "Check for Homogeneity of Variance",
       x = "Fitted Values", 
       y = "Residuals")

# 2. 정규성: Q-Q Plot
p_qq <- ggplot(diag_df, aes(sample = resid)) +
  stat_qq(alpha = 0.3, color = "#2C3E50") +
  stat_qq_line(color = "#E74C3C", linetype = "dashed", linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(title = "Normal Q-Q Plot", 
       subtitle = "Check for Normality of Residuals",
       x = "Theoretical Quantiles", 
       y = "Sample Quantiles")

# 두 그래프 나란히 배치 (patchwork 패키지 활용)
p_final <- p_resid + p_qq

# 저장
plot_path <- file.path(output_dir, paste0("0_pc_interaction_diagnostics", out_suffix, ".png"))
ggsave(filename = plot_path, plot = p_final, width = 12, height = 6, dpi = 300, bg = "white")

cat("- Saved diagnostic plot:", plot_path, "\n")
# ---------------------------------------------------------

# -------------------------
# STEP 1) Fit Stats
# -------------------------
fit_stats <- data.frame(
  Model = "Full_Model_lmer_interaction_ver3",
  AIC = AIC(mF),
  BIC = BIC(mF),
  logLik = as.numeric(logLik(mF)),
  stringsAsFactors = FALSE
)

r2m <- NA_real_
r2c <- NA_real_
tryCatch({
  vc <- as.data.frame(VarCorr(mF))
  sigma2_re <- sum(vc$vcov)
  sigma2_res <- sigma(mF)^2
  sigma2_fe <- var(predict(mF, re.form = NA))
  r2m <- sigma2_fe / (sigma2_fe + sigma2_re + sigma2_res)
  r2c <- (sigma2_fe + sigma2_re) / (sigma2_fe + sigma2_re + sigma2_res)
}, error = function(e) {
  cat("[WARN] R2 calculation failed:", conditionMessage(e), "\n")
})

r2_stats <- data.frame(
  Model = "Full_Model_lmer_interaction_ver3",
  R2_Marginal = r2m,
  R2_Conditional = r2c,
  stringsAsFactors = FALSE
)

write.csv(fit_stats, out_file("1_pc_interaction_lmer_fit_stats"), row.names = FALSE)
write.csv(r2_stats, out_file("2_pc_interaction_lmer_r2_stats"), row.names = FALSE)

# -------------------------
# STEP 2) Coefficients + FDR
# -------------------------
fe <- fixef(mF)
V_unsc <- as.matrix(mF@pp$unsc())
se <- sigma(mF) * sqrt(diag(V_unsc))
tvals <- fe / se
pvals <- 2 * pnorm(abs(tvals), lower.tail = FALSE)

coef_df <- data.frame(
  Term = names(fe),
  Estimate = unname(fe),
  StdError = se,
  DF = NA,
  tValue = tvals,
  pValue = pvals,
  row.names = NULL,
  check.names = FALSE
)
cat("[INFO] Coefficients extracted directly (normal approx p-values, N=", nrow(dfm), ").\n")

coef_df$p_fdr_bh <- NA_real_
p_raw <- coef_df$pValue
idx <- which(coef_df$Term != "(Intercept)" & !is.na(p_raw))
if (length(idx) > 0) coef_df$p_fdr_bh[idx] <- p.adjust(p_raw[idx], method = "BH")

coef_df$sig_p_0.05 <- !is.na(p_raw) & (p_raw < 0.05)
coef_df$sig_fdr_0.05 <- !is.na(coef_df$p_fdr_bh) & (coef_df$p_fdr_bh < 0.05)

write.csv(coef_df, out_file("3_pc_interaction_lmer_coefficients_FDR"), row.names = FALSE)
cat("- Saved:", basename(out_file("3_pc_interaction_lmer_coefficients_FDR")), "\n")

# -------------------------
# STEP 3) LOSO Validation
# -------------------------
run_loso_cv_lmer <- function(formula_str, data, subject_col) {
  subjects <- unique(data[[subject_col]])
  n_subs <- length(subjects)

  res_sub <- character(0)
  res_rmse <- numeric(0)
  res_mae <- numeric(0)
  res_cor <- numeric(0)

  pb <- txtProgressBar(min = 0, max = n_subs, style = 3)
  for (i in seq_along(subjects)) {
    sid <- subjects[i]
    test_data <- data[data[[subject_col]] == sid, ]
    train_data <- data[data[[subject_col]] != sid, ]
    if (nrow(test_data) == 0) {
      setTxtProgressBar(pb, i)
      next
    }

    m_loso <- tryCatch({
      lme4::lmer(
        as.formula(formula_str),
        data = train_data,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))
      )
    }, error = function(e) NULL)

    if (!is.null(m_loso)) {
      preds <- tryCatch({
        predict(m_loso, newdata = test_data, allow.new.levels = TRUE)
      }, error = function(e) NULL)

      if (!is.null(preds)) {
        actual <- test_data$log_rt
        rmse <- sqrt(mean((preds - actual)^2))
        mae <- mean(abs(preds - actual))
        corv <- suppressWarnings(cor(preds, actual))

        res_sub <- c(res_sub, as.character(sid))
        res_rmse <- c(res_rmse, rmse)
        res_mae <- c(res_mae, mae)
        res_cor <- c(res_cor, corv)
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)

  data.frame(
    Subject_ID = res_sub,
    RMSE = res_rmse,
    MAE = res_mae,
    Corr = res_cor,
    stringsAsFactors = FALSE
  )
}

cat("\n[INFO] Running LOSO CV...\n")
loso_detail <- run_loso_cv_lmer(fF_str, dfm, "subject_id")
write.csv(loso_detail, out_file("4_pc_interaction_lmer_loso_by_subject"), row.names = FALSE)

loso_summary <- loso_detail %>%
  summarise(
    RMSE_mean = mean(RMSE, na.rm = TRUE),
    MAE_mean = mean(MAE, na.rm = TRUE),
    Corr_mean = mean(Corr, na.rm = TRUE)
  ) %>%
  as.data.frame()

write.csv(loso_summary, out_file("5_pc_interaction_lmer_loso_summary"), row.names = FALSE)

# # -------------------------
# # STEP 4) Subject Bootstrap
# # -------------------------
# run_subject_bootstrap_lmer <- function(formula_str, data, subject_col, B = 300, seed = 42) {
#   set.seed(seed)
#   subjects <- unique(data[[subject_col]])
#   n_sub <- length(subjects)

#   coef_list <- list()
#   pb <- txtProgressBar(min = 0, max = B, style = 3)

#   for (b in 1:B) {
#     boot_sub <- sample(subjects, size = n_sub, replace = TRUE)
#     boot_data <- bind_rows(lapply(boot_sub, function(s) data[data[[subject_col]] == s, ]))
#     boot_data <- droplevels(boot_data)

#     m_boot <- tryCatch({
#       lme4::lmer(
#         as.formula(formula_str),
#         data = boot_data,
#         REML = FALSE,
#         control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
#       )
#     }, error = function(e) NULL)

#     if (!is.null(m_boot)) {
#       fe_b <- fixef(m_boot)
#       V_b <- as.matrix(m_boot@pp$unsc())
#       se_b <- sigma(m_boot) * sqrt(diag(V_b))
#       t_b <- fe_b / se_b
#       p_vals <- 2 * pnorm(abs(t_b), lower.tail = FALSE)

#       tmp <- data.frame(
#         iter = b,
#         Term = names(fe_b),
#         Estimate = unname(fe_b),
#         p = p_vals,
#         row.names = NULL
#       )

#       p_fdr_iter <- rep(NA_real_, nrow(tmp))
#       idx2 <- which(tmp$Term != "(Intercept)" & !is.na(tmp$p))
#       if (length(idx2) > 0) p_fdr_iter[idx2] <- p.adjust(tmp$p[idx2], method = "BH")

#       tmp$p_fdr_bh <- p_fdr_iter
#       tmp$sig_p_0.05 <- !is.na(tmp$p) & (tmp$p < 0.05)
#       tmp$sig_fdr_0.05 <- !is.na(tmp$p_fdr_bh) & (tmp$p_fdr_bh < 0.05)
#       coef_list[[length(coef_list) + 1]] <- tmp
#     }

#     setTxtProgressBar(pb, b)
#   }
#   close(pb)

#   if (length(coef_list) == 0) return(data.frame())
#   bind_rows(coef_list)
# }

# cat("\n[INFO] Running Bootstrap...\n")
# boot_B <- 300
# boot_coefs <- run_subject_bootstrap_lmer(fF_str, dfm, "subject_id", B = boot_B, seed = 42)
# write.csv(boot_coefs, out_file("6_pc_interaction_lmer_bootstrap_coefficients_long"), row.names = FALSE)

# if (nrow(boot_coefs) > 0) {
#   boot_summary <- boot_coefs %>%
#     group_by(Term) %>%
#     summarise(
#       mean_Est = mean(Estimate, na.rm = TRUE),
#       ci_low = quantile(Estimate, 0.025, na.rm = TRUE),
#       ci_high = quantile(Estimate, 0.975, na.rm = TRUE),
#       prop_sig_p_0.05 = mean(sig_p_0.05, na.rm = TRUE),
#       prop_sig_fdr_0.05 = mean(sig_fdr_0.05, na.rm = TRUE)
#     ) %>%
#     as.data.frame()
# } else {
#   boot_summary <- data.frame()
# }

# write.csv(boot_summary, out_file("9_pc_interaction_lmer_bootstrap_summary"), row.names = FALSE)
# cat("- Saved:", basename(out_file("9_pc_interaction_lmer_bootstrap_summary")), "\n")

# cat("\n=== ALL DONE ===\n")

# ---------------------------------------------------------
# STEP 5) Interaction Plot (Marginal Effects)
# ---------------------------------------------------------
cat("\n[INFO] Generating Marginal Effects Interaction Plots...\n")
# install.packages("patchwork")

# 필요한 패키지 로드 (없다면 install.packages("ggeffects"), install.packages("patchwork") 실행 필요)
suppressPackageStartupMessages({
  library(ggeffects)
  library(ggplot2)
  library(patchwork)
})

# 시각화를 위한 함수 정의
# ggpredict를 사용해 LME 모델의 수식(beta)을 바탕으로 예측된 반응 시간을 계산합니다.
plot_interaction_marginal <- function(model, pc_var, pc_label, panel_title) {
  
  # 1. ggpredict로 예측값 생성
  pred_terms <- c(paste0(pc_var, " [all]"), "trial_ear_norm_baseline_sd [-1sd, 1sd]")
  pred_df <- ggpredict(model, terms = pred_terms) %>% as.data.frame()
  
  # 2. [해결 핵심] 크기를 비교하여 안전하게 그룹을 분리합니다.
  pred_df <- pred_df %>%
    mutate(
      # group 변수를 실제 숫자로 변환
      group_num = as.numeric(as.character(group)),
      
      # 작은 값(min)은 Alert(-1 SD), 큰 값은 Unstable(+1 SD)로 명확히 지정
      group_label = ifelse(group_num == min(group_num), 
                           "Alert (-1 SD of EAR SD)", 
                           "Unstable (+1 SD of EAR SD)")
    )
  
  # 3. 반응시간(RT)을 초(second) 단위로 되돌리기
  pred_df <- pred_df %>%
    mutate(
      predicted_rt = exp(predicted),
      conf.low_rt = exp(conf.low),
      conf.high_rt = exp(conf.high)
    )
  
  # 4. 그래프 그리기 (group = group_label 추가로 선 꼬임 방지)
  p <- ggplot(pred_df, aes(x = x, y = predicted_rt, 
                           color = group_label, 
                           fill = group_label,
                           group = group_label)) +  
    geom_ribbon(aes(ymin = conf.low_rt, ymax = conf.high_rt), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("Alert (-1 SD of EAR SD)" = "#457B9D", "Unstable (+1 SD of EAR SD)" = "#C62828")) +
    scale_fill_manual(values = c("Alert (-1 SD of EAR SD)" = "#457B9D", "Unstable (+1 SD of EAR SD)" = "#C62828")) +
    labs(
      title = panel_title,
      x = pc_label,
      y = "Predicted Reaction Time (s)",
      color = "Driver State",
      fill = "Driver State"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 15),
      axis.title = element_text(face = "bold")
    )
  
  return(p)
}

# 1. Panel A: PC2 (Negative/Tedious) x EAR_SD
p_pc2 <- plot_interaction_marginal(mF, "snd_pc2_z", "Sound PC2 z-score", "(A) Sound PC2 × EAR SD Interaction")
# 2. Panel B: PC4 (Coldness) x EAR_SD
# 코드 상에서 PC 변수들이 snd_pc4_z 처럼 z-스케일링 되어 있다고 가정합니다.
p_pc4 <- plot_interaction_marginal(mF, "snd_pc4_z", "Sound PC4 z-score", "(B) Sound PC4 × EAR SD Interaction")



# 3. 두 패널을 나란히 병합 (Patchwork 패키지)
# guides = 'collect'를 통해 범례를 하나로 합쳐 하단에 배치합니다.
final_interaction_plot <- p_pc2 + p_pc4 + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

# 4. 고해상도 이미지 및 PDF 저장
plot_path_png <- file.path(output_dir, paste0("10_interaction_marginal_effects", out_suffix, ".png"))
plot_path_pdf <- file.path(output_dir, paste0("10_interaction_marginal_effects", out_suffix, ".pdf"))

ggsave(plot_path_png, final_interaction_plot, width = 14, height = 6, dpi = 600)
ggsave(plot_path_pdf, final_interaction_plot, width = 14, height = 6)

cat("- Saved Marginal Effects plots to:", plot_path_png, "\n")