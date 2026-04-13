#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(MuMIn)
  library(car)      # VIF 계산용
  library(lmerTest)
})

# -------------------------
# Paths
# -------------------------
in_csv <- "/users/kcl/code/ES/dfm_for_nlme_fixed900ms.csv"
window_src_csv <- "/users/kcl/data/ES/merged_6.csv"
false_count_src_csv <- "/users/kcl/code/ES/dfm_for_lme4_final_with_drowsy_fixed900ms.csv"
output_dir <- getwd()
output_prefix <- "pc_lmer_wo_demo_ver3"

cat("Results will be saved to:", output_dir, "\n")
cat("Window source CSV:", window_src_csv, "\n")
cat("False-count source CSV:", false_count_src_csv, "\n")

parse_hms_to_sec <- function(x) {
  if (is.na(x)) return(NA_real_)
  if (is.numeric(x)) return(as.numeric(x))
  s <- trimws(as.character(x))
  if (!nzchar(s)) return(NA_real_)
  parts <- strsplit(s, ":", fixed = TRUE)[[1]]
  if (length(parts) != 3) return(NA_real_)
  h <- suppressWarnings(as.numeric(parts[1]))
  m <- suppressWarnings(as.numeric(parts[2]))
  sec <- suppressWarnings(as.numeric(parts[3]))
  if (!all(is.finite(c(h, m, sec)))) return(NA_real_)
  h * 3600 + m * 60 + sec
}

wallclock_diff_sec <- function(later_sec, earlier_sec) {
  if (!is.finite(later_sec) || !is.finite(earlier_sec)) return(NA_real_)
  diff <- later_sec - earlier_sec
  if (diff > 43200) {
    diff <- diff - 86400
  } else if (diff < -43200) {
    diff <- diff + 86400
  }
  diff
}

build_trial_window_map <- function(csv_path) {
  if (!file.exists(csv_path)) {
    cat("[WARN] trial window source not found:", csv_path, "\n")
    return(NULL)
  }

  w <- read.csv(csv_path, stringsAsFactors = FALSE)
  req <- c("subject_id", "trial", "cue_time", "sprite_switch_time", "response_time")
  miss <- setdiff(req, names(w))
  if (length(miss) > 0) {
    cat("[WARN] trial window source missing columns:", paste(miss, collapse = ", "), "\n")
    return(NULL)
  }

  w$subject_id_chr <- as.character(w$subject_id)
  w$trial_num <- suppressWarnings(as.integer(as.character(w$trial)))
  cue_sec <- vapply(w$cue_time, parse_hms_to_sec, numeric(1))
  sprite_sec <- vapply(w$sprite_switch_time, parse_hms_to_sec, numeric(1))
  rt_sec <- suppressWarnings(as.numeric(as.character(w$response_time)))
  stim_dur_sec <- mapply(wallclock_diff_sec, sprite_sec, cue_sec)
  w$trial_window_sec <- stim_dur_sec + rt_sec

  out <- w %>%
    select(subject_id_chr, trial_num, trial_window_sec) %>%
    filter(!is.na(subject_id_chr), !is.na(trial_num)) %>%
    group_by(subject_id_chr, trial_num) %>%
    summarise(
      trial_window_sec = if (all(!is.finite(trial_window_sec))) NA_real_
      else median(trial_window_sec[is.finite(trial_window_sec)]),
      .groups = "drop"
    )

  as.data.frame(out)
}

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

# [Fix] log_rt 생성
if (!"log_rt" %in% names(dfm)) {
  if ("response_time" %in% names(dfm)) {
    cat("[INFO] 'log_rt' not found. Creating from 'response_time'...\n")
    dfm$log_rt <- ifelse(dfm$response_time > 0, log(dfm$response_time), NA)
  } else {
    stop("[ERROR] Neither 'log_rt' nor 'response_time' columns found.")
  }
}

# Join keys for covariates
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

# trial window covariate from merged_6.csv
window_map <- build_trial_window_map(window_src_csv)
if (!is.null(window_map)) {
  if ("trial_window_sec" %in% names(dfm)) {
    dfm$trial_window_sec <- suppressWarnings(as.numeric(as.character(dfm$trial_window_sec)))
  }
  dfm <- dfm %>%
    left_join(window_map, by = c("subject_id_chr", "trial_num"), suffix = c("", "_from_src"))

  if ("trial_window_sec_from_src" %in% names(dfm)) {
    if ("trial_window_sec" %in% names(dfm)) {
      idx_fill <- (!is.finite(dfm$trial_window_sec)) & is.finite(dfm$trial_window_sec_from_src)
      dfm$trial_window_sec[idx_fill] <- dfm$trial_window_sec_from_src[idx_fill]
      dfm$trial_window_sec_from_src <- NULL
    } else {
      names(dfm)[names(dfm) == "trial_window_sec_from_src"] <- "trial_window_sec"
    }
  }
  cat("[INFO] trial_window_sec ready (non-NA):", sum(is.finite(dfm$trial_window_sec)), "/", nrow(dfm), "\n")
} else {
  cat("[WARN] trial_window_sec could not be created from source CSV.\n")
}

# -------------------------
# 2) Type casting & Scaling
# -------------------------
# Factor 변환
dfm$subject_id <- factor(dfm$subject_id)

# [중요] BGM ID 변수 설정 (sound_norm을 BGM ID로 가정)
if ("sound_norm" %in% names(dfm)) {
  dfm$sound_norm <- factor(dfm$sound_norm)
  bgm_id_col <- "sound_norm" 
} else {
  # 만약 sound_norm이 없다면, is_bgm dummies에서 역추적하거나 에러 처리
  # 여기서는 sound_norm이 있다고 가정합니다.
  cat("[WARN] 'sound_norm' column not found! Assuming 'is_bgm' implies structure.\n")
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

fac_cols <- c("cueType","arrowDirection","congruency","targetDirection","drive_exp")
for (cc in fac_cols) if (cc %in% names(dfm)) dfm[[cc]] <- factor(dfm[[cc]])

# PC scores (PC1, PC2 등)
snd_pc_cols <- grep("^snd_pc[0-9]+$", names(dfm), value = TRUE)

# Psychoacoustic vars
psycho_vars <- c("Sharpness", "Roughness", "Tonality")

# [추가] 심리 변수 및 기타 연속형 변수 정의
mental_vars <- c("BDI_II", "BAI", "PSS")

# Continuous vars to scale
# avg_ear_norm is already normalized upstream (subject-wise z-score),
# so we keep it as-is and do not re-scale here.
# drowsy_cont_vars <- intersect(c("eye_status_ratio"), names(dfm))
# [수정] cont_to_scale에 mental_vars 추가
cont_to_scale <- c("age", "sex", "satisfaction_filled", 
                   snd_pc_cols, psycho_vars, mental_vars,
                   "block_order")

# # [추가] eye_status를 Factor로 변환 (0/1인 경우)
# if ("eye_status" %in% names(dfm)) {
#   dfm$eye_status <- factor(as.integer(dfm$eye_status))
# }
cont_to_scale <- intersect(cont_to_scale, names(dfm))

# Scaling (Standardization)
for (cc in cont_to_scale) {
  dfm[[paste0(cc, "_z")]] <- as.numeric(scale(dfm[[cc]]))
}

# NA Removal
vars_needed <- c("log_rt", "subject_id", bgm_id_col, "cueType", "arrowDirection", "congruency", "targetDirection")
# if ("eye_status" %in% names(dfm)) vars_needed <- c(vars_needed, "eye_status")

mental_vars_z <- paste0(mental_vars, "_z")

# (PC, EAR 등이 NA인 행도 제거해야 함)
check_vars_exist <- c(
  paste0(snd_pc_cols, "_z"),
  # paste0(drowsy_cont_vars, "_z"),
  "trial_ear_norm_baseline",
  "trial_ear_norm_baseline_sd",
  "satisfaction_filled_z",
  mental_vars_z,
  "block_order_z"
)
vars_needed <- c(vars_needed, intersect(check_vars_exist, names(dfm)))

dfm <- dfm[complete.cases(dfm[, vars_needed]), ]
dfm <- droplevels(dfm)
cat("After final NA removal:", nrow(dfm), "rows\n")


# -------------------------
# 3) Formula Construction (lme4 style)
# -------------------------
# Task
task_part <- "cueType + arrowDirection + congruency + targetDirection"

# Demo
demo_terms <- c()
if ("age_z" %in% names(dfm))   demo_terms <- c(demo_terms, "age_z")
if ("sex_z" %in% names(dfm))   demo_terms <- c(demo_terms, "sex_z")
if ("drive_exp" %in% names(dfm)) demo_terms <- c(demo_terms, "drive_exp")

for (mv in c("BDI_II_z", "BAI_z", "PSS_z")) {
  if (mv %in% names(dfm)) demo_terms <- c(demo_terms, mv)
}
demo_part <- paste(demo_terms, collapse = " + ")

# Satisfaction
satis_part <- if ("satisfaction_filled_z" %in% names(dfm)) " + satisfaction_filled_z" else ""

# PC (Emotion)
snd_part <- if (length(snd_pc_cols) > 0) paste(paste0(snd_pc_cols, "_z"), collapse = " + ") else ""

# Drowsy (EAR or Label)
drowsy_part <- ""
if ("trial_ear_norm_baseline" %in% names(dfm)) drowsy_part <- paste0(drowsy_part, " + trial_ear_norm_baseline")
if ("trial_ear_norm_baseline_sd" %in% names(dfm)) drowsy_part <- paste0(drowsy_part, " + trial_ear_norm_baseline_sd")

# Additional confounds
confound_terms <- c()
if ("block_order_z" %in% names(dfm)) {
  confound_terms <- c(confound_terms, "block_order_z")
} else if ("block_order" %in% names(dfm)) {
  confound_terms <- c(confound_terms, "block_order")
}
confound_part <- paste(confound_terms, collapse = " + ")

# if ("eye_status" %in% names(dfm)) drowsy_part <- paste0(drowsy_part, " + eye_status")
# 필요시 label도 추가 가능: if ("drowsiness_label" %in% names(dfm)) drowsy_part <- paste0(drowsy_part, " + drowsiness_label")

# [핵심] Random Effects Part: (1|subject) + (1|BGM)
# sound_norm이 BGM ID라고 가정
# random_part <- "+ (1 | subject_id)"
# if (!is.null(bgm_id_col)) {
#   # random_part <- paste0(random_part, " + (1 | ", bgm_id_col, ") + (1 | subject_id:", bgm_id_col, ")")
#   random_part <- paste0("+ (1 | subject_id) + (1 | subject_id:", bgm_id_col, ")")
# }
if (!is.null(bgm_id_col) && (bgm_id_col %in% names(dfm))) {
  random_part <- paste0(
    "+ (1 | subject_id) ",
    "+ (1 | subject_id:", bgm_id_col, ")"
  )
} else {
  random_part <- "+ (1 | subject_id)"
}

# Fixed Effects 합치기
fixed_terms <- c(task_part)

# # BGM ID를 고정효과(독립변수)로 추가
# if (!is.null(bgm_id_col) && bgm_id_col %in% names(dfm)) {
#   fixed_terms <- c(fixed_terms, bgm_id_col)
# }

# if (nchar(demo_part) > 0) fixed_terms <- c(fixed_terms, demo_part)
if (nchar(satis_part) > 0) fixed_terms <- c(fixed_terms, gsub("^ \\+ ", "", satis_part))
if (nchar(snd_part) > 0) fixed_terms <- c(fixed_terms, snd_part)
if (nchar(drowsy_part) > 0) fixed_terms <- c(fixed_terms, gsub("^ \\+ ", "", drowsy_part))
if (nchar(confound_part) > 0) fixed_terms <- c(fixed_terms, confound_part)

fixed_terms <- unique(fixed_terms)
fixed_formula <- paste(fixed_terms, collapse = " + ")

# 최종 수식 문자열
fF_str <- paste("log_rt ~", fixed_formula, random_part)

cat("\n=== FULL MODEL FORMULA (lmer) ===\n", fF_str, "\n")


# -------------------------
# 4) Fit Model (lmer)
# -------------------------
# REML = FALSE (ML) is generally better for comparing models with different fixed effects via AIC/BIC
# But for reporting coefficients, REML = TRUE is often standard. 
# Here we use REML = FALSE to match your previous 'method="ML"' setting.
# 3. 모델 피팅: lme4::lmer 사용 (lmerTest::lmer 아님)
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

cat("\n[OK] Full model fitted using lmer (p-values via normal approx).\n")

# -------------------------
# [VIF Check] — computed from lm (vcov on merMod segfaults in Matrix C code)
# -------------------------
cat("\n[INFO] Calculating VIF from fixed-effects lm proxy...\n")
fixed_only_str <- paste("log_rt ~", fixed_formula) 
mF_lm <- lm(as.formula(fixed_only_str), data = dfm)
vif_res <- tryCatch({ car::vif(mF_lm) }, error = function(e) NULL)

if (!is.null(vif_res)) {
  if (is.matrix(vif_res)) {
    vif_df <- data.frame(Term = rownames(vif_res), vif_res)
  } else {
    vif_df <- data.frame(Term = names(vif_res), VIF = vif_res)
  }
  write.csv(vif_df, file.path(output_dir, paste0("0_", output_prefix, "_model_VIF.csv")), row.names = FALSE)
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
plot_path <- file.path(output_dir, paste0("0_", output_prefix, "_diagnostics_v1.png"))
ggsave(filename = plot_path, plot = p_final, width = 12, height = 6, dpi = 300, bg = "white")

cat("- Saved diagnostic plot:", plot_path, "\n")
# ---------------------------------------------------------

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
plot_path <- file.path(output_dir, paste0("0_", output_prefix, "_diagnostics.png"))
ggsave(filename = plot_path, plot = p_final, width = 12, height = 6, dpi = 300, bg = "white")

cat("- Saved diagnostic plot:", plot_path, "\n")
# ---------------------------------------------------------

# -------------------------
# STEP 1) Fit Stats
# -------------------------
fit_stats <- data.frame(
  Model = "Full_Model_lmer",
  AIC = AIC(mF),
  BIC = BIC(mF),
  logLik = as.numeric(logLik(mF)),
  stringsAsFactors = FALSE
)

# Manual R2 calculation (Nakagawa & Schielzeth) to avoid vcov() segfault
r2m <- NA_real_
r2c <- NA_real_
tryCatch({
  vc <- as.data.frame(VarCorr(mF))
  sigma2_re <- sum(vc$vcov)              # total random effect variance
  sigma2_res <- sigma(mF)^2              # residual variance
  # fixed-effect variance: var of fitted values from fixed effects only
  sigma2_fe <- var(predict(mF, re.form = NA))
  r2m <- sigma2_fe / (sigma2_fe + sigma2_re + sigma2_res)
  r2c <- (sigma2_fe + sigma2_re) / (sigma2_fe + sigma2_re + sigma2_res)
}, error = function(e) {
  cat("[WARN] R2 calculation failed:", conditionMessage(e), "\n")
})
r2_stats <- data.frame(
  Model = "Full_Model_lmer",
  R2_Marginal = r2m,
  R2_Conditional = r2c,
  stringsAsFactors = FALSE
)

write.csv(fit_stats, file.path(output_dir, paste0("1_", output_prefix, "_fit_stats.csv")), row.names = FALSE)
write.csv(r2_stats,  file.path(output_dir, paste0("2_", output_prefix, "_r2_stats.csv")),  row.names = FALSE)


# -------------------------
# STEP 2) Coefficients + FDR
# -------------------------
# NOTE: summary() and vcov() segfault on this model (Matrix C bug).
# Extract coefficients directly from model slots instead.
fe <- fixef(mF)
V_unsc <- as.matrix(mF@pp$unsc())   # unscaled vcov, bypass dpoMatrix conversion
se <- sigma(mF) * sqrt(diag(V_unsc))
tvals <- fe / se
pvals <- 2 * pnorm(abs(tvals), lower.tail = FALSE)  # normal approx (N=17262)

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

# FDR Calculation
coef_df$p_fdr_bh <- NA_real_
p_raw <- coef_df$pValue
idx <- which(coef_df$Term != "(Intercept)" & !is.na(p_raw))

if (length(idx) > 0) {
  coef_df$p_fdr_bh[idx] <- p.adjust(p_raw[idx], method = "BH")
}

coef_df$sig_p_0.05   <- !is.na(p_raw) & (p_raw < 0.05)
coef_df$sig_fdr_0.05 <- !is.na(coef_df$p_fdr_bh) & (coef_df$p_fdr_bh < 0.05)

coef_out_name <- paste0("3_", output_prefix, "_coefficients_FDR.csv")
write.csv(coef_df, file.path(output_dir, coef_out_name), row.names = FALSE)
cat("- Saved:", coef_out_name, "\n")


# -------------------------
# STEP 3) LOSO Validation (Updated for lmer)
# -------------------------
run_loso_cv_lmer <- function(formula_str, data, subject_col) {
  subjects <- unique(data[[subject_col]])
  n_subs <- length(subjects)
  
  res_sub <- character(0); res_rmse <- numeric(0); res_mae  <- numeric(0); res_cor  <- numeric(0)
  
  pb <- txtProgressBar(min = 0, max = n_subs, style = 3)
  
  for (i in seq_along(subjects)) {
    sid <- subjects[i]
    test_data  <- data[data[[subject_col]] == sid, ]
    train_data <- data[data[[subject_col]] != sid, ]
    
    # Check levels (Skip if levels mismatch too much, or try catch)
    if(nrow(test_data) == 0) { setTxtProgressBar(pb, i); next }
    
    m_loso <- tryCatch({
      lme4::lmer(as.formula(formula_str), data = train_data, REML = FALSE,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=10000)))
    }, error = function(e) NULL)
    
    if (!is.null(m_loso)) {
      # allow.new.levels = TRUE 필수 (Test set에 새로운 subject가 있으므로... 
      # 사실 subject는 random effect이므로 0으로 수렴하거나 population mean 사용)
      preds <- tryCatch({
         predict(m_loso, newdata = test_data, allow.new.levels = TRUE)
      }, error = function(e) NULL)
      
      if(!is.null(preds)){
          actual <- test_data$log_rt
          rmse <- sqrt(mean((preds - actual)^2))
          mae  <- mean(abs(preds - actual))
          corv <- suppressWarnings(cor(preds, actual))
          
          res_sub  <- c(res_sub, as.character(sid))
          res_rmse <- c(res_rmse, rmse)
          res_mae  <- c(res_mae, mae)
          res_cor  <- c(res_cor, corv)
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  data.frame(Subject_ID = res_sub, RMSE = res_rmse, MAE = res_mae, Corr = res_cor, stringsAsFactors = FALSE)
}

cat("\n[INFO] Running LOSO CV...\n")
loso_detail <- run_loso_cv_lmer(fF_str, dfm, "subject_id")
write.csv(loso_detail, file.path(output_dir, paste0("4_", output_prefix, "_loso_by_subject.csv")), row.names = FALSE)

loso_summary <- loso_detail %>%
  summarise(
    RMSE_mean = mean(RMSE, na.rm=TRUE),
    MAE_mean = mean(MAE, na.rm=TRUE),
    Corr_mean = mean(Corr, na.rm=TRUE)
  ) %>% as.data.frame()
write.csv(loso_summary, file.path(output_dir, paste0("5_", output_prefix, "_loso_summary.csv")), row.names = FALSE)


# # -------------------------
# # STEP 4) Subject Bootstrap (Updated for lmer)
# # -------------------------
# run_subject_bootstrap_lmer <- function(formula_str, data, subject_col, B=300, seed=42) {
#   set.seed(seed)
#   subjects <- unique(data[[subject_col]])
#   n_sub <- length(subjects)
  
#   coef_list <- list()
  
#   pb <- txtProgressBar(min=0, max=B, style=3)
  
#   for (b in 1:B) {
#     boot_sub <- sample(subjects, size=n_sub, replace=TRUE)
#     # 데이터 복원 추출
#     boot_data <- bind_rows(lapply(boot_sub, function(s) data[data[[subject_col]] == s, ]))
#     boot_data <- droplevels(boot_data)
    
#     m_boot <- tryCatch({
#       lme4::lmer(as.formula(formula_str), data = boot_data, REML = FALSE,
#                  control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)) # 속도 최적화
#     }, error = function(e) NULL)
    
#     if (!is.null(m_boot)) {
#       # Extract coefficients directly (avoid summary/vcov segfault)
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
      
#       # FDR Calculation per iter
#       p_fdr_iter <- rep(NA_real_, nrow(tmp))
#       idx2 <- which(tmp$Term != "(Intercept)" & !is.na(tmp$p))
#       if(length(idx2)>0) p_fdr_iter[idx2] <- p.adjust(tmp$p[idx2], method = "BH")
      
#       tmp$p_fdr_bh <- p_fdr_iter
#       tmp$sig_p_0.05 <- !is.na(tmp$p) & (tmp$p < 0.05)
#       tmp$sig_fdr_0.05 <- !is.na(tmp$p_fdr_bh) & (tmp$p_fdr_bh < 0.05)
      
#       coef_list[[length(coef_list)+1]] <- tmp
#     }
#     setTxtProgressBar(pb, b)
#   }
#   close(pb)
  
#   bind_rows(coef_list)
# }

# cat("\n[INFO] Running Bootstrap...\n")
# boot_B <- 300
# boot_coefs <- run_subject_bootstrap_lmer(fF_str, dfm, "subject_id", B=boot_B, seed=42)

# write.csv(boot_coefs, file.path(output_dir, "6_pc_lmer_bootstrap_coefficients_long.csv"), row.names=FALSE)

# # Bootstrap Summary
# boot_summary <- boot_coefs %>%
#   group_by(Term) %>%
#   summarise(
#     mean_Est = mean(Estimate, na.rm=TRUE),
#     ci_low = quantile(Estimate, 0.025, na.rm=TRUE),
#     ci_high= quantile(Estimate, 0.975, na.rm=TRUE),
#     prop_sig_p_0.05   = mean(sig_p_0.05, na.rm=TRUE),
#     prop_sig_fdr_0.05 = mean(sig_fdr_0.05, na.rm=TRUE)
#   ) %>% as.data.frame()

# write.csv(boot_summary, file.path(output_dir, "9_pc_lmer_bootstrap_summary.csv"), row.names=FALSE)
# cat("- Saved: 9_pc_lmer_bootstrap_summary.csv\n")

# cat("\n=== ALL DONE ===\n")
