library(tidyverse)

edge_file <- list.files(
  ".",
  pattern = "spat21_aim2_merged_data_with_weights_5MAR2020\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(edge_file) != 1) {
  stop("Could not uniquely identify the weighted P(TE) edge file.")
}

first_nonmissing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA else x[1]
}

edges <- read_csv(edge_file, show_col_types = FALSE)

required <- c(
  "sample_id_human",
  "p_te_all_csp",
  "p_te_all_ama",
  "pfr364Q_std_combined",
  "aim2_exposure",
  "mean_moi",
  "total_num_mosq_in_hh",
  "age_cat_baseline",
  "HH_ID_human",
  "unq_memID"
)

missing_vars <- setdiff(required, names(edges))

if (length(missing_vars) > 0) {
  stop("Missing required columns: ", paste(missing_vars, collapse = ", "))
}

edges <- edges %>%
  mutate(
    p_te_all_csp = as.numeric(p_te_all_csp),
    p_te_all_ama = as.numeric(p_te_all_ama)
  )

coverage <- tibble(
  edges_total = nrow(edges),
  csp_observed = sum(!is.na(edges$p_te_all_csp)),
  ama_observed = sum(!is.na(edges$p_te_all_ama)),
  both_observed = sum(!is.na(edges$p_te_all_csp) & !is.na(edges$p_te_all_ama))
)

cat("\nMarker coverage:\n")
print(coverage, width = Inf)

# Restrict to the same human-mosquito pairs so marker choice is the only difference.
common_edges <- edges %>%
  filter(!is.na(p_te_all_csp), !is.na(p_te_all_ama))

infection_data <- common_edges %>%
  group_by(sample_id_human) %>%
  summarise(
    n_edges = n(),
    sum_csp = sum(p_te_all_csp),
    sum_ama = sum(p_te_all_ama),
    avg_csp = mean(p_te_all_csp),
    avg_ama = mean(p_te_all_ama),
    density = first_nonmissing(pfr364Q_std_combined),
    exposure = first_nonmissing(aim2_exposure),
    moi = first_nonmissing(mean_moi),
    hh_mosq = first_nonmissing(total_num_mosq_in_hh),
    age = first_nonmissing(age_cat_baseline),
    hh = first_nonmissing(HH_ID_human),
    memid = first_nonmissing(unq_memID),
    .groups = "drop"
  ) %>%
  mutate(
    density = as.numeric(density),
    moi = as.numeric(moi),
    hh_mosq = as.numeric(hh_mosq),
    ldens = log1p(density),
    lmosq = log1p(hh_mosq),
    asym = as.integer(exposure == "asymptomatic infection"),
    agecat = factor(age)
  ) %>%
  filter(
    is.finite(ldens),
    is.finite(lmosq),
    is.finite(moi),
    !is.na(agecat),
    !is.na(hh),
    !is.na(memid)
  )

fit_marker <- function(data, marker, min_edges) {
  outcome <- paste0("avg_", marker)

  d <- data %>%
    filter(n_edges >= min_edges)

  fit <- lm(
    reformulate(
      c("ldens", "asym", "moi", "agecat", "lmosq"),
      response = outcome
    ),
    data = d
  )

  d <- d %>%
    mutate(
      predicted = predict(fit),
      residual = .data[[outcome]] - predicted,
      surprising = residual >= quantile(residual, 0.90, na.rm = TRUE),
      marker = toupper(marker),
      min_edges = min_edges
    )

  coefs <- summary(fit)$coefficients

  summary_row <- tibble(
    marker = toupper(marker),
    min_edges = min_edges,
    infections = nrow(d),
    people = n_distinct(d$memid),
    households = n_distinct(d$hh),
    corr_avg_pte_edges = cor(d[[outcome]], d$n_edges),
    r_squared = summary(fit)$r.squared,
    surprising = sum(d$surprising),
    moi_estimate = coefs["moi", "Estimate"],
    moi_p = coefs["moi", "Pr(>|t|)"],
    mosquito_estimate = coefs["lmosq", "Estimate"],
    mosquito_p = coefs["lmosq", "Pr(>|t|)"]
  )

  list(
    data = d,
    fit = fit,
    summary = summary_row
  )
}

cutoffs <- c(1, 5, 10)

fits <- list()

for (cutoff in cutoffs) {
  fits[[paste0("csp_", cutoff)]] <- fit_marker(
    infection_data,
    "csp",
    cutoff
  )

  fits[[paste0("ama_", cutoff)]] <- fit_marker(
    infection_data,
    "ama",
    cutoff
  )
}

model_summary <- bind_rows(
  lapply(fits, `[[`, "summary")
)

compare_markers <- function(cutoff) {
  csp <- fits[[paste0("csp_", cutoff)]]$data
  ama <- fits[[paste0("ama_", cutoff)]]$data

  joined <- csp %>%
    select(
      sample_id_human,
      avg_csp,
      residual_csp = residual,
      surprising_csp = surprising
    ) %>%
    inner_join(
      ama %>%
        select(
          sample_id_human,
          avg_ama,
          residual_ama = residual,
          surprising_ama = surprising
        ),
      by = "sample_id_human"
    )

  csp_ids <- joined %>%
    filter(surprising_csp) %>%
    pull(sample_id_human)

  ama_ids <- joined %>%
    filter(surprising_ama) %>%
    pull(sample_id_human)

  overlap <- length(intersect(csp_ids, ama_ids))
  union_n <- length(union(csp_ids, ama_ids))

  tibble(
    min_edges = cutoff,
    infections = nrow(joined),
    corr_avg_pte_csp_ama = cor(joined$avg_csp, joined$avg_ama),
    corr_residual_csp_ama = cor(joined$residual_csp, joined$residual_ama),
    surprising_csp = length(csp_ids),
    surprising_ama = length(ama_ids),
    surprise_overlap = overlap,
    jaccard = overlap / union_n
  )
}

marker_comparison <- bind_rows(
  lapply(cutoffs, compare_markers)
)

csp5 <- fits[["csp_5"]]$data %>%
  filter(surprising) %>%
  pull(sample_id_human)

ama5 <- fits[["ama_5"]]$data %>%
  filter(surprising) %>%
  pull(sample_id_human)

csp10 <- fits[["csp_10"]]$data %>%
  filter(surprising) %>%
  pull(sample_id_human)

ama10 <- fits[["ama_10"]]$data %>%
  filter(surprising) %>%
  pull(sample_id_human)

robust_ids <- Reduce(
  intersect,
  list(csp5, ama5, csp10, ama10)
)

robust_core <- infection_data %>%
  filter(sample_id_human %in% robust_ids) %>%
  arrange(sample_id_human)

cat("\nModel summaries:\n")
print(model_summary, n = Inf, width = Inf)

cat("\nCSP vs AMA stability:\n")
print(marker_comparison, n = Inf, width = Inf)

cat(
  "\nSurprising in CSP and AMA at both >=5 and >=10 links:",
  length(robust_ids),
  "\n"
)

dir.create(
  file.path("data", "processed"),
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  model_summary,
  file.path(
    "data",
    "processed",
    "pte-csp-ama-model-summary.csv"
  )
)

write_csv(
  marker_comparison,
  file.path(
    "data",
    "processed",
    "pte-csp-ama-stability.csv"
  )
)

write_csv(
  robust_core,
  file.path(
    "data",
    "processed",
    "pte-csp-ama-robust-core.csv"
  )
)

for (name in names(fits)) {
  write_csv(
    fits[[name]]$data %>%
      filter(surprising) %>%
      arrange(desc(residual)),
    file.path(
      "data",
      "processed",
      paste0("pte-surprising-", name, ".csv")
    )
  )
}
