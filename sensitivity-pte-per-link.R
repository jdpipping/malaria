library(tidyverse)

input_file <- file.path("data", "processed", "pte-infection-level-model-sample.csv")

if (!file.exists(input_file)) {
  input_file <- "C:/Users/Alayna Veeck/Documents/SUMR2/malaria/pte-infection-level-model-sample.csv"
}

if (!file.exists(input_file)) {
  stop("Could not find pte-infection-level-model-sample.csv.")
}

dat <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    avg_pte = sum_pte / n_edges,
    old_surprising = surprising %in% TRUE
  )

fit_rate_model <- function(dat, min_edges) {
  d <- dat %>%
    filter(n_edges >= min_edges) %>%
    mutate(agecat = factor(agecat))

  fit <- lm(
    avg_pte ~ ldens + asym + moi + agecat + lmosq,
    data = d
  )

  d <- d %>%
    mutate(
      pred_avg_pte = predict(fit),
      resid_avg_pte = avg_pte - pred_avg_pte,
      surprising_avg_pte = resid_avg_pte >= quantile(resid_avg_pte, 0.90)
    )

  coefs <- summary(fit)$coefficients

  summary_row <- tibble(
    min_edges = min_edges,
    infections = nrow(d),
    people = n_distinct(d$memid),
    households = n_distinct(d$hh),
    corr_avg_pte_edges = cor(d$avg_pte, d$n_edges),
    r_squared = summary(fit)$r.squared,
    surprising = sum(d$surprising_avg_pte),
    overlap_original_58 = sum(d$surprising_avg_pte & d$old_surprising),
    moi_estimate = coefs["moi", "Estimate"],
    moi_p = coefs["moi", "Pr(>|t|)"],
    mosquito_estimate = coefs["lmosq", "Estimate"],
    mosquito_p = coefs["lmosq", "Pr(>|t|)"]
  )

  list(data = d, fit = fit, summary = summary_row)
}

original_corr <- cor(dat$sum_pte, dat$n_edges)

fits <- lapply(c(1, 5, 10), function(x) fit_rate_model(dat, x))

results <- bind_rows(lapply(fits, `[[`, "summary"))

all_ids <- fits[[1]]$data %>%
  filter(surprising_avg_pte) %>%
  pull(sample_id_human)

five_ids <- fits[[2]]$data %>%
  filter(surprising_avg_pte) %>%
  pull(sample_id_human)

ten_ids <- fits[[3]]$data %>%
  filter(surprising_avg_pte) %>%
  pull(sample_id_human)

stability <- tibble(
  comparison = c(
    "All vs >=5 links",
    "All vs >=10 links",
    ">=5 vs >=10 links",
    "Surprising in all three"
  ),
  overlap = c(
    length(intersect(all_ids, five_ids)),
    length(intersect(all_ids, ten_ids)),
    length(intersect(five_ids, ten_ids)),
    length(Reduce(intersect, list(all_ids, five_ids, ten_ids)))
  )
)

cat("\nOriginal summed P(TE) correlation with link count:",
    round(original_corr, 3), "\n\n")

print(results, n = Inf, width = Inf)

cat("\nSurprise-set stability:\n")
print(stability, n = Inf)

dir.create(file.path("data", "processed"), recursive = TRUE, showWarnings = FALSE)

write_csv(results, file.path("data", "processed", "pte-per-link-sensitivity-summary.csv"))
write_csv(stability, file.path("data", "processed", "pte-per-link-surprise-stability.csv"))

for (i in seq_along(fits)) {
  cutoff <- c("all", "min5", "min10")[i]

  write_csv(
    fits[[i]]$data %>%
      filter(surprising_avg_pte) %>%
      arrange(desc(resid_avg_pte)),
    file.path("data", "processed", paste0("pte-surprising-per-link-", cutoff, ".csv"))
  )
}
