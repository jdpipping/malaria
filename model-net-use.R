library(tidyverse)
library(lubridate)

cohort_path <- "data/raw/malaria/OnceBittenMalariaCohort_jul2020-sep2021_withhaps.csv"
mosquito_path <- "data/raw/malaria/OnceBittenMosquito-Immediately processed summary.csv"
score_path <- "data/processed/drop-ssnet-net-use-surprise-scores.csv"
summary_path <- "data/processed/drop-ssnet-net-use-surprise-summary.csv"

cv_folds <- 5
cv_seed <- 20260831
surprise_cutoff <- 0.30
rare_level_min <- 20

to_num <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

to_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}

parse_date <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_

  suppressWarnings(
    as.Date(
      parse_date_time(
        x,
        orders = c("ymd", "mdy", "dmy", "Ymd", "mdY", "dmY"),
        quiet = TRUE
      )
    )
  )
}

first_num <- function(x) {
  x <- to_num(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else x[1]
}

first_chr <- function(x) {
  x <- to_chr(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_character_ else x[1]
}

mean_na <- function(x) {
  x <- to_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

auc <- function(y, p) {
  keep <- y %in% c(0, 1) & is.finite(p)
  y <- y[keep]
  p <- p[keep]

  n1 <- sum(y == 1)
  n0 <- sum(y == 0)

  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }

  r <- rank(p, ties.method = "average")

  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

is_checked <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1", "true", "checked", "yes")
}

choice_number <- function(x) {
  out <- sub("^.*___", "", x)

  if (!identical(out, x)) {
    return(suppressWarnings(as.integer(sub("[^0-9].*$", "", out))))
  }

  hit <- regmatches(x, regexpr("[0-9]+(?=[^0-9]*$)", x, perl = TRUE))

  if (length(hit) == 0 || hit == "") {
    return(NA_integer_)
  }

  suppressWarnings(as.integer(hit))
}

ensure_columns <- function(data, vars) {
  for (v in setdiff(vars, names(data))) {
    data[[v]] <- NA
  }

  data
}

make_code_factor <- function(x, prefix) {
  x <- to_chr(x)
  x[is.na(x)] <- "Missing"
  factor(paste0(prefix, "_", x))
}

collapse_rare <- function(x, min_n = rare_level_min) {
  x <- as.character(x)
  counts <- table(x)
  rare <- names(counts)[counts < min_n]
  x[x %in% rare] <- "OtherRare"
  factor(x)
}

make_group_folds <- function(person_id, k, seed) {
  people <- unique(as.character(person_id))

  set.seed(seed)
  assignment <- sample(rep_len(seq_len(k), length(people)))
  names(assignment) <- people

  unname(assignment[as.character(person_id)])
}

impute_train_test <- function(train, test, vars) {
  for (v in vars) {
    value <- median(train[[v]], na.rm = TRUE)

    if (!is.finite(value)) {
      value <- 0
    }

    train[[v]][is.na(train[[v]])] <- value
    test[[v]][is.na(test[[v]])] <- value
  }

  list(train = train, test = test)
}

varying_predictors <- function(data, vars) {
  vars[
    vapply(
      vars,
      function(v) {
        x <- data[[v]]
        x <- x[is.finite(x)]
        length(unique(x)) > 1
      },
      logical(1)
    )
  ]
}

estimable_predictors <- function(data, vars) {
  vars <- varying_predictors(data, vars)

  if (length(vars) == 0) {
    return(character(0))
  }

  x <- cbind(
    `(Intercept)` = 1,
    as.matrix(data[, vars, drop = FALSE])
  )

  q <- qr(x)
  keep <- colnames(x)[q$pivot[seq_len(q$rank)]]

  intersect(vars, keep)
}

d <- read_csv(
  cohort_path,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

mosquito <- read_csv(
  mosquito_path,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

d <- ensure_columns(
  d,
  c(
    "memid", "gender",
    "birth_year", "birth_month", "birth_day",
    "reported_age", "roster_first_form_date",
    "educ_level", "employment", "relationship",
    "village_name", "village_name_2", "village_name_3",
    "bldg_num", "net", "ss_slpspace_added", "ss_revision_date"
  )
)

monthly <- d %>%
  filter(
    !is.na(visit_date_monthly),
    trimws(as.character(visit_date_monthly)) != "",
    !is.na(member_monthly),
    trimws(as.character(member_monthly)) != ""
  ) %>%
  transmute(
    household_id = trimws(as.character(hh_id)),
    member_id = trimws(as.character(member_monthly)),
    person_id = paste(household_id, member_id, sep = "_"),
    visit_date = parse_date(visit_date_monthly),
    net_use = to_num(slept_net),
    pcr = to_num(pf_pcr_infection_status)
  ) %>%
  arrange(person_id, visit_date) %>%
  group_by(person_id) %>%
  mutate(
    next_visit_date = lead(visit_date),
    next_pcr = lead(pcr),
    days_to_next = as.numeric(next_visit_date - visit_date),
    prior_pcr_positive_count = lag(
      cumsum(coalesce(pcr == 1, FALSE)),
      default = 0
    )
  ) %>%
  ungroup() %>%
  mutate(
    year_month = format(visit_date, "%Y-%m"),
    month_start = floor_date(visit_date, "month"),
    month_of_year = month(visit_date)
  )

eligible <- monthly %>%
  filter(
    pcr == 0,
    net_use %in% c(0, 1),
    between(days_to_next, 20, 45),
    next_pcr %in% c(0, 1)
  ) %>%
  mutate(row_id = row_number())

stopifnot(
  nrow(monthly) == 5825,
  nrow(eligible) == 3459,
  sum(eligible$net_use == 1) == 2458,
  sum(eligible$net_use == 0) == 1001
)

roster <- d %>%
  filter(
    !is.na(memid),
    trimws(as.character(memid)) != "",
    !is.na(hh_id),
    trimws(as.character(hh_id)) != ""
  ) %>%
  mutate(
    household_id = trimws(as.character(hh_id)),
    member_id = trimws(as.character(memid))
  ) %>%
  group_by(household_id, member_id) %>%
  summarise(
    gender_raw = first_chr(gender),
    birth_year = first_num(birth_year),
    birth_month = first_num(birth_month),
    birth_day = first_num(birth_day),
    reported_age = first_num(reported_age),
    roster_date_raw = first_chr(roster_first_form_date),
    education_raw = first_chr(educ_level),
    employment_raw = first_chr(employment),
    relationship_raw = first_chr(relationship),
    .groups = "drop"
  ) %>%
  mutate(
    sex = case_when(
      gender_raw == "1" ~ "Male",
      gender_raw == "2" ~ "Female",
      TRUE ~ NA_character_
    ),
    roster_date = parse_date(roster_date_raw),
    dob_month = if_else(birth_month %in% 1:12, birth_month, 7),
    dob_day = if_else(birth_day %in% 1:31, birth_day, 15),
    dob = make_date(
      as.integer(birth_year),
      as.integer(dob_month),
      as.integer(dob_day)
    ),

    education = make_code_factor(education_raw, "education"),
    employment = make_code_factor(employment_raw, "employment"),
    household_role = make_code_factor(relationship_raw, "role")
  )

household_village <- d %>%
  mutate(
    household_id = trimws(as.character(hh_id)),
    village_raw = coalesce(
      to_chr(village_name),
      to_chr(village_name_2),
      to_chr(village_name_3)
    )
  ) %>%
  filter(
    !is.na(household_id),
    household_id != ""
  ) %>%
  group_by(household_id) %>%
  summarise(
    village_raw = first_chr(village_raw),
    .groups = "drop"
  ) %>%
  mutate(
    village = case_when(
      village_raw == "1" ~ "Sitabicha",
      village_raw == "2" ~ "Lurare",
      village_raw == "3" ~ "Nangili",
      village_raw == "4" ~ "Maruti",
      village_raw == "5" ~ "Kinesamo",
      TRUE ~ village_raw
    ),
    village = factor(village)
  ) %>%
  select(
    household_id,
    village
  )

household_size <- roster %>%
  count(household_id, name = "household_size")

eligible <- eligible %>%
  left_join(
    roster %>%
      select(
        household_id, member_id, sex, dob, reported_age, roster_date,
        education, employment, household_role
      ),
    by = c("household_id", "member_id")
  ) %>%
  left_join(
    household_size,
    by = "household_id"
  ) %>%
  left_join(
    household_village,
    by = "household_id"
  ) %>%
  mutate(
    age_from_dob = as.numeric(visit_date - dob) / 365.2425,
    age_from_report = if_else(
      !is.na(reported_age) & !is.na(roster_date),
      reported_age + as.numeric(visit_date - roster_date) / 365.2425,
      reported_age
    ),
    age = coalesce(age_from_dob, age_from_report)
  )

stopifnot(
  sum(eligible$sex == "Male") == 1492,
  sum(eligible$sex == "Female") == 1967
)

asset_candidates <- grep("^household_items", names(d), value = TRUE)
asset_choice <- vapply(asset_candidates, choice_number, integer(1))
asset_cols <- asset_candidates[
  !is.na(asset_choice) & asset_choice %in% 1:8
]

stopifnot(length(asset_cols) == 8)

asset_work <- d %>%
  mutate(
    source_row = row_number(),
    household_id = trimws(as.character(hh_id))
  )

asset_matrix <- sapply(
  asset_work[asset_cols],
  function(x) as.numeric(is_checked(x))
)

asset_observed <- sapply(
  asset_work[asset_cols],
  function(x) as.numeric(!is.na(to_chr(x)))
)

asset_work$asset_count <- rowSums(asset_matrix, na.rm = TRUE)
asset_work$asset_observed <- rowSums(asset_observed, na.rm = TRUE)

household_assets <- asset_work %>%
  filter(
    !is.na(household_id),
    household_id != "",
    asset_observed > 0
  ) %>%
  arrange(household_id, desc(asset_observed), source_row) %>%
  group_by(household_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    household_id,
    household_assets = asset_count
  )

eligible <- eligible %>%
  left_join(
    household_assets,
    by = "household_id"
  )

building_parts <- vector("list", 10)

for (b in 1:10) {
  gap_col <- paste0("gap", b)
  netting_col <- paste0("netting", b)

  d <- ensure_columns(d, c(gap_col, netting_col))

  building_parts[[b]] <- d %>%
    transmute(
      household_id = trimws(as.character(hh_id)),
      building_num = b,
      gap = to_num(.data[[gap_col]]),
      netting = to_num(.data[[netting_col]])
    ) %>%
    filter(
      !is.na(household_id),
      household_id != ""
    ) %>%
    group_by(household_id, building_num) %>%
    summarise(
      gap = first_num(gap),
      netting = first_num(netting),
      .groups = "drop"
    ) %>%
    mutate(
      closed_eaves = case_when(
        gap == 0 ~ 1,
        gap == 1 ~ 0,
        TRUE ~ NA_real_
      ),
      screened_windows = case_when(
        netting == 0 ~ 0,
        netting == 1 ~ 0.5,
        netting == 2 ~ 1,
        TRUE ~ NA_real_
      ),

      housing_protection = case_when(
        is.na(closed_eaves) & is.na(screened_windows) ~ NA_real_,
        TRUE ~ rowMeans(
          cbind(closed_eaves, screened_windows),
          na.rm = TRUE
        ) * 2
      )
    )
}

building_lookup <- bind_rows(building_parts)

household_protection <- building_lookup %>%
  group_by(household_id) %>%
  summarise(
    household_protection = mean_na(housing_protection),
    .groups = "drop"
  )

member_candidates <- grep("^mem_inspace", names(d), value = TRUE)
member_choice <- vapply(member_candidates, choice_number, integer(1))
keep_member <- !is.na(member_choice) & member_choice %in% 1:30
member_cols <- member_candidates[keep_member]
member_ids <- member_choice[keep_member]

stopifnot(length(member_cols) == 30)

space_base <- d %>%
  transmute(
    source_row = row_number(),
    household_id = trimws(as.character(hh_id)),
    building_num = to_num(bldg_num),
    net_available = case_when(
      to_num(net) == 1 ~ 1,
      to_num(net) == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    date_added = parse_date(ss_slpspace_added),
    date_revised = parse_date(ss_revision_date)
  )

membership_matrix <- sapply(
  d[member_cols],
  function(x) as.numeric(is_checked(x))
)

space_base$sleeping_space_crowding <- if_else(
  rowSums(membership_matrix, na.rm = TRUE) > 0,
  rowSums(membership_matrix, na.rm = TRUE),
  NA_real_
)

space_members <- map2(
  member_cols,
  member_ids,
  function(member_col, member_id) {
    keep <- is_checked(d[[member_col]])

    space_base[keep, , drop = FALSE] %>%
      mutate(
        member_id = as.character(member_id),
        effective_date = coalesce(date_revised, date_added)
      )
  }
) %>%
  bind_rows() %>%
  filter(
    !is.na(household_id),
    household_id != "",
    !is.na(member_id),
    member_id != ""
  ) %>%
  left_join(
    building_lookup %>%
      select(
        household_id,
        building_num,
        housing_protection
      ),
    by = c("household_id", "building_num")
  )

space_at_t <- eligible %>%
  select(
    row_id,
    household_id,
    member_id,
    visit_date
  ) %>%
  inner_join(
    space_members,
    by = c("household_id", "member_id"),
    relationship = "many-to-many"
  ) %>%
  filter(
    is.na(effective_date) |
      effective_date <= visit_date
  ) %>%
  mutate(
    rank_date = coalesce(
      effective_date,
      as.Date("1900-01-01")
    )
  ) %>%
  arrange(
    row_id,
    desc(rank_date),
    desc(source_row)
  ) %>%
  group_by(row_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    row_id,
    net_available,
    sleeping_space_crowding,
    housing_protection
  )

eligible <- eligible %>%
  left_join(
    space_at_t,
    by = "row_id"
  ) %>%
  left_join(
    household_protection,
    by = "household_id"
  ) %>%
  mutate(
    housing_protection = coalesce(
      housing_protection,
      household_protection
    )
  ) %>%
  select(-household_protection)

mosquito_month <- mosquito %>%
  transmute(
    household_id = trimws(as.character(household_id)),
    collection_date = parse_date(collection_date_1a),
    female_anopheles = to_num(femalea_total)
  ) %>%
  filter(
    !is.na(household_id),
    household_id != "",
    !is.na(collection_date),
    !is.na(female_anopheles)
  ) %>%
  mutate(
    month_start = floor_date(collection_date, "month")
  ) %>%
  group_by(household_id, month_start) %>%
  summarise(
    mosquito_abundance = mean(female_anopheles),
    .groups = "drop"
  )

eligible <- eligible %>%
  mutate(
    previous_month_start = month_start %m-% months(1)
  ) %>%
  left_join(
    mosquito_month %>%
      rename(
        previous_month_start = month_start
      ),
    by = c(
      "household_id",
      "previous_month_start"
    )
  ) %>%
  mutate(
    prior_mosquito_abundance = log1p(mosquito_abundance),
    education = collapse_rare(education),
    employment = collapse_rare(employment),
    household_role = collapse_rare(household_role),
    sleeping_space_net = factor(
      case_when(
        net_available == 0 ~ "No",
        net_available == 1 ~ "Yes",
        TRUE ~ "Missing"
      )
    ),
    month_season = factor(
      month_of_year,
      levels = 1:12
    )
  )

numeric_covariates <- c(
  "age",
  "household_size",
  "household_assets",
  "housing_protection",
  "sleeping_space_crowding",
  "prior_pcr_positive_count",
  "prior_mosquito_abundance"
)

factor_covariates <- c(
  "education",
  "employment",
  "household_role",
  "month_season",
  "village"
)

stopifnot(
  length(numeric_covariates) +
    length(factor_covariates) == 12
)

factor_matrix <- model.matrix(
  reformulate(factor_covariates),
  data = eligible
)

factor_matrix <- factor_matrix[
  ,
  colnames(factor_matrix) != "(Intercept)",
  drop = FALSE
]

colnames(factor_matrix) <- make.names(
  colnames(factor_matrix),
  unique = TRUE
)

dummy_covariates <- colnames(factor_matrix)

model_data <- bind_cols(
  eligible,
  as.data.frame(
    factor_matrix,
    check.names = FALSE
  )
)

model_covariates <- c(
  numeric_covariates,
  dummy_covariates
)

fit_sex <- function(data, sex_name, seed) {
  s <- data %>%
    filter(sex == sex_name)

  full_data <- s

  for (v in numeric_covariates) {
    value <- median(
      full_data[[v]],
      na.rm = TRUE
    )

    if (!is.finite(value)) {
      value <- 0
    }

    full_data[[v]][is.na(full_data[[v]])] <- value
  }

  full_covariates <- estimable_predictors(
    full_data,
    model_covariates
  )

  full_fit <- glm(
    reformulate(
      full_covariates,
      response = "net_use"
    ),
    family = binomial(),
    data = full_data
  )

  pscore <- as.numeric(
    predict(
      full_fit,
      newdata = full_data,
      type = "response"
    )
  )

  fold <- make_group_folds(
    s$person_id,
    cv_folds,
    seed
  )

  pscore_oof <- rep(
    NA_real_,
    nrow(s)
  )

  for (k in seq_len(cv_folds)) {
    train <- s[
      fold != k,
      ,
      drop = FALSE
    ]

    test <- s[
      fold == k,
      ,
      drop = FALSE
    ]

    imputed <- impute_train_test(
      train,
      test,
      numeric_covariates
    )

    train <- imputed$train
    test <- imputed$test

    fold_covariates <- estimable_predictors(
      train,
      model_covariates
    )

    fit <- glm(
      reformulate(
        fold_covariates,
        response = "net_use"
      ),
      family = binomial(),
      data = train
    )

    pscore_oof[fold == k] <- as.numeric(
      predict(
        fit,
        newdata = test,
        type = "response"
      )
    )
  }

  if (any(!is.finite(pscore_oof))) {
    stop(
      sex_name,
      ": nonfinite out-of-fold propensity scores"
    )
  }

  s %>%
    mutate(
      pscore = pscore,
      pscore_oof = pscore_oof,
      surprise = abs(net_use - pscore_oof),
      surprise_full = abs(net_use - pscore),
      surprising_net_user = net_use == 1 &
        pscore_oof <= surprise_cutoff,
      surprising_control = net_use == 0 &
        pscore_oof >= 1 - surprise_cutoff,
      surprising_net_user_full = net_use == 1 &
        pscore <= surprise_cutoff,
      surprising_control_full = net_use == 0 &
        pscore >= 1 - surprise_cutoff
    )
}

scored <- bind_rows(
  fit_sex(
    model_data,
    "Male",
    cv_seed + 1
  ),
  fit_sex(
    model_data,
    "Female",
    cv_seed + 2
  )
)

auc_summary <- scored %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    net_users = sum(net_use == 1),
    auc_full = auc(net_use, pscore),
    auc_oof = auc(net_use, pscore_oof),
    .groups = "drop"
  )

summarise_surprise <- function(data, flag, label, score_type) {
  x <- data %>%
    filter(.data[[flag]])

  tibble(
    score = score_type,
    group = label,
    total = nrow(x),
    male = sum(x$sex == "Male"),
    female = sum(x$sex == "Female"),
    people = n_distinct(x$person_id),
    households = n_distinct(x$household_id),
    calendar_months = n_distinct(x$year_month)
  )
}

surprise_summary <- bind_rows(
  summarise_surprise(
    scored,
    "surprising_net_user",
    "net user <= 0.30",
    "person-grouped OOF"
  ),
  summarise_surprise(
    scored,
    "surprising_control",
    "control >= 0.70",
    "person-grouped OOF"
  ),
  summarise_surprise(
    scored,
    "surprising_net_user_full",
    "net user <= 0.30",
    "full fitted"
  ),
  summarise_surprise(
    scored,
    "surprising_control_full",
    "control >= 0.70",
    "full fitted"
  )
)

primary_surprise <- surprise_summary %>%
  filter(score == "person-grouped OOF")

model_summary <- bind_rows(
  auc_summary %>%
    transmute(
      result = paste0("AUC ", sex),
      value = auc_oof
    ),
  primary_surprise %>%
    select(
      group,
      total,
      male,
      female,
      people,
      households,
      calendar_months
    ) %>%
    pivot_longer(
      cols = -group,
      names_to = "measure",
      values_to = "value"
    ) %>%
    transmute(
      result = paste(group, measure),
      value
    )
)

dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  scored %>%
    select(
      row_id,
      household_id,
      member_id,
      person_id,
      visit_date,
      year_month,
      sex,
      net_use,
      next_pcr,
      age,
      education,
      employment,
      household_role,
      household_size,
      household_assets,
      housing_protection,
      sleeping_space_net,
      sleeping_space_crowding,
      prior_pcr_positive_count,
      prior_mosquito_abundance,
      month_season,
      village,
      pscore,
      pscore_oof,
      surprise,
      surprise_full,
      surprising_net_user,
      surprising_control,
      surprising_net_user_full,
      surprising_control_full
    ),
  score_path
)

write_csv(
  model_summary,
  summary_path
)

print(
  auc_summary %>%
    mutate(
      auc_full = round(auc_full, 3),
      auc_oof = round(auc_oof, 3)
    )
)

print(
  surprise_summary,
  n = Inf,
  width = Inf
)

proxy_monthly <- d %>%
  filter(
    !is.na(visit_date_monthly),
    trimws(as.character(visit_date_monthly)) != "",
    !is.na(member_monthly),
    trimws(as.character(member_monthly)) != ""
  ) %>%
  transmute(
    household_id = trimws(as.character(hh_id)),
    member_id    = trimws(as.character(member_monthly)),
    person_id    = paste(household_id, member_id, sep = "_"),
    visit_date   = parse_date(visit_date_monthly),
    month_start  = floor_date(visit_date, "month"),
    net_use      = to_num(slept_net),
    slept_times  = ifelse(to_num(slept_times) %in% 0:7, to_num(slept_times), NA_real_)
  ) %>%
  arrange(person_id, visit_date) %>%
  group_by(person_id) %>%
  mutate(
    prior_personal_net_use = lag(net_use),
    prior_net_frequency    = lag(slept_times)
  ) %>%
  ungroup()

hh_prior_rate <- proxy_monthly %>%
  group_by(household_id, month_start) %>%
  summarise(rate = mean(net_use, na.rm = TRUE), .groups = "drop") %>%
  mutate(month_start = month_start %m+% months(1)) %>%
  rename(prior_household_net_rate = rate)

proxy_monthly <- proxy_monthly %>%
  left_join(hh_prior_rate, by = c("household_id", "month_start"))

screen_data <- eligible %>%
  left_join(
    proxy_monthly %>%
      select(person_id, visit_date,
             prior_personal_net_use, prior_net_frequency, prior_household_net_rate),
    by = c("person_id", "visit_date")
  )

screen_covariates <- c(
  "age", "education", "employment", "household_role", "household_size",
  "household_assets", "housing_protection", "sleeping_space_net",
  "sleeping_space_crowding", "prior_pcr_positive_count", "prior_mosquito_abundance",
  "month_season", "village",
  "prior_personal_net_use", "prior_net_frequency", "prior_household_net_rate"
)
retained_covariates <- setdiff(
  screen_covariates,
  c("sleeping_space_net", "prior_personal_net_use",
    "prior_net_frequency", "prior_household_net_rate")
)

auc_from_rank <- function(p, yy) {
  ok <- is.finite(p); p <- p[ok]; yy <- yy[ok]
  n1 <- sum(yy == 1); n0 <- sum(yy == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[yy == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

univ_oof_auc <- function(nm) {
  x <- screen_data[[nm]]; yy <- screen_data$net_use; person <- screen_data$person_id
  ok <- !is.na(x) & !is.na(yy)
  dd <- data.frame(y = yy[ok], person = as.character(person[ok]))
  dd$x <- if (is.numeric(x[ok])) x[ok] else factor(x[ok])
  set.seed(cv_seed)
  ppl  <- unique(dd$person)
  fold <- setNames(sample(rep(seq_len(cv_folds), length.out = length(ppl))), ppl)
  dd$fold <- fold[dd$person]; dd$p <- NA_real_
  for (k in seq_len(cv_folds)) {
    te  <- dd$fold == k
    fit <- try(glm(y ~ x, family = binomial(), data = dd[!te, , drop = FALSE]), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      pr <- try(predict(fit, dd[te, , drop = FALSE], type = "response"), silent = TRUE)
      if (!inherits(pr, "try-error")) dd$p[te] <- as.numeric(pr)
    }
  }
  auc_from_rank(dd$p, dd$y)
}

cramers_v <- function(x, yy) {
  t  <- table(x, yy); c2 <- suppressWarnings(chisq.test(t)$statistic)
  n  <- sum(t); k <- min(nrow(t), ncol(t))
  as.numeric(sqrt(c2 / (n * max(1, k - 1))))
}
treat_assoc <- function(nm) {
  x <- screen_data[[nm]]; yy <- screen_data$net_use; ok <- !is.na(x)
  if (is.numeric(x)) abs(cor(x[ok], yy[ok])) else cramers_v(as.character(x[ok]), yy[ok])
}

covariate_screen <- tibble(
  covariate     = screen_covariates,
  oof_auc       = round(vapply(screen_covariates, univ_oof_auc, numeric(1)), 3),
  cor_treatment = round(vapply(screen_covariates, treat_assoc,  numeric(1)), 3),
  included      = ifelse(screen_covariates %in% retained_covariates, "Y", "N")
) %>%
  arrange(desc(oof_auc))

print(covariate_screen, n = Inf)
