library(tidyverse)
library(ranger)

CSV        <- "data/processed/drop-ssnet-net-use-surprise-scores.csv"
OUT_SCORES <- "data/processed/rf-net-use-surprise-scores.csv"
OUT_PNG    <- "data/processed/rf_net_use_two_zone.png"

d <- read_csv(CSV, show_col_types = FALSE)

covars <- c("age","education","employment","household_role","household_size",
            "household_assets","housing_protection","sleeping_space_crowding",
            "prior_pcr_positive_count","prior_mosquito_abundance","month_season","village")
num    <- c("age","household_size","household_assets","housing_protection",
            "sleeping_space_crowding","prior_pcr_positive_count","prior_mosquito_abundance")
catv   <- setdiff(covars, num)

d <- d %>%
  mutate(across(all_of(num),  ~ { x <- suppressWarnings(as.numeric(.)); ifelse(is.na(x), median(x, na.rm = TRUE), x) }),
         across(all_of(catv), as.factor),
         net_use_f = factor(net_use, levels = c(0, 1)))

auc_rank <- function(p, y) {                     # AUC = Mann-Whitney, no packages
  ok <- is.finite(p); p <- p[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

set.seed(1)
d$pscore_rf <- NA_real_
for (sx in c("Male", "Female")) {
  idx  <- which(d$sex == sx)
  sub  <- d[idx, ]
  ppl  <- unique(sub$person_id)                  # fold at the PERSON level (no leakage)
  fold <- setNames(sample(rep(1:5, length.out = length(ppl))), ppl)
  sub$fold <- fold[sub$person_id]
  pr <- rep(NA_real_, nrow(sub))
  for (k in 1:5) {
    tr <- sub[sub$fold != k, ]
    te <- which(sub$fold == k)
    rf <- ranger(x = tr[, covars], y = tr$net_use_f,
                 probability = TRUE, num.trees = 500, min.node.size = 10, seed = 1)
    pr[te] <- predict(rf, data = sub[te, covars])$predictions[, "1"]
  }
  d$pscore_rf[idx] <- pr
  cat(sprintf("%s: random-forest OOF AUC = %.3f  (n = %d)\n", sx, auc_rank(pr, sub$net_use), length(idx)))
}

write_csv(d, OUT_SCORES)

# ---- same two-zone plot, x = random-forest OOF probability ----
COL_NET <- "#9D0C0C"; COL_NON <- "#041851"; A <- 0.65
LAB_NET <- "Net Users (n=2,458)"; LAB_NON <- "Non-Users (n=1,001)"
brks <- seq(0, 1, by = 0.05)
bins_tbl <- d %>%
  mutate(grp = ifelse(net_use == 1, LAB_NET, LAB_NON),
         bin = cut(pscore_rf, breaks = brks, include.lowest = TRUE, right = FALSE)) %>%
  group_by(grp, bin, .drop = FALSE) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(grp) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(xmid = brks[as.integer(bin)] + 0.025)

p <- ggplot() +
  geom_col(data = filter(bins_tbl, grp == LAB_NON), aes(xmid, pct, fill = LAB_NON), width = 0.05, alpha = A) +
  geom_col(data = filter(bins_tbl, grp == LAB_NET), aes(xmid, pct, fill = LAB_NET), width = 0.05, alpha = A) +
  scale_fill_manual(values = setNames(c(COL_NET, COL_NON), c(LAB_NET, LAB_NON)),
                    breaks = c(LAB_NET, LAB_NON), name = NULL,
                    guide = guide_legend(override.aes = list(alpha = c(A, A)))) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = seq(0, 40, 10), labels = function(v) paste0(v, "%"),
                     limits = c(0, 40), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Predicted Net-Use Probability (Random Forest), Colored by Actual Use",
       x = "Predicted Probability of Net Use (random forest, out-of-fold)", y = "Percent of Group") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 8)),
        legend.position = "top", legend.key = element_rect(fill = "white", color = NA),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "#D9D9D9", linewidth = 0.4),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
        panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white", color = NA))

ggsave(OUT_PNG, p, width = 10.85, height = 6.43, dpi = 200)
cat(sprintf("\nsaved -> %s\n", OUT_PNG))
