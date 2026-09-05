suppressMessages({library(dplyr); library(ggplot2)})

CSV <- "data/processed/drop-ssnet-net-use-surprise-scores.csv"
OUT <- "data/processed/drop-ssnet_net_use_two_zone.png"

df <- read.csv(CSV, stringsAsFactors = FALSE)

n_net <- sum(df$net_use == 1); n_non <- sum(df$net_use == 0)
n_surp_t <- sum(df$net_use == 1 & df$pscore_oof <= 0.30)
n_surp_c <- sum(df$net_use == 0 & df$pscore_oof >= 0.70)
rng_ok <- min(df$pscore_oof) >= 0 && max(df$pscore_oof) <= 1
cat(sprintf("net-users: %d (expect 2458) | non-users: %d (expect 1001)\n", n_net, n_non))
cat(sprintf("surprising treated (<=0.30): %d (expect 122) | surprising controls (>=0.70): %d (expect 366)\n",
            n_surp_t, n_surp_c))
if (!(n_net == 2458 && n_non == 1001 && n_surp_t == 122 && n_surp_c == 366 && rng_ok))
  stop("CSV checks failed — inspect the file before plotting.")

COL_NET  <- "#9D0C0C"; COL_NON <- "#041851"
COL_ZONE <- "#DADADA"; ZA <- 0.55
COL_LINE <- "#8A8A8A"
A <- 0.65
LAB_NET <- "Net Users (n=2,458)"; LAB_NON <- "Non-Users (n=1,001)"; LAB_ZONE <- "Surprise zones"

zones <- data.frame(xmin = c(0, 0.70), xmax = c(0.30, 1.00),
                    ymin = 0, ymax = Inf, zone = LAB_ZONE)

brks <- seq(0, 1, by = 0.05)
bins_tbl <- df %>%
  mutate(grp = ifelse(net_use == 1, LAB_NET, LAB_NON),
         bin = cut(pscore_oof, breaks = brks, include.lowest = TRUE, right = FALSE)) %>%
  group_by(grp, bin, .drop = FALSE) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(grp) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(xmid = brks[as.integer(bin)] + 0.025)

p <- ggplot() +
  geom_rect(data = zones, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = zone),
            alpha = ZA) +
  geom_col(data = filter(bins_tbl, grp == LAB_NON),
           aes(x = xmid, y = pct, fill = LAB_NON), width = 0.05, alpha = A) +
  geom_col(data = filter(bins_tbl, grp == LAB_NET),
           aes(x = xmid, y = pct, fill = LAB_NET), width = 0.05, alpha = A) +
  geom_vline(xintercept = c(0.30, 0.70), linetype = "dashed", color = COL_LINE, linewidth = 0.6) +
  scale_fill_manual(
    values = setNames(c(COL_NET, COL_NON, COL_ZONE), c(LAB_NET, LAB_NON, LAB_ZONE)),
    breaks = c(LAB_NET, LAB_NON, LAB_ZONE), name = NULL,
    guide  = guide_legend(override.aes = list(alpha = c(A, A, ZA)))) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = seq(0, 40, 10), labels = function(v) paste0(v, "%"),
                     limits = c(0, 40), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Predicted Net-Use Probability, Colored by Actual Use",
       x = "Predicted Probability of Net Use", y = "Percent of Group") +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 8)),
    legend.position = "top",
    legend.key = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#D9D9D9", linewidth = 0.4),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA))

ggsave(OUT, p, width = 10.85, height = 6.43, dpi = 200)
cat(sprintf("\nsaved -> %s\n", OUT))
