# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(tidyverse)
library(scales)
library(viridis)
library(scico)
library(cowplot)
library(ggpubr)
library(rstatix)

# read in raw data
df <- read_tsv("./raw_quantifications_Fig2to5_panelBtoE.tsv")

# tidy the table
df_mod <- df %>%
  pivot_longer(
    cols = starts_with("DIV"),
    names_to = "timepoint_gene",
    values_to = "proportion"
  ) %>%
  separate(timepoint_gene, into = c("timepoint", "gene"), sep = "_") %>%
  mutate(diagnosis = str_replace(sample, "_.*", "")) %>%
  select(diagnosis, sample, replicate, timepoint, gene, proportion)
write_tsv(df_mod, "./raw_quantifications_Fig2to5_panelBtoE_tidy.tsv")

# get mean, median, SD and SEM across replicates by sample, timepoint, and gene
sem <- function(x) {
  N <- sum(!is.na(x))
  sd(x, na.rm = T) / sqrt(N)
}

ci <- function(x, conf.interval = 0.95) {
  N <- sum(!is.na(x))
  sem <- sd(x, na.rm = T) / sqrt(N)
  qt(1 - (1 - conf.interval) / 2, N - 1) * sem
}

df_summarized <- df_mod %>%
  group_by(diagnosis, sample, timepoint, gene) %>%
  summarise(
    proportion_mean = mean(proportion, na.rm = TRUE),
    proportion_median = median(proportion, na.rm = TRUE),
    proportion_sd = sd(proportion, na.rm = TRUE),
    proportion_sem = sem(proportion),
    proportion_ci = ci(proportion)
  ) %>%
  ungroup()
write_tsv(df_summarized, "./quantifications_Fig2to5_panelBtoE_summarized_by_sample_timepoint_gene.tsv")

# set up colors
diagnosis_colors <- c(
  "Ctrl" = "#9DC08B",
  "EPI" = "#00A0B0",
  "MIC" = "#655DBB",
  "ID" = "#CC333F",
  "PMG" = "#EDC951"
)

# define a function to make plots and save statistics test results
run_plot_and_test <- function(x_gene, x_timepoint, y_gene, y_timepoint, diagnosis_1, diagnosis_2) {
  message(str_glue("Running plot for {x_gene} at {x_timepoint} vs {y_gene} at {y_timepoint} in {diagnosis_1} vs {diagnosis_2}"))

  x_to_plot <- paste(x_timepoint, x_gene, sep = "_")
  y_to_plot <- paste(y_timepoint, y_gene, sep = "_")

  df_mod_to_plot <- df_mod %>%
    filter(
      diagnosis %in% c(diagnosis_1, diagnosis_2),
      (gene == x_gene & timepoint == x_timepoint) | (gene == y_gene & timepoint == y_timepoint)
    ) %>%
    mutate(timepoint_gene = paste(timepoint, gene, sep = "_")) %>%
    select(-timepoint, -gene) %>%
    pivot_wider(names_from = timepoint_gene, values_from = proportion)

  df_summarized_to_plot <- df_summarized %>%
    filter(
      diagnosis %in% c(diagnosis_1, diagnosis_2),
      (gene == x_gene & timepoint == x_timepoint) | (gene == y_gene & timepoint == y_timepoint)
    ) %>%
    mutate(timepoint_gene = paste(timepoint, gene, sep = "_")) %>%
    select(-timepoint, -gene) %>%
    pivot_wider(
      names_from = timepoint_gene,
      values_from = c(proportion_mean, proportion_median, proportion_sd, proportion_sem, proportion_ci)
    )

  # main plot, with each datapoint from each replicate underneath
  # and the average for each individual sample on top
  p <- df_mod_to_plot %>%
    ggplot(aes(x = .data[[x_to_plot]], y = .data[[y_to_plot]]))
  p_main <- p +
    geom_point(aes(color = diagnosis), alpha = 0.2, size = 1) +
    geom_point(
      data = df_summarized_to_plot,
      aes(
        x = .data[[paste("proportion_mean", x_to_plot, sep = "_")]],
        y = .data[[paste("proportion_mean", y_to_plot, sep = "_")]],
        fill = diagnosis
      ),
      color = "black",
      shape = 21,
      size = 2
    ) +
    xlab(str_glue("Proportion of cells expressing {x_gene} at {x_timepoint}")) +
    ylab(str_glue("Proportion of cells expressing {y_gene} at {y_timepoint}")) +
    scale_fill_manual(values = diagnosis_colors) +
    scale_color_manual(values = diagnosis_colors) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )

  ggsave(
    str_glue(
      "./plots/plot_Fig2to5_panelBtoE/cell_count_proportion_",
      "{x_to_plot}_vs_{y_to_plot}_in_{diagnosis_1}_vs_{diagnosis_2}_",
      "scatter_plus_average_per_sample.pdf"
    ),
    plot = p_main,
    device = cairo_pdf(),
    width = 4,
    height = 3,
    useDingbats = FALSE
  )
  dev.off()


  # add marginal density plot
  x_marginal_density <- axis_canvas(p_main, axis = "x") +
    geom_density(
      data = df_summarized_to_plot,
      aes(
        x = .data[[paste("proportion_mean", x_to_plot, sep = "_")]],
        y = after_stat(scaled),
        fill = diagnosis
      ),
      size = 0.5,
      alpha = 0.5
    ) +
    scale_fill_manual(values = diagnosis_colors)

  y_marginal_density <- axis_canvas(p_main, axis = "y", coord_flip = TRUE) +
    geom_density(
      data = df_summarized_to_plot,
      aes(
        x = .data[[paste("proportion_mean", y_to_plot, sep = "_")]],
        y = after_stat(scaled),
        fill = diagnosis
      ),
      size = 0.5,
      alpha = 0.5
    ) +
    coord_flip() +
    scale_fill_manual(values = diagnosis_colors)

  p1 <- insert_xaxis_grob(p_main, x_marginal_density, grid::unit(.1, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, y_marginal_density, grid::unit(.1, "null"), position = "right")
  p_density <- ggdraw(p2)

  ggsave(
    str_glue(
      "./plots/plot_Fig2to5_panelBtoE/cell_count_proportion_",
      "{x_to_plot}_vs_{y_to_plot}_in_{diagnosis_1}_vs_{diagnosis_2}_",
      "scatter_plus_average_per_sample_add_marginal_density.pdf"
    ),
    plot = p_density,
    device = cairo_pdf(),
    width = 4,
    height = 3,
    useDingbats = FALSE
  )
  dev.off()

  # add marginal box plot
  stat_test_x <- df_summarized_to_plot %>%
    t_test(reformulate("diagnosis", str_glue("proportion_mean_{x_to_plot}")), detailed = TRUE) %>%
    add_xy_position(x = str_glue("proportion_mean_{x_to_plot}"), dodge = 0.5) %>%
    mutate(gene = x_gene, timepoint = x_timepoint)
  stat_test_y <- df_summarized_to_plot %>%
    t_test(reformulate("diagnosis", str_glue("proportion_mean_{y_to_plot}")), detailed = TRUE) %>%
    add_xy_position(x = str_glue("proportion_mean_{y_to_plot}"), dodge = 0.5) %>%
    mutate(gene = y_gene, timepoint = y_timepoint)

  x_marginal_boxplot <- axis_canvas(p_main, axis = "x") +
    geom_boxplot(
      data = df_summarized_to_plot,
      aes(
        x = .data[[paste("proportion_mean", x_to_plot, sep = "_")]],
        color = diagnosis
      ),
      outlier.shape = NA,
      position = position_dodge(width = 0.5),
      width = 0.4,
      size = 0.5
    ) +
    annotate(
      "text",
      x = stat_test_x$y.position,
      y = 0,
      hjust = -0.1,
      label = str_glue("p = {stat_test_x$p}"),
      size = 1.5
    ) +
    annotate(
      "segment",
      x = stat_test_x$y.position,
      xend = stat_test_x$y.position,
      y = 0.125,
      yend = -0.125,
      size = 0.25
    ) +
    annotate(
      "segment",
      x = stat_test_x$y.position,
      xend = stat_test_x$y.position - 0.01,
      y = 0.125,
      yend = 0.125,
      size = 0.25
    ) +
    annotate(
      "segment",
      x = stat_test_x$y.position,
      xend = stat_test_x$y.position - 0.01,
      y = -0.125,
      yend = -0.125,
      size = 0.25
    ) +
    scale_color_manual(values = diagnosis_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

  y_marginal_boxplot <- axis_canvas(p_main, axis = "y", coord_flip = TRUE) +
    geom_boxplot(
      data = df_summarized_to_plot,
      aes(
        x = .data[[paste("proportion_mean", y_to_plot, sep = "_")]],
        color = diagnosis
      ),
      outlier.shape = NA,
      position = position_dodge(width = 0.5),
      width = 0.4,
      size = 0.5
    ) +
    annotate(
      "text",
      x = stat_test_y$y.position,
      y = 0,
      hjust = -0.1,
      label = str_glue("p = {stat_test_y$p}"),
      size = 1.5,
      angle = 90
    ) +
    annotate(
      "segment",
      x = stat_test_y$y.position,
      xend = stat_test_y$y.position,
      y = 0.125,
      yend = -0.125,
      size = 0.25
    ) +
    annotate(
      "segment",
      x = stat_test_y$y.position,
      xend = stat_test_y$y.position - 0.01,
      y = 0.125,
      yend = 0.125,
      size = 0.25
    ) +
    annotate(
      "segment",
      x = stat_test_y$y.position,
      xend = stat_test_y$y.position - 0.01,
      y = -0.125,
      yend = -0.125,
      size = 0.25
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    coord_flip() +
    scale_color_manual(values = diagnosis_colors)

  p3 <- insert_xaxis_grob(p_main, x_marginal_boxplot, grid::unit(.1, "null"), position = "top")
  p4 <- insert_yaxis_grob(p3, y_marginal_boxplot, grid::unit(.1, "null"), position = "right")
  p_boxplot <- ggdraw(p4)

  ggsave(
    str_glue(
      "./plots/plot_Fig2to5_panelBtoE/cell_count_proportion_",
      "{x_to_plot}_vs_{y_to_plot}_in_{diagnosis_1}_vs_{diagnosis_2}_",
      "scatter_plus_average_per_sample_add_marginal_boxplot.pdf"
    ),
    plot = p_boxplot,
    device = cairo_pdf(),
    width = 4,
    height = 3,
    useDingbats = FALSE
  )
  dev.off()

  # scatter of the mean by sample, add error bars to show SEM from the replicates
  p <- df_summarized_to_plot %>%
    ggplot(
      aes(
        x = .data[[paste("proportion_mean", x_to_plot, sep = "_")]],
        y = .data[[paste("proportion_mean", y_to_plot, sep = "_")]],
      )
    )
  p_SEM <- p +
    geom_errorbar(
      aes(
        ymin = .data[[paste("proportion_mean", y_to_plot, sep = "_")]] - .data[[paste("proportion_sem", y_to_plot, sep = "_")]],
        ymax = .data[[paste("proportion_mean", y_to_plot, sep = "_")]] + .data[[paste("proportion_sem", y_to_plot, sep = "_")]],
        color = diagnosis
      ),
      width = 0
    ) +
    geom_errorbarh(
      aes(
        xmin = .data[[paste("proportion_mean", x_to_plot, sep = "_")]] - .data[[paste("proportion_sem", x_to_plot, sep = "_")]],
        xmax = .data[[paste("proportion_mean", x_to_plot, sep = "_")]] + .data[[paste("proportion_sem", x_to_plot, sep = "_")]],
        color = diagnosis
      ),
      height = 0
    ) +
    geom_point(
      aes(fill = diagnosis),
      color = "black",
      shape = 21,
      size = 2
    ) +
    xlab(str_glue("Proportion of cells expressing {x_gene} at {x_timepoint}")) +
    ylab(str_glue("Proportion of cells expressing {y_gene} at {y_timepoint}")) +
    scale_fill_manual(values = diagnosis_colors) +
    scale_color_manual(values = diagnosis_colors) +
    theme_bw(base_size = 8, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "lines")
    )

  ggsave(
    str_glue(
      "./plots/plot_Fig2to5_panelBtoE/cell_count_proportion_",
      "{x_to_plot}_vs_{y_to_plot}_in_{diagnosis_1}_vs_{diagnosis_2}_",
      "average_per_sample_and_SEM.pdf"
    ),
    plot = p_SEM,
    device = cairo_pdf(),
    width = 4,
    height = 3,
    useDingbats = FALSE
  )
  dev.off()

  x_SEM_marginal_boxplot <- axis_canvas(p_SEM, axis = "x") +
    geom_boxplot(
      data = df_summarized_to_plot,
      aes(
        x = .data[[paste("proportion_mean", x_to_plot, sep = "_")]],
        color = diagnosis
      ),
      outlier.shape = NA,
      position = position_dodge(width = 0.5),
      width = 0.4,
      size = 0.5
    ) +
    annotate(
      "text",
      x = stat_test_x$y.position,
      y = 0,
      hjust = -0.1,
      label = str_glue("p = {stat_test_x$p}"),
      size = 1.5
    ) +
    annotate(
      "segment",
      x = stat_test_x$y.position,
      xend = stat_test_x$y.position,
      y = 0.125,
      yend = -0.125,
      size = 0.25
    ) +
    annotate(
      "segment",
      x = stat_test_x$y.position,
      xend = stat_test_x$y.position - 0.01,
      y = 0.125,
      yend = 0.125,
      size = 0.25
    ) +
    annotate(
      "segment",
      x = stat_test_x$y.position,
      xend = stat_test_x$y.position - 0.01,
      y = -0.125,
      yend = -0.125,
      size = 0.25
    ) +
    scale_color_manual(values = diagnosis_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

  y_SEM_marginal_boxplot <- axis_canvas(p_SEM, axis = "y", coord_flip = TRUE) +
    geom_boxplot(
      data = df_summarized_to_plot,
      aes(
        x = .data[[paste("proportion_mean", y_to_plot, sep = "_")]],
        color = diagnosis
      ),
      outlier.shape = NA,
      position = position_dodge(width = 0.5),
      width = 0.4,
      size = 0.5
    ) +
    annotate(
      "text",
      x = stat_test_y$y.position,
      y = 0,
      hjust = -0.1,
      label = str_glue("p = {stat_test_y$p}"),
      size = 1.5
    ) +
    annotate(
      "segment",
      x = stat_test_y$y.position,
      xend = stat_test_y$y.position,
      y = 0.125,
      yend = -0.125,
      size = 0.25
    ) +
    annotate(
      "segment",
      x = stat_test_y$y.position,
      xend = stat_test_y$y.position - 0.01,
      y = 0.125,
      yend = 0.125,
      size = 0.25
    ) +
    annotate(
      "segment",
      x = stat_test_y$y.position,
      xend = stat_test_y$y.position - 0.01,
      y = -0.125,
      yend = -0.125,
      size = 0.25
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    coord_flip() +
    scale_color_manual(values = diagnosis_colors)

  p5 <- insert_xaxis_grob(p_SEM, x_SEM_marginal_boxplot, grid::unit(.1, "null"), position = "top")
  p6 <- insert_yaxis_grob(p5, y_SEM_marginal_boxplot, grid::unit(.1, "null"), position = "right")
  p_SEM_boxplot <- ggdraw(p6)

  ggsave(
    str_glue(
      "./plots/plot_Fig2to5_panelBtoE/cell_count_proportion_",
      "{x_to_plot}_vs_{y_to_plot}_in_{diagnosis_1}_vs_{diagnosis_2}_",
      "average_per_sample_and_SEM_add_maginal_boxplot.pdf"
    ),
    plot = p_SEM_boxplot,
    device = cairo_pdf(),
    width = 4,
    height = 3,
    useDingbats = FALSE
  )
  dev.off()

  stat_res <- bind_rows(stat_test_x, stat_test_y) %>%
    dplyr::select(
      gene, timepoint, .y.,
      method, alternative, group1, group2, n1, n2,
      estimate1, estimate2, estimate, df, statistic, p, conf.low, conf.high
    )

  stat_res
}

# run through all combinations of genes and timepoints of interest
params <- tibble(
  x_gene = rep(c("SOX2", "KI67", "CC3", "CTIP2"), 4),
  x_timepoint = rep(c("DIV28", "DIV28", "DIV28", "DIV52"), 4),
  y_gene = rep(c("SOX2", "KI67", "CC3", "TBR2"), 4),
  y_timepoint = rep("DIV52", 16),
  diagnosis_1 = rep(c("MIC", "PMG", "EPI", "ID"), each = 4),
  diagnosis_2 = "Ctrl"
)

merged_stat_res <- params %>%
  pmap_dfr(run_plot_and_test)

write_tsv(merged_stat_res, "./statistics_test_results_Fig2to5_panelBtoE.tsv")

## log sessionInfo
sessionInfo()
