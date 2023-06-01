# Load libraries
source("src/libraries.R")


# Extended data Figure 5.a ----

df_fav <- read_csv("data/processed_data/MCGRAW_Humann_CPM_TargetedPathways_Beneficial.txt")

df_fav %>%
    mutate(Visit = factor(Visit, levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1"))) %>%
    filter(!is.na(Visit)) %>%
    filter(!is.na(Treatment.Name)) %>%
    group_by(Visit, Pathway_Category) %>%
    do(broom::tidy(wilcox.test(Value ~ Treatment.Name, data = .))) %>%
    data.frame() %>%
    mutate(adj_p_value = p.adjust(p.value, method = "fdr")) %>%
    mutate(adj_p_value = format.pval(adj_p_value, digits = 2)) %>%
    mutate(p.value = format.pval(p.value, digits = 2)) %>%
    gt::gt() %>%
    tab_style(
        style = list(
            cell_fill(color = "#e41a1d2c")
        ),
        locations = cells_body(rows = p.value < 0.05)
    ) %>%
    tab_header(
        title = md(paste0("Favorable Pathways  "))
    ) %>%
    gtsave("output/figures/figure_5_ext_data/wilcox_Test_Favorable_Pathways_Table.html")



df_fav %>%
    mutate(Visit = factor(Visit, levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1"))) %>%
    filter(!is.na(Visit)) %>%
    filter(!is.na(Treatment.Name)) %>%
    mutate(Treatment.Name = factor(Treatment.Name, levels = c("SER-401", "Placebo"))) %>%
    ggplot(., aes(x = Visit, y = Value, color = Treatment.Name)) +
    geom_point(alpha = 0.4) +
    geom_line(aes(group = Subject), alpha = 0.1) +
    stat_summary(
        fun = median, fun.min = median, fun.max = median,
        geom = "crossbar", width = 0.5
    ) +
    facet_wrap(. ~ Pathway_Category, ncol = 4, scales = "free_y") +
    scale_y_log10() +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave("output/figures/figure_5_ext_data/Line_plot_Favorable_Pathways.pdf", width = 20, height = 10)



# Extended data Figure 5.b ----
# Unfavorable pathways of interest ----

df_unfav <- read_csv("data/processed_data/MCGRAW_Humann_CPM_TargetedPathways_Adverse.txt")

df_unfav %>%
    filter(!is.na(Visit)) %>%
    dplyr::select(Category, Visit, Value, Treatment.Name) %>%
    mutate(Visit = factor(Visit, levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1"))) %>%
    filter(!is.na(Visit)) %>%
    filter(!is.na(Treatment.Name)) %>%
    group_by(Visit, Category) %>%
    do(broom::tidy(wilcox.test(Value ~ Treatment.Name, data = .))) %>%
    data.frame() %>%
    mutate(adj_p_value = p.adjust(p.value, method = "fdr")) %>%
    mutate(adj_p_value = format.pval(adj_p_value, digits = 2)) %>%
    mutate(p.value = format.pval(p.value, digits = 2)) %>%
    gt::gt() %>%
    tab_style(
        style = list(
            cell_fill(color = "#e41a1d2c")
        ),
        locations = cells_body(rows = p.value < 0.05)
    ) %>%
    tab_header(
        title = md(paste0("Unfavorable Pathways  "))
    ) %>%
    gtsave("output/figures/figure_5_ext_data/wilcox_Test_UNFavorable_Pathways_Table_grouped.html")



df_unfav %>%
    mutate(Visit = factor(Visit, levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1"))) %>%
    filter(!is.na(Visit)) %>%
    filter(!is.na(Treatment.Name)) %>%
    mutate(Treatment.Name = factor(Treatment.Name, levels = c("SER-401", "Placebo"))) %>%
    ggplot(., aes(x = Visit, y = Value, color = Treatment.Name)) +
    geom_point(alpha = 0.4) +
    geom_line(aes(group = Subject), alpha = 0.1) +
    stat_summary(
        fun = median, fun.min = median, fun.max = median,
        geom = "crossbar", width = 0.5
    ) +
    facet_wrap(. ~ Category, ncol = 4, scales = "free_y") +
    scale_y_log10() +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave("output/figures/figure_5_ext_data/lineplot_Unfavorable_Pathways.pdf", width = 20, height = 10)
