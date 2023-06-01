# Load libraries
source("src/libraries.R")

# Load Unrarefied data
load("data/processed_data/final_physeq_met4_cnts.RData")
final_physeq_met4_cnts


# Metadata
df_metadata <- data.frame(sample_data(final_physeq_met4_cnts)) %>%
    rownames_to_column(var = "sample_id") %>%
    mutate(visit = factor(visit, levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT")))

# Taxa palette family level -----
myPal_family <- tax_palette(final_physeq_met4_cnts, rank = "Family", pal = "brewerPlus", n = 12)

# Flip the color for Rumminococcaceae to red
myPal_family["f__Ruminococcaceae"] <- "#E31A1C"
myPal_family["f__Clostridia_unclassified"] <- "#B2DF8A"


# Figure 3.a ----
df_rumino_family <- final_physeq_met4_cnts %>%
    speedyseq::tax_glom(taxrank = "Family") %>%
    speedyseq::psmelt() %>%
    data.frame() %>%
    dplyr::select(subject, visit, Family, Abundance) %>%
    group_by(subject, visit, Family, .drop = FALSE) %>%
    dplyr::summarize(cnt = sum(Abundance)) %>%
    dplyr::mutate(pct = cnt / sum(cnt) * 100) %>%
    inner_join(df_metadata, by = c("subject", "visit")) %>%
    mutate(arm = case_when(
        arm == "SER-401 Active" ~ "SER-401",
        arm == "SER-401 Placebo" ~ "Placebo + Nivo",
    )) %>%
    mutate(arm = factor(arm, levels = c("SER-401", "Placebo + Nivo"))) %>%
    ungroup() %>%
    filter(Family == "f__Ruminococcaceae") %>%
    mutate(visit = factor(visit, levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT"))) %>%
    arrange(visit) %>%
    group_by(subject, arm) %>%
    mutate(pct = pct - pct[visit == "BL"])

df_rumino_family %>%
    ggplot(., aes(x = visit, y = pct, color = arm)) +
    geom_point(size = 3, shape = 21, stroke = 1.5) +
    geom_line(aes(group = subject), alpha = 0.3) +
    stat_summary(
        fun.y = median,
        fun.ymin = median,
        fun.ymax = median,
        geom = "crossbar",
        size = 1,
        width = 0.5
    ) +
    scale_color_manual(values = c("SER-401" = "#F8766D", "Placebo + Nivo" = "#00BFC4")) +
    facet_grid(arm ~ ., switch = "y") +
    theme_bw() +
    ylab("Percent change w.r.t BL") +
    theme(text = element_text(size = 16)) +
    labs(color = "Arm") +
    ylim(c(-25, 25)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggtitle("Ruminococcaceae change over time")
ggsave("output/figures/figure_3/Ruminococcaceae_over_time_wrt_Baseline.pdf", width = 12, height = 7)


## Stat for Figure 3.a between BL and D-7

df_stat_rumino_family <- NULL
for (i in c("D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT")) {
    for (j in c("Placebo + Nivo", "SER-401")) {
        df_res <- df_rumino_family %>%
            filter(arm == j) %>%
            filter(visit %in% c("BL", i)) %>%
            select(subject, visit, pct) %>%
            ungroup() %>%
            pivot_wider(names_from = visit, values_from = pct) %>%
            do(tidy(wilcox.test(.$BL, .data[[i]], data = ., paired = TRUE))) %>%
            mutate(arm = j, timepoint = i) %>%
            dplyr::select(arm, timepoint, everything())
        df_stat_rumino_family <- rbind(df_stat_rumino_family, df_res)
    }
}
# write the results to csv
df_stat_rumino_family %>%
    mutate(FDR_p_value = p.adjust(p.value, method = "fdr")) %>%
    write_csv(., file = "output/figures/figure_3/df_stat_rumino_family.csv")





# Figure 3.b ----
# Ruminococcaceae with respect to baseline ----
final_physeq_met4_cnts %>%
    speedyseq::tax_glom(taxrank = "Family") %>%
    speedyseq::psmelt() %>%
    data.frame() %>%
    dplyr::select(subject, visit, Family, Abundance) %>%
    group_by(subject, visit, Family, .drop = FALSE) %>%
    dplyr::summarize(cnt = sum(Abundance)) %>%
    dplyr::mutate(pct = cnt / sum(cnt) * 100) %>%
    inner_join(df_metadata, by = c("subject", "visit")) %>%
    mutate(arm = case_when(
        arm == "SER-401 Active" ~ "SER-401",
        arm == "SER-401 Placebo" ~ "Placebo + Nivo",
    )) %>%
    mutate(arm = factor(arm, levels = c("SER-401", "Placebo + Nivo"))) %>%
    ungroup() %>%
    filter(Family == "f__Ruminococcaceae") %>%
    mutate(visit = factor(visit, levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT"))) %>%
    arrange(visit) %>%
    group_by(subject, arm) %>%
    mutate(pct = pct - pct[visit == "BL"]) %>%
    mutate(best_overall_response.y = factor(best_overall_response.y, levels = c(
        "Complete Response",
        "Partial Response",
        "Stable Disease",
        "Progressive Disease"
    ))) %>%
    mutate(rumino_scrn = factor(rumino_scrn, levels = c(
        "Low",
        "High"
    ))) %>%
    ggplot(., aes(x = visit, y = pct, color = best_overall_response.y)) +
    geom_point(size = 3) +
    geom_line(aes(group = subject), alpha = 0.8, size = 2) +
    scale_color_manual(values = c(
        "Stable Disease" = "#fdae61",
        "Progressive Disease" = "#d7191c",
        "Partial Response" = "#a6d96a",
        "Complete Response" = "#1a9850"
    )) +
    labs(color = "Best overall response") +
    ggrepel::geom_text_repel(
        data = . %>% group_by(subject) %>% filter(as.numeric(visit) == max(as.numeric(visit))),
        aes(label = as.character(subject), color = best_overall_response.y),
        show.legend = FALSE,
        max.overlaps = 5
    ) +
    ggnewscale::new_scale_color() +
    stat_summary(
        aes(x = visit, y = median_val, color = arm),
        data = . %>% group_by(visit, arm, rumino_scrn) %>% summarize(median_val = median(pct)),
        fun = "median",
        geom = "crossbar",
        width = 0.5,
        alpha = 0.9, inherit.aes = FALSE
    ) +
    scale_color_manual(values = c("SER-401" = "#F8766D", "Placebo + Nivo" = "#00BFC4")) +
    facet_grid(arm ~ rumino_scrn, switch = "y") +
    theme_bw() +
    ylab("Percent change w.r.t BL") +
    theme(text = element_text(size = 16)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(color = "Arm") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggtitle("Ruminococcaceae change over time \n by Rumino screening status") +
    labs(caption = "")
ggsave("output/figures/figure_3/Ruminococcaceae_over_time_by_arm_and_rumino_screen.pdf", width = 14, height = 7)


# Figure 3.c ----
p_trt_arm <- final_physeq_met4_cnts %>%
    speedyseq::mutate_sample_data(arm = case_when(
        arm == "SER-401 Active" ~ "SER-401",
        arm == "SER-401 Placebo" ~ "Placebo + Nivo",
    )) %>%
    speedyseq::mutate_sample_data(
        subject = factor(
            subject,
            levels = c(
                "113-S0001",
                "107-S0002",
                "102-S0001",
                "102-S0002",
                "102-S0006",
                "113-S0002",
                "108-S0002",
                "102-S0005"
            )
        )
    ) %>%
    speedyseq::mutate_sample_data(arm = factor(arm, levels = c("SER-401", "Placebo + Nivo"))) %>%
    speedyseq::filter_sample_data(arm == "SER-401") %>%
    speedyseq::mutate_sample_data(visit = factor(
        visit,
        levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT")
    )) %>%
    ps_arrange(arm, visit) %>%
    tax_fix(unknowns = c("g__")) %>%
    tax_fix(unknowns = c("unclassified")) %>%
    comp_barplot(tax_level = "Family", label = "visit", sample_order = "asis", n_taxa = 10, palette = myPal_family) +
    facet_wrap(subject ~ ., scales = "free", switch = "y", ncol = 1) +
    ggtitle("Top 10 Family") +
    geom_text(
        data = . %>% filter(
            visit == "D-7"
        ) %>%
            group_by(subject) %>%
            filter(Abundance == max(Abundance)),
        aes(
            label = visit, fontface = "bold"
        ), size = 5, hjust = 0, vjust = 0, nudge_y = -0.35, nudge_x = -0.3
    ) +
    theme(axis.text.x = element_blank()) +
    theme(legend.position = "none") +
    scale_y_continuous(position = "right")


p_trt_placebo <- final_physeq_met4_cnts %>%
    speedyseq::mutate_sample_data(arm = case_when(
        arm == "SER-401 Active" ~ "SER-401",
        arm == "SER-401 Placebo" ~ "Placebo + Nivo",
    )) %>%
    speedyseq::mutate_sample_data(
        subject = factor(
            subject,
            levels = c(
                "102-S0003",
                "102-S0009",
                "107-S0003",
                "107-S0006",
                "108-S0004",
                "112-S0001"
            )
        )
    ) %>%
    speedyseq::mutate_sample_data(arm = factor(arm, levels = c("SER-401", "Placebo + Nivo"))) %>%
    speedyseq::filter_sample_data(arm == "Placebo + Nivo") %>%
    speedyseq::mutate_sample_data(visit = factor(
        visit,
        levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT")
    )) %>%
    ps_arrange(arm, visit) %>%
    tax_fix(unknowns = c("g__")) %>%
    tax_fix(unknowns = c("unclassified")) %>%
    comp_barplot(tax_level = "Family", label = "visit", sample_order = "asis", n_taxa = 10, palette = myPal_family) +
    facet_wrap(subject ~ ., scales = "free", switch = "y", ncol = 1) +
    ggtitle("Top 10 Family") +
    theme(axis.text.x = element_blank()) +
    scale_y_continuous(position = "right")


# Combine the two plots using patchwork package
p_trt_arm + p_trt_placebo
ggsave("output/figures/figure_3/Top_10_Family_by_arm.pdf", width = 12, height = 10)
