# Load libraries
source("src/libraries.R")

# Load Unrarefied data
load("data/processed_data/final_physeq_met4_cnts.RData")
final_physeq_met4_cnts

# Metadata
df_metadata <- data.frame(sample_data(final_physeq_met4_cnts)) %>%
    rownames_to_column(var = "sample_id") %>%
    mutate(visit = factor(visit, levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT")))


# All those spore forming families combined -----

list_of_spore_formers <- c(
    "Ruminococcaceae",
    "Lachnospiraceae",
    "Erysipelotrichaceae",
    "Clostridiaceae",
    "Eubacteriaceae",
    "Peptostreptococcaceae"
)


# Figure 3.a ----
final_physeq_met4_cnts %>%
    speedyseq::tax_glom(taxrank = "Family") %>%
    speedyseq::psmelt() %>%
    data.frame() %>%
    dplyr::select(subject, visit, Family, Abundance) %>%
    mutate(Family = gsub("f__", "", Family)) %>%
    mutate(Family = case_when(
        Family %in% list_of_spore_formers ~ "Spore Formers",
        TRUE ~ Family
    )) %>%
    group_by(subject, visit, Family, .drop = FALSE) %>%
    dplyr::summarize(cnt = sum(Abundance)) %>%
    dplyr::mutate(pct = cnt / sum(cnt) * 100) %>%
    inner_join(df_metadata, by = c("subject", "visit")) %>%
    filter(Family == "Spore Formers") %T>%
    write_csv(paste0("data/processed_data/", "spore_formers", "_family_counts.csv")) %>%
    mutate(visit = factor(visit, levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT"))) %>%
    ggplot(., aes(x = visit, y = pct, color = arm)) +
    geom_point(size = 3) +
    geom_line(aes(group = subject), alpha = 0.3) +
    stat_summary(
        fun.y = median,
        fun.ymin = median,
        fun.ymax = median,
        geom = "crossbar",
        size = 1,
        width = 0.5
    ) +
    scale_color_manual(values = c("SER-401 Active" = "#F8766D", "SER-401 Placebo" = "#00BFC4")) +
    facet_grid(arm ~ ., switch = "y") +
    theme_bw() +
    theme(text = element_text(size = 16)) +
    labs(color = "Arm") +
    ggtitle(paste0("Proportion of ", "spore_formers", " over time"))
ggsave(paste0("output/figures/figure_3_ext_data/", "spore_formers", "_over_time.pdf"), width = 12, height = 7)


final_physeq_met4_cnts %>%
    speedyseq::tax_glom(taxrank = "Family") %>%
    speedyseq::psmelt() %>%
    data.frame() %>%
    dplyr::select(subject, visit, Family, Abundance) %>%
    mutate(Family = gsub("f__", "", Family)) %>%
    mutate(Family = case_when(
        Family %in% list_of_spore_formers ~ "Spore Formers",
        TRUE ~ Family
    )) %>%
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
    filter(Family == "Spore Formers") %>%
    mutate(visit = factor(visit, levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT"))) %>%
    arrange(visit) %>%
    group_by(subject, arm) %>%
    mutate(pct = pct - pct[visit == "BL"]) %>%
    ggplot(., aes(x = visit, y = pct, color = arm)) +
    geom_point(size = 3) +
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
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggtitle(paste0("Spore formers (Family Level)", " change over time"))
ggsave(paste0(
    "output/figures/figure_3_ext_data/",
    "spore_formers",
    "_over_time_wrt_Baseline.pdf"
), width = 12, height = 7)



# Figure 3.b ----
# Ruminococcaceae over time ----
final_physeq_met4_cnts %>%
    speedyseq::tax_glom(taxrank = "Family") %>%
    speedyseq::psmelt() %>%
    data.frame() %>%
    dplyr::select(subject, visit, Family, Abundance) %>%
    group_by(subject, visit, Family, .drop = FALSE) %>%
    dplyr::summarize(cnt = sum(Abundance)) %>%
    dplyr::mutate(pct = cnt / sum(cnt) * 100) %>%
    inner_join(df_metadata, by = c("subject", "visit")) %>%
    filter(Family == "f__Ruminococcaceae") %T>%
    write_csv("data/processed_data/Ruminococcaceae_family_counts.csv") %>%
    mutate(visit = factor(visit, levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT"))) %>%
    mutate(best_overall_response.y = factor(best_overall_response.y, levels = c(
        "Complete Response",
        "Partial Response",
        "Stable Disease",
        "Progressive Disease"
    ))) %>%
    ungroup() %>%
    mutate(rumino_scrn = factor(rumino_scrn, levels = c(
        "Low",
        "High"
    ))) %>%
    ggplot(., aes(x = visit, y = pct, color = best_overall_response.y)) +
    geom_point(size = 3) +
    geom_line(aes(group = subject), alpha = 0.3) +
    stat_summary(
        fun.y = median,
        fun.ymin = median,
        fun.ymax = median,
        geom = "crossbar",
        size = 1,
        width = 0.5
    ) +
    # geom_text_repel(
    #     aes(label = subject),
    #     size = 3,
    #     point.padding = unit(0.5, "lines"),
    #     segment.color = "grey50"
    # ) +
    scale_color_manual(values = c("SER-401" = "#F8766D", "Placebo + Nivo" = "#00BFC4")) +
    scale_color_manual(values = c(
        "Stable Disease" = "#fdae61",
        "Progressive Disease" = "#d7191c",
        "Partial Response" = "#a6d96a",
        "Complete Response" = "#1a9850"
    )) +
    facet_grid(arm ~ rumino_scrn, switch = "y") +
    theme_bw() +
    theme(text = element_text(size = 16)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(color = "Arm") +
    ggtitle("Proportion of Ruminococcaceae over time")
ggsave("output/figures/figure_3_ext_data/Ruminococcaceae_over_time.pdf", width = 12, height = 7)


# Figure 3.c ----
# Alpha diversity ----
df_meta_richness <- data.frame(estimate_richness(final_physeq_met4_cnts)) %>%
    rownames_to_column(var = "sample_id") %>%
    inner_join(df_metadata, by = c("sample_id")) %>%
    view()



df_meta_richness %>%
    dplyr::select(InvSimpson, visit, arm) %>%
    group_by(arm) %>%
    wilcox_test(InvSimpson ~ visit)


df_meta_richness %>%
    dplyr::select(InvSimpson, visit, arm, subject) %>%
    ggplot(., aes(x = visit, y = InvSimpson, group = subject)) +
    geom_point() +
    geom_line() +
    facet_grid(subject ~ arm)



# Change in InvSimpson Diversity by Arm
df_meta_richness %>%
    mutate(arm = case_when(
        arm == "SER-401 Active" ~ "SER-401",
        arm == "SER-401 Placebo" ~ "Placebo + Nivo",
    )) %>%
    mutate(arm = factor(arm, levels = c("SER-401", "Placebo + Nivo"))) %>%
    dplyr::select(InvSimpson, visit, arm, subject) %>%
    arrange(visit) %>%
    group_by(subject, arm) %>%
    mutate(change = InvSimpson - InvSimpson[visit == "BL"]) %>%
    ggplot(., aes(x = visit, y = change, color = arm)) +
    geom_point(size = 3) +
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
    theme(text = element_text(size = 16)) +
    labs(color = "Arm") +
    ggtitle("Change in InvSimpson Diversity (%)")
ggsave("output/figures/figure_3_ext_data/change_in_InvSimpson_over_time.pdf", height = 7, width = 12)



# Get table of p values for change in InvSimpson Diversity by Arm ----
df_all_richness_invsimpson <- df_meta_richness %>%
    mutate(arm = case_when(
        arm == "SER-401 Active" ~ "SER-401",
        arm == "SER-401 Placebo" ~ "Placebo + Nivo",
    )) %>%
    mutate(arm = factor(arm, levels = c("SER-401", "Placebo + Nivo"))) %>%
    dplyr::select(InvSimpson, visit, arm, subject) %>%
    arrange(visit) %>%
    group_by(subject, arm) %>%
    mutate(change = InvSimpson - InvSimpson[visit == "BL"]) %>%
    ungroup()


df_all_richness <- NULL
for (i in c("D-7", "C1D1", "C1D8", "C2D1")) {
    for (j in c("SER-401", "Placebo + Nivo")) {
        df_all_richness_temp <- df_all_richness_invsimpson %>%
            filter(arm == j) %>%
            filter(visit %in% c("BL", i)) %>%
            dplyr::select(visit, subject, change) %>%
            pivot_wider(names_from = visit, values_from = change) %>%
            do(tidy(wilcox.test(.$BL, .data[[i]], data = ., paired = TRUE))) %>%
            mutate(arm = j, timepoint = i) %>%
            dplyr::select(arm, timepoint, everything()) %>%
            mutate(FDR_pvalue = p.adjust(p.value, method = "fdr"))

        df_all_richness <- rbind(df_all_richness, df_all_richness_temp)
    }
}
write_csv(df_all_richness, "output/figures/figure_3_ext_data/change_in_InvSimpson_over_time_stat.csv")












# Figure 3.d ----
# Plot bray distance by sample -----
dist_final_phyloseq_before_rarefaction <- phyloseq::distance(final_physeq_met4_cnts,
    method = "bray", weighted = FALSE
)


# Get the distance dataframe by sample ----
df_dist_by_sample <- metagMisc::dist2list(dist_final_phyloseq_before_rarefaction) %>%
    inner_join(df_metadata, by = c("col" = "sample_id")) %>%
    dplyr::select(col, row, value, sgpid, subject, visit, arm) %>%
    inner_join(df_metadata, by = c("row" = "sample_id")) %>%
    dplyr::select(
        col,
        row,
        value,
        sgpid.x,
        subject.x,
        visit.x,
        arm.x,
        sgpid.y,
        subject.y,
        visit.y,
        arm.y,
        best_overall_response.y,
        rumino_scrn
    ) %>%
    filter(subject.x == subject.y) %>%
    dplyr::select(
        subject = subject.x,
        sgpid.x,
        value,
        visit.x,
        visit.y,
        arm = arm.x,
        best_overall_response.y,
        rumino_scrn
    ) %>%
    mutate(visit.x = factor(
        visit.x,
        levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT")
    )) %>%
    mutate(visit.y = factor(
        visit.y,
        levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT")
    )) %>%
    mutate(visit1 = case_when(
        as.numeric(visit.x) < as.numeric(visit.y) ~ visit.x,
        TRUE ~ visit.y
    )) %>%
    mutate(visit2 = case_when(
        as.numeric(visit.x) > as.numeric(visit.y) ~ visit.x,
        TRUE ~ visit.y
    )) %>%
    mutate(visit = paste0(visit1, "-", visit2))

# By Timepoint with respect to baseline ----
df_dist_by_sample %>%
    mutate(arm = case_when(
        arm == "SER-401 Active" ~ "SER-401",
        arm == "SER-401 Placebo" ~ "Placebo + Nivo",
    )) %>%
    mutate(arm = factor(arm, levels = c("SER-401", "Placebo + Nivo"))) %>%
    filter(visit %in% c("BL-D-7", "BL-C1D1", "BL-C1D8", "BL-C2D1", "BL-C3D1", "BL-C4D1", "BL-C7D1", "BL-EOT")) %>%
    mutate(visit = factor(
        visit,
        levels = c("BL-D-7", "BL-C1D1", "BL-C1D8", "BL-C2D1", "BL-C3D1", "BL-C4D1", "BL-C7D1", "BL-EOT")
    )) %>%
    dplyr::select(value, arm, visit, subject, best_overall_response.y, rumino_scrn) %>%
    mutate(rumino_scrn = factor(rumino_scrn, levels = c(
        "Low",
        "High"
    ))) %>%
    ggplot(., aes(x = visit, y = value, color = best_overall_response.y)) +
    geom_point() +
    geom_line(aes(group = subject)) +
    facet_grid(arm ~ rumino_scrn) +
    scale_color_manual(values = c("SER-401" = "#F8766D", "Placebo + Nivo" = "#00BFC4")) +
    theme_classic() +
    scale_color_manual(values = c(
        "Stable Disease" = "#fdae61",
        "Progressive Disease" = "#d7191c",
        "Partial Response" = "#a6d96a",
        "Complete Response" = "#1a9850"
    )) +
    ggrepel::geom_text_repel(
        data = . %>% group_by(subject) %>% filter(as.numeric(visit) == max(as.numeric(visit))),
        aes(label = as.character(subject), color = best_overall_response.y),
        show.legend = FALSE,
        max.overlaps = 5
    ) +
    labs(color = "Best overall response") +
    theme(text = element_text(size = 16)) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
    ggtitle("Change in Bray distance by sample between timepoints and Baseline") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave("output/figures/figure_3_ext_data/bray_distance_by_sample_with_Baseline.pdf", width = 15, height = 8)


# Figure 3.e ----
# PCoA plot -----
dist_final_phyloseq_before_rarefaction <- phyloseq::distance(final_physeq_met4_cnts,
    method = "bray", weighted = FALSE
)

mydf <- as.data.frame(as.matrix(dist_final_phyloseq_before_rarefaction))
write.csv(mydf, file = "dist_final_phyloseq_before_rarefaction_mcgraw.csv")



ord_final_physeq_met4_cnts_all <- ordinate(final_physeq_met4_cnts,
    method = "PCoA", distance = dist_final_phyloseq_before_rarefaction
)


df_ord <- plot_ordination(final_physeq_met4_cnts, ord_final_physeq_met4_cnts_all, justDF = TRUE)

# Get Percent explained valued from the ordination
percent_explained <- 100 * ord_final_physeq_met4_cnts_all$values$Eigenvalues / sum(ord_final_physeq_met4_cnts_all$values$Eigenvalues)




df_ord %>%
    mutate(visit = factor(visit, levels = c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1", "C4D1", "C7D1", "EOT"))) %>%
    select(Axis.1, Axis.2, visit, arm, subject) %>%
    mutate(arm = case_when(
        arm == "SER-401 Active" ~ "SER-401",
        arm == "SER-401 Placebo" ~ "Placebo + Nivo",
    )) %>%
    mutate(arm = factor(arm, levels = c("SER-401", "Placebo + Nivo"))) %>%
    rownames_to_column(var = "sample_id") %>%
    filter(sample_id != "SGP033335") %>%
    filter(visit != "EOT") %>%
    group_by(subject) %>%
    arrange(subject, visit) %>%
    dplyr::mutate(Axis.1_end = case_when(
        visit == "BL" ~ lead(Axis.1),
        visit == "D-7" ~ lead(Axis.1),
        visit == "C1D1" ~ lead(Axis.1),
        visit == "C1D8" ~ lead(Axis.1),
        visit == "C2D1" ~ lead(Axis.1),
        visit == "C3D1" ~ lead(Axis.1),
        visit == "C4D1" ~ lead(Axis.1),
        visit == "C7D1" ~ Axis.1
    )) %>%
    dplyr::mutate(Axis.2_end = case_when(
        visit == "BL" ~ lead(Axis.2),
        visit == "D-7" ~ lead(Axis.2),
        visit == "C1D1" ~ lead(Axis.2),
        visit == "C1D8" ~ lead(Axis.2),
        visit == "C2D1" ~ lead(Axis.2),
        visit == "C3D1" ~ lead(Axis.2),
        visit == "C4D1" ~ lead(Axis.2),
        visit == "C7D1" ~ Axis.2
    )) %>%
    ggplot(., aes(
        x = Axis.1, y = Axis.2,
        color = `arm`,
        label = visit,
        group = subject
    )) +
    geom_point(size = 3) +
    geom_arrowsegment(aes(x = Axis.1, xend = Axis.1_end, y = Axis.2, yend = Axis.2_end, fill = arm),
        arrow_positions = 0.5,
        arrows = arrow(),
        size = 1,
        show.legend = FALSE
    ) +
    theme_bw() +
    facet_wrap(arm ~ subject, nrow = 7) +
    xlab(glue::glue("PCoA1: {round(percent_explained[1], digits = 1)}%")) +
    ylab(glue::glue("PCoA2: {round(percent_explained[2], digits = 1)}%")) +
    scale_fill_manual(values = c("SER-401" = "#F8766D", "Placebo + Nivo" = "#00BFC4")) +
    scale_color_manual(values = c("SER-401" = "#F8766D", "Placebo + Nivo" = "#00BFC4")) +
    geom_text_repel(show_guide = FALSE) +
    ggtitle("PCoA using Bray distance") +
    theme(legend.position = "top") +
    theme_bw() +
    labs(color = "Arm") +
    theme(text = element_text(size = 12), axis.text.x = element_text(size = 8))
ggsave("output/figures/figure_3_ext_data/PcoA_by_arm_and_sample_with_direction_arrow.pdf", width = 9, height = 15)
