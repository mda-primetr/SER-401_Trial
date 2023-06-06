# Load libraries
source("src/libraries.R")


# Set the seed for reproducibility
RNGkind(sample.kind = "Rounding")
set.seed(711)

# Load Unrarefied data
load("data/processed_data/final_physeq_met4_cnts.RData")
final_physeq_met4_cnts

tax_table_asv_species <- data.frame(tax_table(final_physeq_met4_cnts)) %>%
    rownames_to_column(var = "tax_id") %>%
    dplyr::select(1, Species) %>%
    view()

# Save taxa table for network analysis metadata such as taxonomy and abundance
tax_table_for_network <- data.frame(tax_table(final_physeq_met4_cnts)) %>%
    mutate(Species = gsub("s__", "", Species)) %>%
    distinct(Species, .keep_all = TRUE) %>%
    rownames_to_column(var = "tax_id")



# Fig 4.a

# **This code will take some time to run so please be patient and do not interrupt the code**
# You are welcome to skip the execution of this code and use the network output from data/processed_data folder

# Plot network analysis figures for each arm and timepoint
for (ARM in c("SER-401 Placebo", "SER-401 Active")) {
    for (tp in c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1")) {
        load(file = paste0("data/processed_data/network_out_", ARM, "_", tp, ".RData"))

        # That load file was stored as "secom_linear" variable

        # Save the results for plotting degree density plot which is not filtered for co-occurrence
        corr_dist_degree <- secom_linear$corr_fl


        # Store the secom sparse correlation distance matrix
        corr_dist <- secom_linear$corr_fl
        # Store the secom cooccurance matrix
        cooccur_dist <- secom_linear$mat_coocc


        # Filter by co-occurrence
        ## Get the Ns per arm and timepoint
        # This code filters the phyloseq object final_physeq_met4_cnts by the ARM and timepoint
        # and then saves the number of samples in overlap_N.
        # The code is used to ensure that there are enough samples for the network analysis

        final_physeq_met4_cnts_n <- final_physeq_met4_cnts %>%
            speedyseq::filter_sample_data(arm == ARM) %>%
            speedyseq::filter_sample_data(visit == tp)
        overlap_N <- round(NROW(sample_data(final_physeq_met4_cnts_n)) * 0.60)

        corr_dist[cooccur_dist < overlap_N] <- 0


        # make a dataframe from the correlation matrix
        df1 <- corr_dist %>%
            data.frame() %>%
            rownames_to_column(var = "term") %>%
            pivot_longer(-term) %>%
            inner_join(tax_table_asv_species, by = c("term" = "tax_id")) %>%
            mutate(name = gsub("X", "", name)) %>%
            inner_join(tax_table_asv_species, by = c("name" = "tax_id")) %>%
            dplyr::select(-c(term, name)) %>%
            dplyr::rename(term = Species.x, name = Species.y) %>%
            # keep only the upper triangle of the correlation matrix
            filter(term > name) %>%
            # keep only correlations that are not zero
            filter(abs(value) > 0) %>%
            # remove the s__ prefix from the terms
            dplyr::mutate(term = gsub("s__", "", term)) %>%
            dplyr::mutate(name = gsub("s__", "", name)) %>%
            dplyr::select(term, name, value) %>%
            filter(!is.na(value))

        # make a graph from the dataframe
        g1 <- as_tbl_graph(df1) %>%
            dplyr::mutate(degree_measure = centrality_degree()) %>%
            dplyr::mutate(centrality = centrality_pagerank(weight = value)) %>%
            # convert the graph to an undirected graph
            to_undirected() %>%
            activate(edges) %>%
            # add a sign column to the edges
            dplyr::mutate(cor_sign = case_when(
                value > 0 ~ "P",
                value < 0 ~ "N"
            )) %>%
            activate(nodes) %>%
            dplyr::mutate(community = as.factor(group_louvain()))

        data <- toVisNetworkData(g1)

        if (ARM == "SER-401 Placebo") {
            data$edges$color[data$edges$cor_sign == "N"] <- "#f8766da1"
            data$nodes$color <- "#00bec4b0"
            bar_color <- "#00BFC4"
        }

        if (ARM == "SER-401 Active") {
            data$edges$color[data$edges$cor_sign == "N"] <- "#67a9cf"
            data$nodes$color <- "#f8766da1"
            bar_color <- "#F8766D"
        }

        visNet <- visNetwork(
            nodes = data$nodes,
            edges = data$edges,
            width = "100%",
            height = "100vh",
            main = (title <- paste0(ARM, " ", tp, " network")),
        ) %>%
            visOptions(
                selectedBy = "community",
                highlightNearest = TRUE,
                nodesIdSelection = TRUE
            ) %>%
            visPhysics(solver = "repulsion", stabilization = FALSE) # to steady the network graphic when opening the html file
        visSave(visNet, file = paste0("output/figures/figure_4_ext_data/networks/visNet_", ARM, "_", tp, ".html"), selfcontained = TRUE)



        # This code is used to calculate the Page Rank centrality of the species in the
        # network. It is used to determine which species are most central to the network
        # and thus which species are most likely to be important for the study.
        #
        # The code does this by first creating a graph object from the network, then
        # calculating the Page Rank centrality for each node in the network. The code
        # then creates a data frame from the graph object, and plots the data frame using
        # ggplot. The output is a bar chart showing the Page Rank centrality for each
        # species.
        #
        # The code uses the following identifiers:
        # g1: the graph object created from the network
        # centrality: the vector containing the Page Rank centrality for each node
        # name: the species name for each node
        # community: the community to which each node belongs
        # ARM: the ARM to which the network belongs
        # tp: the time point to which the network belongs

        # Plot Page Rank centrality and communities as shown in that visnetwork figure
        # *This plot is filtered for co-occurrence *
        g1 %>%
            activate(nodes) %>%
            as.data.frame() %>%
            ggplot(., aes(x = reorder(name, centrality), y = centrality)) +
            geom_bar(stat = "identity", fill = bar_color) +
            coord_flip() +
            xlab("Species") +
            ylab("Page Rank centrality") +
            facet_wrap(community ~ ., scales = "free_y") +
            theme_bw() +
            theme(text = element_text(size = 16)) +
            theme(legend.position = "none") +
            ggtitle(paste0(ARM, "-", tp, " Network centrality and communities"))
        ggsave(paste0(
            "output/figures/figure_4_ext_data/networks/Page_rank_centrality_",
            ARM,
            "_",
            tp,
            ".pdf"
        ), width = 25, height = 20, units = "in", dpi = 300)




        #   *This plot DOES NOT filter for co-occurrence*
        df2 <- corr_dist_degree %>%
            data.frame() %>%
            rownames_to_column(var = "term") %>%
            pivot_longer(-term) %>%
            inner_join(tax_table_asv_species, by = c("term" = "tax_id")) %>%
            mutate(name = gsub("X", "", name)) %>%
            inner_join(tax_table_asv_species, by = c("name" = "tax_id")) %>%
            dplyr::select(-c(term, name)) %>%
            dplyr::rename(term = Species.x, name = Species.y) %>%
            dplyr::mutate(term = gsub("s__", "", term)) %>%
            dplyr::mutate(name = gsub("s__", "", name)) %>%
            dplyr::select(term, name, value) %>%
            filter(!is.na(value))


        g2 <- as_tbl_graph(df2) %>%
            to_undirected() %>%
            dplyr::mutate(degree_measure = centrality_degree()) %>%
            dplyr::mutate(centrality = centrality_pagerank(weight = value)) %>%
            to_undirected() %>%
            activate(edges) %>%
            dplyr::mutate(cor_sign = case_when(
                value > 0 ~ "P",
                value < 0 ~ "N"
            ))

        # This code creates a dataframe of the degree distribution for a given graph, and plots it with the fraction of nodes with that degree on the y-axis and the degree on the x-axis.
        # The graph is created by creating a new graph object, adding the nodes and edges, and then activating the nodes so that the functions can be applied to the nodes.
        # The dataframe is created by converting the graph object to a dataframe, ordering the dataframe by the degree measure, and then grouping by the degree measure and summarizing the number of nodes with that degree.
        # The dataframe is then plotted with the fraction of nodes with that degree on the y-axis and the degree on the x-axis.



        df_g3 <- g2 %>%
            activate(nodes) %>%
            as.data.frame() %>%
            arrange(desc(degree_measure))

        # Plot fraction of degree measure and number of nodes
        df_g3 %>%
            inner_join(tax_table_for_network, by = c("name" = "Species")) %>%
            write_csv(paste0(
                "output/figures/figure_4_ext_data/networks/degree_nodes_connections_",
                ARM,
                "_",
                tp,
                ".csv"
            )) %>%
            group_by(degree_measure) %>%
            ungroup() %>%
            ggplot(., aes(x = degree_measure, y = centrality)) +
            geom_point() +
            geom_hex(bins = 25) +
            scale_fill_continuous(type = "viridis") +
            scale_x_log10() +
            theme_bw() +
            geom_label_repel(
                data = . %>% filter(dense_rank(desc(degree_measure)) <= 10), aes(label = name),
                size = 3,
                nudge_x = 0.1, nudge_y = 0.01,
                segment.size = 0.2, segment.color = "grey50", point.padding = 0.2, force = 10, max.overlaps = Inf, segment.alpha = 0.5
            ) +
            theme(text = element_text(size = 18)) +
            ylab("Centrality") +
            xlab("Degree: Number of connections") +
            ggtitle(paste0(ARM, "-", tp, "Degree and centrality")) +
            labs(caption = "Labelled taxa are the top 10 most connected taxa in the network")
        ggsave(paste0(
            "output/figures/figure_4_ext_data/networks/degree_dist_",
            ARM,
            "_",
            tp,
            ".pdf"
        ), width = 10, height = 10, units = "in", dpi = 300)
    }
}
