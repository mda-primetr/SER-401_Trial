# Load libraries
source("src/libraries.R")

# Set the seed for reproducibility
RNGkind(sample.kind = "Rounding")
set.seed(711)

# Load Unrarefied data
load("data/processed_data/final_physeq_met4_cnts.RData")
final_physeq_met4_cnts

sample_data(final_physeq_met4_cnts) %>% view()

# For each arm and visit
for (ARM in c("SER-401 Placebo", "SER-401 Active")) {
    for (tp in c("BL", "D-7", "C1D1", "C1D8", "C2D1", "C3D1")) {
        # Subset the phyloseq object to the arm and visit and rename the samples
        final_phyloseq_for_netwrk <- final_physeq_met4_cnts %>%
            speedyseq::filter_sample_data(arm == ARM) %>%
            # speedyseq::mutate_sample_data(sam_acc = paste0("S", acc_number)) %>%
            speedyseq::filter_sample_data(visit == tp)

        sample_names(final_phyloseq_for_netwrk) <- data.frame(sample_data(final_phyloseq_for_netwrk))$subject

        # Create a secom_linear object
        secom_linear <- secom_linear(
            pseqs = list(c(final_phyloseq_for_netwrk, final_phyloseq_for_netwrk)),
            method = c("spearman"),
            prv_cut = 0.20, lib_cut = 0, n_cl = 16, max_p = 0.05, thresh_hard = 0.70
        )

        # Save the secom_linear object
        save(secom_linear, file = paste0("data/processed_data/network_out_", ARM, "_", tp, ".RData"))
    }
}
