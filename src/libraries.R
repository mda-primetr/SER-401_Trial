# Load the libraries for the analysis ----
library(ANCOMBC)
library(broom)
library(ggarchery)
library(ggprism)
library(ggraph)
library(ggrepel)
library(ggtext)
library(gt)
library(metagMisc)
library(microViz)
library(patchwork)
library(phyloseq)
library(rstatix)
library(scales)
library(speedyseq)
library(tidygraph)
library(tidyverse)
library(visNetwork)


# Here is the session info ----

# > sessionInfo()
# R version 4.3.0 (2023-04-21)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.4

# Matrix products: default
# BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

# Random number generation:
#  RNG:     Mersenne-Twister
#  Normal:  Inversion
#  Sample:  Rounding

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# time zone: America/Chicago
# tzcode source: internal

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
#  [1] ggarchery_0.4.2      patchwork_1.1.2      broom_1.0.4
#  [4] visNetwork_2.1.2     tidygraph_1.2.3      ggraph_2.1.0
#  [7] metagMisc_0.5.0      ANCOMBC_1.6.2        ggtext_0.1.2
# [10] ggrepel_0.9.3        scales_1.2.1         rstatix_0.7.2
# [13] ggprism_1.0.4        microViz_0.10.8      lubridate_1.9.2
# [16] forcats_1.0.0        stringr_1.5.0        dplyr_1.1.2
# [19] purrr_1.0.1          readr_2.1.4          tidyr_1.3.0
# [22] tibble_3.2.1         ggplot2_3.4.2        tidyverse_2.0.0
# [25] speedyseq_0.5.3.9018 phyloseq_1.44.0
