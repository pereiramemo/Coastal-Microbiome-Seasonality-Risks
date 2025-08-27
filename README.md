# Coastal-Microbiome-Seasonality-Risks

This repository contains the data and scripts to reproduce the main findings shown in the study ["Seasonal Dynamics and Thermal Vulnerability of Marine Microbial Communities in a Coastal South Atlantic Warming Hotspot" by E. Pereira, Zanetti J, Griffero L,  Martínez A, Amann R, Alonso C.](https://doi.org/10.1101/2025.06.19.660397).


Dependencies:
doParallel
fANCOVA
ggpubr
indicspecies
maps
tidyverse
vegan
rmarkdown

Repo structure and files
data
├── ABUND_dist22_shared_list.rds
├── asvs_non_rare_long.tsv.gz
├── asvs_workable.tsv.gz
├── asv_table_nbcandem_annot_clean_long.tsv.gz
├── date2season2community.tsv
├── metagenomic_sample_name2date.tsv
├── opus2ko_non_rare_long.tsv.gz
├── opus_workable.tsv.gz
├── samo_metadata_workable.tsv
├── samo_vs_tara_workable.tsv.gz
├── tara_metadata_workable.tsv
├── kofam
│   ├── ko2desc.tsv
│   ├── KO_classification.csv
│   ├── pathway2desc.tsv
│   └── pathway2ko.tsv
└── tree_analysis
    ├── bNTI_weighted_26_1000_asvs.csv
    └── NTI_weighted_26_1000_asvs.csv

scripts
├── Figure1A_analysis.Rmd
├── Figure1B_analysis.Rmd
├── Figure2_analysis.Rmd
├── Figure3A_analysis.Rmd
├── Figure3B_analysis.Rmd
├── Figure3C_analysis.Rmd
├── Figure3D_analysis.Rmd
├── Figure4A_analysis.Rmd
├── Figure4B_analysis.Rmd
├── Figure5A_analysis.Rmd
├── Figure5B_analysis.Rmd
├── Figure6A_analysis.Rmd
├── Figure6B_analysis.Rmd
├── Figure6C_analysis.Rmd
├── Figure6D_analysis.Rmd
├── render_commands.R
└── resources
    ├── custom_bray_curtis.R
    ├── increasing_temp_simulation.R
    └── rare_indval.R
