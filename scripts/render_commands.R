if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x

###############################################################################
# Render Figure 1 
###############################################################################

rmarkdown::render("scripts/Figure1A_analysis.Rmd", output_file = "../results/Figure1A_analysis.pdf")
rmarkdown::render("scripts/Figure1B_analysis.Rmd", output_file = "../results/Figure1B_analysis.pdf")

###############################################################################
# Render Figure 2 
###############################################################################

rmarkdown::render("scripts/Figure2_analysis.Rmd", output_file = "../results/Figure2_analysis.pdf")

###############################################################################
# Render Figure 3
###############################################################################

rmarkdown::render("scripts/Figure3A_analysis.Rmd", output_file = "../results/Figure3A_analysis.pdf")
rmarkdown::render("scripts/Figure3B_analysis.Rmd", output_file = "../results/Figure3B_analysis.pdf")
rmarkdown::render("scripts/Figure3C_analysis.Rmd", output_file = "../results/Figure3C_analysis.pdf")
rmarkdown::render("scripts/Figure3D_analysis.Rmd", output_file = "../results/Figure3D_analysis.pdf")

###############################################################################
# Render Figure 4
###############################################################################

rmarkdown::render("scripts/Figure4A_analysis.Rmd", output_file = "../results/Figure4A_analysis.pdf")
rmarkdown::render("scripts/Figure4B_analysis.Rmd", output_file = "../results/Figure4B_analysis.pdf")

###############################################################################
# Render Figure 5
###############################################################################

rmarkdown::render("scripts/Figure5A_analysis.Rmd", output_file = "../results/Figure5A_analysis.pdf")
rmarkdown::render("scripts/Figure5B_analysis.Rmd", output_file = "../results/Figure5B_analysis.pdf")

###############################################################################
# Render Figure 6
###############################################################################

rmarkdown::render("scripts/Figure6A_analysis.Rmd", output_file = "../results/Figure6A_analysis.pdf")
rmarkdown::render("scripts/Figure6B_analysis.Rmd", output_file = "../results/Figure6B_analysis.pdf")
rmarkdown::render("scripts/Figure6C_analysis.Rmd", output_file = "../results/Figure6C_analysis.pdf")
rmarkdown::render("scripts/Figure6D_analysis.Rmd", output_file = "../results/Figure6D_analysis.pdf")
