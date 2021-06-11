# The analysis pipeline supports multiple samples as input
# However, it will write outputs in separate directories based on sample names

# Define general paths
RData_dir <- c("~/Google Drive/Hallu/codes/ckgroup/spatial_data/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/V1_Mouse_Kidney/V1_Mouse_Kidney.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/V1_Breast_Cancer_Block_A_Section_2/V1_Breast_Cancer_Block_A_Section_2.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/V1_Human_Lymph_Node/V1_Human_Lymph_Node.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/Targeted_Visium_Human_SpinalCord_Neuroscience/Targeted_Visium_Human_SpinalCord_Neuroscience.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/spatialLIBD/sp_151507/sp_151507.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/spatialLIBD/sp_151508/sp_151508.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/spatialLIBD/sp_151669/sp_151669.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/spatialLIBD/sp_151670/sp_151670.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/spatialLIBD/sp_151673/sp_151673.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/spatialLIBD/sp_151674/sp_151674.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/1_S1_manual/1_S1_manual.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/2_S2_manual/2_S2_manual.RData",
               "~/Google Drive/Hallu/codes/ckgroup/spatial_data/3_S3_manual/3_S3_manual.RData")

tif_dir <- c("~/Google Drive/Hallu/codes/ckgroup/spatial_data/manuscript/figures/Rplots/contamination")

###############################
# Generate visualization report
# Specify analysis outputs from `spatial_analysis_pipeline.R`

sample_names <- basename(dirname(RData_dir))
html_title <- paste0("Spatial data visualization: ",sample_names)
out_html <- gsub(".RData\\>","_fig.html",RData_dir)

for(i in seq_along(sample_names)){
    
    rmarkdown::render("~/Google Drive/Hallu/codes/ckgroup/spatial_data/spatial_visualization_pipeline_for_manuscript.Rmd",
                      params=list(html_title=html_title[i], RData_dir=RData_dir[i],
                                  tif_dir=tif_dir,rm_border=TRUE),
                      output_dir = dirname(RData_dir[i]),
                      output_file = out_html[i])
    
}
