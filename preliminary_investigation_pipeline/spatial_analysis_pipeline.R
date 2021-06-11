# The analysis pipeline supports multiple samples as input
# However, it will write outputs in separate directories based on sample names

# Define general paths
data_dir <- c("~/Google Drive/Hallu/DATA/10x/Spatial/V1_Adult_Mouse_Brain",
              "~/Google Drive/Hallu/DATA/10x/Spatial/Targeted_Visium_Human_SpinalCord_Neuroscience",
              "~/Google Drive/Hallu/DATA/10x/Spatial/V1_Mouse_Kidney",
              "~/Google Drive/Hallu/DATA/10x/Spatial/V1_Human_Lymph_Node",
              "~/Google Drive/Hallu/DATA/10x/Spatial/V1_Breast_Cancer_Block_A_Section_2",
              "~/Google Drive/Hallu/DATA/10x/Spatial/mix_species/1_S1_manual",
              "~/Google Drive/Hallu/DATA/10x/Spatial/mix_species/2_S2_manual",
              "~/Google Drive/Hallu/DATA/10x/Spatial/mix_species/3_S3_manual",
              "~/Google Drive/Hallu/DATA/10x/Spatial/mix_species/5_S4_manual")

#########
# The following sample names and paths are generated automatically.
# If they generate invalid paths, manually edit them.
#########

# Define samples
sample_names <- basename(data_dir)

image_dir <- numeric(length(data_dir))
for(i in seq_along(data_dir)){
    candidate_spatial_dir <- list.files(data_dir[i],
                                        pattern = "spatial",
                                        full.names = TRUE)
    candidate_spatial_dir <- candidate_spatial_dir[dir.exists(candidate_spatial_dir)]

    if("spatial"%in%list.files(candidate_spatial_dir)){
        candidate_spatial_dir <- file.path(candidate_spatial_dir, "spatial")
    }
    image_dir[i] <- candidate_spatial_dir
}


# Define image-related paths
image_paths <- file.path(image_dir,"tissue_lowres_image.png")
scalefactor_paths <- file.path(image_dir,"scalefactors_json.json")
tissue_paths <- file.path(image_dir,"tissue_positions_list.csv")

# Define gene expression matrix paths
raw_matrix_paths <- numeric(length(data_dir))
for(i in seq_along(data_dir)){
    candidate_raw_dir <- list.files(data_dir[i],
                                        pattern = "raw_feature_bc_matrix",
                                        full.names = TRUE)
    candidate_raw_dir <- candidate_raw_dir[dir.exists(candidate_raw_dir)]

    if("raw_feature_bc_matrix"%in%list.files(candidate_raw_dir)){
        candidate_raw_dir <- file.path(candidate_raw_dir,
                                           "raw_feature_bc_matrix")
    }
    raw_matrix_paths[i] <- candidate_raw_dir
}

# Define output directory
out_dir <- paste0("~/Google Drive/Hallu/codes/ckgroup/spatial_data/",sample_names)
RData_dir <- paste0(out_dir,"/",sample_names,".RData")


#######################
# Check valid input files
all_dirs <- c(data_dir, image_dir, image_paths,
              scalefactor_paths, tissue_paths, raw_matrix_paths)

if(any(!file.exists(all_dirs))){
    bad_dirs <- paste(all_dirs[!file.exists(all_dirs)],collapse = "\n")
    stop("Invalid input directories or files:\n", bad_dirs)
}


#######################
library(ggplot2)
library(Matrix)
library(rjson)
library(grid)
library(readbitmap)
library(Seurat)
library(scran)
library(scuttle)
library(dplyr)
library(hdf5r)
library(data.table)
SEED <- 2021

#####################

# Load downsampled image
images_cl <- list()

for (i in 1:length(sample_names)) {
    images_cl[[i]] <- read.bitmap(image_paths[i])
}

height <- list()

for (i in 1:length(sample_names)) {
    height[[i]] <-  data.frame(height = nrow(images_cl[[i]]))
}

height <- bind_rows(height)

width <- list()

for (i in 1:length(sample_names)) {
    width[[i]] <- data.frame(width = ncol(images_cl[[i]]))
}

width <- bind_rows(width)

grobs <- list()
for (i in 1:length(sample_names)) {
    grobs[[i]] <- rasterGrob(images_cl[[i]], width=unit(1,"npc"), height=unit(1,"npc"))
}

images_tibble <- tibble(sample=sample_names, grob=grobs)
images_tibble$height <- height$height
images_tibble$width <- width$width

scales <- list()

for (i in 1:length(sample_names)) {
    scales[[i]] <- fromJSON(file = scalefactor_paths[i])
}


# Load gene expression matrix
raw_matrix <- list()

for (i in 1:length(sample_names)) {
    raw_matrix[[i]] <- Read10X(raw_matrix_paths[i])
}


# Make summary data.frames
bcs <- list()

for (i in 1:length(sample_names)) {
    bcs[[i]] <- read.csv(tissue_paths[i],
                         col.names=c("barcode","tissue","row",
                                     "col","imagerow","imagecol"),
                         header = FALSE)
    # scale tissue coordinates for lowres image
    bcs[[i]]$imagerow <- bcs[[i]]$imagerow * scales[[i]]$tissue_lowres_scalef
    bcs[[i]]$imagecol <- bcs[[i]]$imagecol * scales[[i]]$tissue_lowres_scalef
    bcs[[i]]$tissue <- as.factor(bcs[[i]]$tissue)
    bcs[[i]]$height <- height$height[i]
    bcs[[i]]$width <- width$width[i]
}

names(bcs) <- sample_names

# total UMI and proportion of MT genes per spot
umi_sum <- list()
prop_mt <- list()
human_sum <- mouse_sum <- human_prop <- list()

for (i in seq_along(sample_names)) {
    mt_gene <- grepl("MT-|mt-",rownames(raw_matrix[[i]]))
    total_counts <- Matrix::colSums(raw_matrix[[i]])

    if(any(mt_gene)){
        mt_counts <- Matrix::colSums(raw_matrix[[i]][mt_gene,])
    }else{
        warning("No MT gene detected in ",sample_names[i])
        mt_counts <- NA
    }

    prop_mt[[i]] <- data.frame(barcode =  colnames(raw_matrix[[i]]),
                               mt_prop= mt_counts/total_counts)

    umi_sum[[i]] <- data.frame(barcode =  colnames(raw_matrix[[i]]),
                               sum_umi = total_counts)

    # for mix-species data

    if(any(grepl("\\<GRCh38_",rownames(raw_matrix[[i]])))){
        hm_gene <- grepl("\\<GRCh38_",rownames(raw_matrix[[i]]))

        human_sum[[i]] <- data.frame(barcode =  colnames(raw_matrix[[i]]),
                                     sum_human_umi = colSums(raw_matrix[[i]][hm_gene,]))
        mouse_sum[[i]] <- data.frame(barcode =  colnames(raw_matrix[[i]]),
                                     sum_mouse_umi = colSums(raw_matrix[[i]][!hm_gene,]))
        human_prop[[i]] <- data.frame(barcode =  colnames(raw_matrix[[i]]),
                                      prop_human = colSums(raw_matrix[[i]][hm_gene,])/colSums(raw_matrix[[i]]))

    }else{
        human_sum[[i]] <- mouse_sum[[i]] <- human_prop[[i]] <- NA
    }

}
names(umi_sum) <- names(prop_mt) <- names(human_sum) <-
    names(mouse_sum) <- names(human_prop) <- sample_names


# total nonzero gene per spot
gene_sum <- list()

for (i in 1:length(sample_names)) {
    gene_sum[[i]] <- data.frame(barcode =  colnames(raw_matrix[[i]]),
                                sum_gene = Matrix::colSums(raw_matrix[[i]] != 0))

}
names(gene_sum) <- sample_names


# Merge summary data together
for (i in 1:length(sample_names)) {
    bcs[[i]] <- bcs[[i]] %>%
        merge(umi_sum[[i]], by = c("barcode")) %>%
        merge(gene_sum[[i]], by = c("barcode")) %>%
        merge(prop_mt[[i]], by = c("barcode"))

    if(is.data.frame(human_sum[[i]])){
        bcs[[i]] <- bcs[[i]] %>%
            merge(human_sum[[i]], by = c("barcode")) %>%
            merge(mouse_sum[[i]], by = c("barcode")) %>%
            merge(human_prop[[i]], by = c("barcode"))
    }
}


# labelling below/above cutoff
cutoffs <- c(0,10,50,100,300,500, 1000, 2000, 5000, 10000, 15000, 20000)
for(i in 1:length(sample_names)){

    for(cutoff in cutoffs){
        bcs[[i]] <- cbind(bcs[[i]],as.integer(bcs[[i]]$sum_umi>cutoff))
    }

    colnames(bcs[[i]])[
        ncol(bcs[[i]])+seq_along(cutoffs)-length(cutoffs)
        ] <- paste0("cutoff_",cutoffs)
}

# scran normalization of raw matrix

norm_matrix <- list()

for(i in seq_along(sample_names)){
    set.seed(SEED)
    sce <- SingleCellExperiment(list(counts=raw_matrix[[i]]))
    clusters <- quickCluster(sce, min.size=200)
    sce <- computeSumFactors(sce, cluster=clusters)
    min_sf <- min(sce$sizeFactor[sce$sizeFactor>0])
    # deal with potential nonpositive size factors
    sce$sizeFactor[sce$sizeFactor<=0] <- min_sf
    sce <- logNormCounts(sce)

    norm_matrix[[i]] <- logcounts(sce)
}

# UMAP and clustering using Seurat

S_obj <- list()

for(i in seq_along(sample_names)){
    set.seed(SEED)
    S_obj[[i]] <- CreateSeuratObject(counts = norm_matrix[[i]], project = "spatial") %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>%
        ScaleData %>%
        RunPCA %>%
        RunUMAP(dims=1:25) %>%
        FindNeighbors(dims = 1:25) %>%
        FindClusters

    S_obj[[i]]$bc_label <- bcs[[i]]$tissue
    bcs[[i]]$clust <- S_obj[[i]]$seurat_clusters
}


# Summarize by sample names and write output as RData

for(i in seq_along(sample_names)){

    analysis_out <- list()
    analysis_out$sample_name <- sample_names[i]
    analysis_out$bcs <- bcs[[i]]
    analysis_out$raw_matrix <- raw_matrix[[i]]
    analysis_out$norm_matrix <- norm_matrix[[i]]
    analysis_out$S_obj <- S_obj[[i]]
    analysis_out$images_tibble <- images_tibble[i,]

    if(!dir.exists(out_dir[i])){
        dir.create(out_dir[i],recursive = T)
    }
    save(analysis_out, file = RData_dir[i])
}



###############################
# Generate visualization report
# Specify analysis outputs from `spatial_analysis_pipeline.R`

html_title <- paste0("Spatial data visualization: ",sample_names)
out_html <- gsub("RData\\>","html",RData_dir)

for(i in seq_along(sample_names)){
    rmarkdown::render("~/Google Drive/Hallu/codes/ckgroup/spatial_data/spatial_visualization_pipeline.Rmd",
                      params=list(html_title=html_title[i], RData_dir=RData_dir[i]),
                      output_dir = out_dir[i],
                      output_file = out_html[i])
}
