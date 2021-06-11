sim_spatial_matrix <- function(slide_obj, umi_per_spot=NULL,
                               bleed_rate=NULL, distal_rate=NULL, 
                               cont_radius=10, gene_tokeep=NULL, 
                               kernel="gaussian",
                               seed=1){
    # Simulate spatial contaminated data based on real data
    # Args:
    #   slide_obj: a slide object created or inherited from CreateSlide().
    #   umi_per_spot (num): target average UMI per tissue spot. This does not 
    #       set the total counts for all tissue spots to be the same value.
    #   bleed_rate (num): bleeding rate
    #   distal_rate (num): distal rate
    #   cont_radius (num): contamination radius
    #   gene_tokeep (vector of chr): genes to use for simulation. 
    #       If null, use top 3000 highest expressed genes.
    #   kernel (chr): name of kernel to use. Supports "gaussian", 
    #       "linear", "laplace", "cauchy". Default: "gaussian".
    #   seed (num): random seed
    # Returns:
    #   (list) 
    #       slide (dataframe): the slide dataframe
    #       true_mat (matrix of num): true expression matrix
    #       obs_mat (matrix of int): observed contaminated matrix
    #       bleed_rate (num): bleeding rate
    #       distal_rate (num): distal rate
    #       cont_radius (num): contamination radius

    set.seed(seed)
    
    # Load slide information
    raw_data <- slide_obj@assays@data$raw  # raw data matrix
    if(is.null(raw_data)){
        stop("Cannot find raw data in input slide object.")
    }
    slide <- slide_obj@metadata$slide  # slide info
    
    slide <- slide %>% select(barcode, tissue, row, col, 
                              imagerow, imagecol, height, width)
    
    ts_bcs <- filter(slide, tissue==1)$barcode
    bg_bcs <- filter(slide, tissue==0)$barcode
    
    # reorder matrix columns to match slide
    raw_data <- raw_data[,slide$barcode]
    
    # Simulate true expression matrix
    if(is.null(gene_tokeep)){
        gene_tokeep <- rowSums(raw_data) %>% 
            sort(decreasing = T) %>% head(3000) %>% names
    }else{
        if(any(!gene_tokeep%in%rownames(raw_data))){
            stop("Some input genes are not found in data matrix.")
        }
    }
    
    true_mat <- raw_data[gene_tokeep,]
    true_mat[,bg_bcs] <- 0 
    
    raw_bg_matrix <- raw_data[gene_tokeep,bg_bcs]
    raw_ts_matrix <- raw_data[gene_tokeep,ts_bcs]
    
    
    # scale true expression to reach target sequencing depth
    if(!is.null(umi_per_spot)){
        target_total_counts <- umi_per_spot*sum(slide$tissue==1)
        ts_ratio <- target_total_counts/sum(true_mat@x)
    }else{
        ts_ratio <- sum(raw_data[gene_tokeep,])/sum(true_mat@x)
    }
    true_mat@x <- round(true_mat@x*ts_ratio)
    
    
    # Estimate bleeding rate, if not specified
    if(is.null(bleed_rate)){
        # estimated total contamination:
        # average background UMI in background spots * total spots 
        # estimated bleeding rate: 
        # total contamination / total counts
        bleed_rate <- sum(raw_bg_matrix)/ncol(raw_bg_matrix)*
            ncol(raw_data)/sum(raw_data[gene_tokeep,])
        
    }
    
    
    # Estimate distal rate, if not specified
    if(is.null(distal_rate)){
        # estimated average uniform contamination:
        # trimmed mean of counts in background spots
        # estimated distal rate:
        # total uniform contamination / total contamination
        bg_sum <- colSums(raw_bg_matrix)
        bg_quant <- quantile(bg_sum, c(0.25, 0.5))
        bg_trim <- bg_sum>=bg_quant[1] & bg_sum<=bg_quant[2]
        uniform_cont <- mean(bg_sum[bg_trim])
        
        distal_rate <- uniform_cont*ncol(raw_data)/
            (bleed_rate*sum(raw_data[gene_tokeep,]))
        
        # manually set upper bound of distal rate
        distal_rate <- min(distal_rate, 0.5) 
    }
    
    
    # estimate spot distance in pixels
    if(abs(cor(slide$imagecol, slide$col))>0.99){
        lm_tmp <- lm(slide$imagecol ~ slide$col)
    }else{
        lm_tmp <- lm(slide$imagerow ~ slide$col)
    }
    
    spot_distance <- abs(coef(lm_tmp)[2])
    
    # Calculate distance matrix and weight matrix
    slide_distance <- .calculate_euclidean_weight(
        select(slide, imagerow, imagecol)
    )
    rownames(slide_distance) <- colnames(slide_distance) <- slide$barcode
    
    local_weight_mat <- .local_kernel(slide_distance,
                                      .points_to_sdv(cont_radius, 
                                                     spot_distance),
                                      kernel)
    local_weight_mat <- local_weight_mat/rowSums(local_weight_mat)
    
    # sum the gaussian and uniform weight matrices
    sum_weight_mat <- local_weight_mat*(1-distal_rate)+
        1/ncol(raw_data)*distal_rate
    
    # Simulate contaminated expression matrix
    cont_mat <- (true_mat*bleed_rate)%*%sum_weight_mat+true_mat*(1-bleed_rate)
    
    # calculate contamination rate per spot
    true_total_counts <- colSums(true_mat[,ts_bcs])
    cont_rate <- 1-(true_total_counts*(1-bleed_rate)+
                        true_total_counts*bleed_rate*
                        diag(sum_weight_mat[ts_bcs,ts_bcs]))/
        colSums(cont_mat[,ts_bcs])
    
    # simulate observed counts from Poisson distribution
    obs_mat <- as(apply(cont_mat, 2, function(x) rpois(nrow(cont_mat),x)), 
                  "dgCMatrix")
    rownames(obs_mat) <- rownames(cont_mat)
    
    return(list(slide=slide,
                true_mat=true_mat,
                obs_mat=obs_mat,
                bleed_rate=bleed_rate,
                distal_rate=distal_rate,
                cont_rate=cont_rate,
                cont_radius=cont_radius))
}


.calculate_euclidean_weight <- function(x){
    # Calculate Euclidean distances matrix of a set of 2-d points
    # Args:
    #   x (matrix of num): Each row is the coordinates of one point.
    # Returns:
    #   (matrix of num) Euclidean distance matrix
    
    if(ncol(x)!=2){
        stop("Incorrect dimension of input coordinates.")
    }
    
    # Calculate Euclidean distances
    d_mat <- dist(x, method = "euclidean", diag = TRUE)
    d_mat <- as.matrix(d_mat)
    return(d_mat)
}


.gaussian_kernel <- function(x, sigma){
    # Gaussian kernel
    # Args:
    #   x (num): euclidean distance
    #   sigma (num): bandwidth
    # Returns:
    #   (num) Gaussian weights
    
    exp(x^2/(-2*sigma^2))
}

.linear_kernel <- function(x, sigma){
    # linear kernel
    # Args:
    #   x (num): euclidean distance
    #   sigma (num): nonzero radius
    # Returns:
    #   (num) linear weights
    abs(sigma-pmin(x,sigma))
}

.laplace_kernel <- function(x, sigma){
    # Laplac kernel
    # Args:
    #   x (num): euclidean distance
    #   sigma (num): bandwidth
    # Returns:
    #   (num) Laplace weights
    exp(-abs(x)/sigma)
}

.cauchy_kernel <- function(x, sigma){
    # Cauchy kernel
    # Args:
    #   x (num): euclidean distance
    #   sigma (num): bandwidth
    # Returns:
    #   (num) Cauchy weights
    1/(1+x^2/sigma^2)
}

.kernel_list <- list(gaussian=.gaussian_kernel,
                     linear=.linear_kernel,
                     laplace=.laplace_kernel,
                     cauchy=.cauchy_kernel)

.local_kernel <- function(x, sigma, kernel){
    if(kernel=="linear"){
        # Linear kernel has different contamination radius 
        # calculation from Gaussian,etc. 
        sigma <- 2*sigma
    }
    .kernel_list[[kernel]](x,sigma)
}

.points_to_sdv <- function(cont_radius, d){
    # Transform points radius to standard deviation (the bandwidth) in
    # Gaussian kernel. This bandwidth assures that around 95% of bled-out
    # expressions are within the circle of specified spots radius
    # Args:
    #   cont_radius (int): number of points as radius where 95% contamination
    #       goes into
    #   d (num): pixel distance between two adjacent spots in a same row
    # Returns:
    #   (num): standard deviation in the Gaussian weight function
    
    r <- cont_radius*d # radius in pixels
    return(r/2) # sigma=r/2 -> radius=2*standard deviation -> mean +- 2SD: 95%
}


run_soupx <- function(raw_mat, slide, seed=1){
    # Make SoupX runnable in spatial data
    # When default parameters not working, try other candidate cutoffs
    # If still not working, return Brah!
    
    set.seed(seed)
    
    raw_ts_mat <- raw_mat[,slide$tissue==1]
    raw_bg_mat <- raw_mat[,slide$tissue==0]
    
    # generate prelim clusters
    set.seed(SEED)
    S_tmp <- CreateSeuratObject(raw_ts_mat) %>% NormalizeData(verbose=F) %>%
        FindVariableFeatures(verbose=F) %>%
        ScaleData(verbose=F) %>%
        RunPCA(verbose=F) %>%
        FindNeighbors(verbose=F) %>%
        FindClusters(verbose=F)
    
    
    # several SoupX cutoffs to make it runnable
    tfcutoffs <- rev(c(0.01, 0.05, 0.1))
    sqcutoffs <- rev(c(0.05, 0.25, 0.5))
    
    # run soupx
    sc = SoupChannel(tod=raw_mat, toc=raw_ts_mat, calcSoupProfile=F)
    sc = estimateSoup(
        sc, soupRange = c(0, max(quantile(colSums(raw_bg_mat),0.1),100))
    )
    sc = setClusters(sc, S_tmp$seurat_clusters)
    sc_tmp = try(autoEstCont(sc,doPlot = F, forceAccept=TRUE))
    
    # try other cutoffs if default not working
    for(cutoff in seq_along(tfcutoffs)){ 
        if(class(sc_tmp)=="try-error"){
            warning("No gene. Reducing cutoff.")
            sc_tmp = try(autoEstCont(sc,doPlot = F,
                                     tfidfMin = tfcutoffs[cutoff],
                                     soupQuantile = sqcutoffs[cutoff],
                                     forceAccept=TRUE))
        }
    }
    
    sc = try(adjustCounts(sc_tmp))
    if(class(sc)=="try-error"){
        return("Brah!")
    }else{
        return(sc)
    }
    
}