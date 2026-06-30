
# subsampling, group assignment, and adding artificial batch effects
subsampling <- function(se, sampleSize, seed){
    set.seed(seed)
    se_sub <- se[,sample(1:ncol(se), sampleSize, replace = FALSE)]

    # we also need the artificial coefficient matrix
    n <- nrow(se_sub)
    group <- factor(rep(c("A","B"),each = sampleSize/2))
    batch1 <- rep(c(0, 1, 0),c((sampleSize/2-sampleSize/6), sampleSize/2, sampleSize/6)) # batch 1 partially correlated to group (2:1)
    batch2 <- c(rnorm(sampleSize/2, 0),rnorm(sampleSize/2, 1))[rep(1:(sampleSize/2),each = 2) + 0:1 * (sampleSize/2)] # linear batch effect uncorrelated to group
# take total samplesize = 24 as an example:
# rnorm(12, 0): 12 samples, mean = 0
# rnorm(12, 1): 12 samples, mean = 1
# rep(1:12,each = 2) + 0 : 1 * 12 ==> interleaves one sample from group A and one from group B. To make it uncorrelated to the group

# start to simulate the effect size of each coefficient
    isDEG <- gco <- bco2 <- bco1 <- rep(0, n)  # initialize them to 0
    bco1[sample.int(n, floor(n)/3)] <- rnorm(floor(n)/3, sd=1.5)  # batch 1 affects 1/3 genes
    bco2[sample.int(n, floor(n)/3)] <- rnorm(floor(n)/3, sd=1.3)  # batch 2 affects 1/3 genes
    # simulate DEGs
    lfcs <- c(0.5, 0.75, 1, 1, 1.5, 2)
    deg_index <- sample.int(n, 500)
    gco[deg_index] <- sample(c(lfcs, -lfcs), 500, replace=TRUE) # randomly choose 500 genes and give them logFC. Other genes' gco = 0
    isDEG[deg_index] <- TRUE

    se_sub$batch1 <- batch1
    se_sub$batch2 <- batch2
    se_sub$group <- group
    rowData(se_sub)$batch1 <- bco1
    rowData(se_sub)$batch2 <- bco2
    rowData(se_sub)$simulated_lfc <- gco
    rowData(se_sub)$isDEG <- isDEG
    return(se_sub)
}

# simulate the count matrix:
simulateData <- function(se, seed, varfactor = 1){
# From se we need:
    # a: expression matrix (raw counts matrix, genes x samples)
    # beta: effect size matrix (coefficient matrix, genes x model terms)
    # original_MM: the design matrix from orignal raw counts, only used to calculate the dispersion. In our case: ~1 (since the original batch effect are all kicked out)
    # predict_MM: the design matrix for the predicted data. Using this to make the predicted data have these characteristics: In our case: ~ batch1+batch2+simulated_lfc
# varfactor: dispersion inflation factor, controls how much you inflate the estimated dispersions beyond what was observed in your real data. (1=no change. above 1 = more inflatoin)

# In summary: the predicted counts will have: (1) same rowMeans and dispersion as the original raw counts; (2) same batch effects as we set beforehand. 
    set.seed(seed)
# =============== 1. Align genes between a and beta ==================
    a <- assay(se, "counts")
    beta <- cbind(batch1=rowData(se)$batch1, batch2=rowData(se)$batch2, groupB=rowData(se)$simulated_lfc)
    rownames(beta) <- rownames(se)

    i <- intersect(row.names(a), row.names(beta))
    a <- a[i,]
    beta <- as.data.frame(beta[i,,drop=FALSE])

    lib <- mean(colSums(a))  # original lib size, used for step5
# =============== 2. Estimate dispersion =============================
    original_MM <- model.matrix(~ 1, data = colData(se))
    disp <- estimateDisp(calcNormFactors(DGEList(a)), original_MM)
    disp <- disp$tagwise.dispersion * as.numeric(varfactor)
    nb_size <- 1 / disp  # convert to negative binomial size parameter

# =============== 3. Prepare expression matrix =======================
    a <- log1p(a) # log transform the counts, stabilizing variance
    beta[["(Intercept)"]] <- rowMeans(a)  
    # the intercept column in beta is set to rowMeans(a) (on the log scale), 
    # so that the baseline expression level for each gene comes from the observed data

# =============== 4. Compute mu (predicted mean expression) ==========
    predict_MM <- model.matrix(~ batch1 + batch2 + group, data = colData(se))
    mu <- t(predict_MM %*% t(beta[,colnames(predict_MM),drop=FALSE])) # mu: the predicted mean expression matrix
    # predict_MM %*% t(beta): matrix multiply to get predicted values
    # %*% : matrix multiplication on two matrices

    mu <- as.matrix(exp(mu)-1) # back-transform from log scale
    mu[mu < 0 | is.na(mu)] <- 0   # clean negatives/NAs

# =============== 5. Simulate nsim replicates =======================
    libsizes <- rnorm(nrow(predict_MM), lib, sd = 0.1 * lib) # library size is randomly drawn with normal distribution, based on the mean lib size from the raw counts matrix.
        # mean: original libsize, sd: 10% of libsize (small variation between samples)

        # scale each sample's mu column to sum to that sample's library size
        # this makes expected counts proportional to library depth
        # In summary: this ensures proportions of counts stay the same across genes, but the total depth matches the library size
    mat <- mu
    for(j in 1:ncol(mat)) mat[ ,j] <- libsizes[j] * mat[ ,j] / sum(mat[ ,j])
        
    nb_size_rep <- rep(nb_size, times = ncol(mat))
        # nb_size: one value per gene. But rnbinom() needs one size value per count. So this just makes the gene dispersions to repeat once per sample

    counts <- rnbinom(n = length(mat), size = nb_size_rep, mu = as.numeric(mat))
        # here for each count, on random NB count is drawn using the value as mu and the corresponding nb_size_rep as size. 
        # higher size = less overdispersoin = counts closer to mu

    e <- matrix(counts, nrow = nrow(mat), ncol = ncol(mat)) # Reshape back into a matrix

# =============== 6. Attach row/column names ========================
    rownames(e) <- rownames(beta)
    colnames(e) <- rownames(predict_MM)
     
    se_pred <- SummarizedExperiment(list(counts = e), colData = colData(se), rowData = rowData(se))
    se_pred
}

# counts filtering
filterLowCounts <- function(se, min_count=20){
    dge <- DGEList(assay(se, "counts"))
    dge <- calcNormFactors(dge)
    design <- model.matrix(~ 0 + group, data = as.data.frame(colData(se)))
    keep <- filterByExpr(dge, min.count = min_count, design = design)
    se <- se[keep, ]
    return(se)
}

# main test:
Lets_test <- function(se, dge, n_SV, MDS_name){
# se: the se without SVA correction (but after filtering)
# dge: from the original se
# n_SV: number of SV
# MDS_name: should be folder/xxsamples_seed_xx

    test_res_list <- list()
    group <- se$group

    assay(se) <- as.matrix(assay(se))

    if(n_SV==0){
        se_sv <- se
        test_res_list[["SE"]] <- se_sv
        # MDS logcpm
        groupcol <- se$group
        groupcol <- c(A = "orange2", B = "navy")[groupcol]
        pdf(paste0(MDS_name, "_logcpm.pdf"), width = 8, height = 6)
        limma::plotMDS(assay(se_sv, "logcpm"), labels = paste0(se_sv$group, "_", se_sv$batch1, "_", se_sv$batch2), col = groupcol)
        dev.off()
        test_res_list[["MDSdata"]] <- plotMDS(assay(se_sv, "logcpm"), labels = paste0(se_sv$group, "_", se_sv$batch1, "_", se_sv$batch2), col = groupcol, plot = FALSE)

        design_0 <- model.matrix(~ 0 + group, data = as.data.frame(colData(se)))
        contrast_0 <- limma::makeContrasts(
            AvsB = groupA - groupB,
            levels = design_0
        )
        y_0 <- estimateDisp(dge, design_0)
        fit_0 <- glmFit(y_0, design_0)
        lrt_0 <- glmLRT(fit_0, contrast = contrast_0)
        test_res_list[["DEA"]] <- c(list(AvsB = as.data.frame(topTags(lrt_0, n = Inf)))) # because when there are SVs, we save DEA results as a list. So here we also save it as a list
    }
    else{
        # SVA correction
        set.seed(96)
        se_sv <- svacor_PL(se, ~ group, n.sv = n_SV)
    
        test_res_list[["SE"]] <- se_sv
        # MDS plot for each n.sv
        groupcol <- se$group
        groupcol <- c(A = "orange2", B = "navy")[groupcol]
        pdf(paste0(MDS_name, "_n.sv_", n_SV, ".pdf"), width = 8, height = 6)
        limma::plotMDS(assay(se_sv, "corrected"), labels = paste0(se_sv$group, "_", se_sv$batch1, "_", se_sv$batch2), col = groupcol, main = paste("MDS: ", basename(MDS_name)," - n.sv =", n_SV, ", corrected value"))
        dev.off()
        test_res_list[["MDSdata"]] <- limma::plotMDS(assay(se_sv, "corrected"), labels = paste0(se_sv$group, "_", se_sv$batch1, "_", se_sv$batch2), col = groupcol, plot = FALSE)

        # DEA
        sv_cols <- paste0("SV", 1:n_SV)
        formula_sv <- as.formula(paste("~0 + group +", paste(sv_cols, collapse = "+")))
        design_sv <- model.matrix(formula_sv, data = as.data.frame(colData(se_sv)))

        # add all the SVs to the contrast so we can extract their contribution to logFC
        contrast_list <- c(list(AvsB = "groupA - groupB"), setNames(as.list(sv_cols), sv_cols)) 
        contrast_sv <- do.call(limma::makeContrasts, c(contrast_list, list(levels = design_sv)))

        y_sv <- estimateDisp(dge, design_sv)
        fit_sv <- glmFit(y_sv, design_sv)
        lrt_sv_AvsB <- glmLRT(fit_sv, contrast = contrast_sv[,"AvsB"])

        # test for all the SVs
        lrt_all_sv <- glmLRT(fit_sv, paste0("SV", 1:n_SV))

        # save DEA results as a list
        test_res_list[["DEA"]] <- c(
            list(AvsB = as.data.frame(topTags(lrt_sv_AvsB, n = Inf))),
            setNames(
                lapply(sv_cols, function(sv){
                    lrt_sv <- glmLRT(fit_sv, contrast = contrast_sv[,sv])
                    as.data.frame(topTags(lrt_sv, n = Inf))
                }),
                sv_cols
            ),
            list(allSV = as.data.frame(topTags(lrt_all_sv, n = Inf)))
        )
    }

    return(test_res_list)
}
