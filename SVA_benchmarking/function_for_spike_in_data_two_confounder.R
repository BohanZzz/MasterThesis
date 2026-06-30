
# Main testing function
Lets_test <- function(se, dge, n_SV, MDS_name){
# se: the original one without SVA correction (but after filtering). It should have a colData called "group", with the group that we want to test
# dge: from the original se
# n_SV: number of SV
# MDS_name: should be folder/xxsamples_seed_xx

    test_res_list <- list()
    group <- se$group

    if(n_SV==0){
        se_sv <- se
        test_res_list[["SE"]] <- se_sv

        # MDS logcpm
        col <- se_sv$group
        col <- c(A = "firebrick3", B = "steelblue2")[col]
        pdf(paste0(MDS_name, "_logcpm.pdf"), width = 8, height = 6)
        limma::plotMDS(assay(se_sv, "logcpm"), labels = paste0(se_sv$group, "_", se_sv$confounders), col = col, main = paste("MDS: ", basename(MDS_name)," - n.sv =", n_SV, ", corrected value"))
        dev.off()
        test_res_list[["MDSdata"]] <- plotMDS(assay(se_sv, "logcpm"), labels = paste0(se_sv$group, "_", se_sv$confounders), col = col, plot = FALSE)

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
        invisible(capture.output({
            set.seed(96)
            se_sv <- svacor_PL(se, ~ group, n.sv = n_SV) }))

        # MDS plot for each n.sv
        test_res_list[["SE"]] <- se_sv
        col <- se_sv$group
        col <- c(A = "firebrick3", B = "steelblue2")[col]
        pdf(paste0(MDS_name, "_n.sv_", n_SV, ".pdf"), width = 8, height = 6)
        limma::plotMDS(assay(se_sv, "corrected"), labels = paste0(se_sv$group, "_", se_sv$confounders), col = col, main = paste("MDS: ", basename(MDS_name)," - n.sv =", n_SV, ", corrected value"))
        dev.off()
        test_res_list[["MDSdata"]] <- plotMDS(assay(se_sv, "corrected"), labels = paste0(se_sv$group, "_", se_sv$confounders), col = col, plot = FALSE)

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

# subsampling function
# to do the subsampling and group assignment
subsampling <- function(se, sampleSize, seed){
    # se: the original dataset
    # sampleSize: the size after sub sampling
    # seed: because it's a random sampling and we want to see how stable the effect is, so we use multiple seeds to do multiple round of sampling
    target <- c(
        Male_InsulRes = sampleSize/6, 
        Male_InsulSens = sampleSize/6,
        Female_InsulRes = sampleSize/6, 
        Female_InsulSens = sampleSize/2) 
    
    set.seed(seed)
    idx_sub <- unlist(lapply(split(seq_len(ncol(se)), se$confounders), \(x) {
    conf_name <- as.character(unique(se$confounders[x]))
    sample(x, target[conf_name])
    }))
    se_sub <- se[, idx_sub]

    # randomly assign them into two groups (but make sure within two confounders, the groups are balanced)
    props <- list(
        Female_InsulSens = c(A = 1, B = 2),
        Female_InsulRes = c(A = 1),
        Male_InsulSens = c(A = 1),
        Male_InsulRes = c(B = 1)
    )
    idx_by_confounder <- split(seq_len(ncol(se_sub)), se_sub$confounders)
    assignment <- dplyr::bind_rows(lapply(names(idx_by_confounder), function(conf){
        x <- idx_by_confounder[[conf]]
        n <- length(x)
        # Scale proportions to actual sample size
        p <- props[[conf]] # the proportion between two groups in this batch
        total_parts <- sum(p)  # overall weight 
        counts <- round(p / total_parts * n) # how many should be A and should be B in the samplesize we choose

        set.seed(seed)
        data.frame(
            i = x, 
            group = sample(rep(names(counts), times = counts))
        )
    }))
    se_sub$group[assignment$i] <- assignment$group
    return(se_sub)
}

filterLowCounts <- function(se, min_count=20){
    dge <- DGEList(assay(se, "counts"))
    dge <- calcNormFactors(dge)
    design <- model.matrix(~ 0 + group, data = as.data.frame(colData(se)))
    keep <- filterByExpr(dge, min.count = min_count, design = design)
    se <- se[keep, ]
    return(se)
}
