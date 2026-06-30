
# DEG number
get_DEG_num <- function(files, path){
  map_dfr(files, function(path){
    fname <- basename(path)
    sample_size = as.integer(str_match(fname, "SampleSize_(\\d+)")[,2])
    seed <- as.integer(str_match(fname, "seed_(\\d+)")[,2])

    obj <- readRDS(path)
    nsv_keys <- names(obj)
    map_dfr(nsv_keys, function(key){
        deg_num <- sum(obj[[key]][["DEA"]]$AvsB$FDR<0.05)
        tibble(
            sample_size = factor(sample_size, levels = c(24,18,12,6)), 
            seed = factor(seed),
            n_sv = key,
            DEG_num = deg_num,
            BE = metadata(obj[[key]][["SE"]])$BE,
            LEEK = metadata(obj[[key]][["SE"]])$LEEK
        )|>
        pivot_longer(cols = c(BE, LEEK), 
             names_to = "sv_method", 
             values_to = "sv_suggested") |>
        mutate(sv_suggested = ifelse(
            as.integer(str_extract(sv_suggested, "\\d+")) > 6, # if the suggested number is larget than 6, change them to 6
            "n.sv_6",
            sv_suggested
        ))
    })
  })
}

# SV value correlation
get_SV_value_correlation <- function(files, path){
  map_dfr(files, function(path) {
    fname       <- basename(path)
    sample_size <- as.integer(str_match(fname, "SampleSize_(\\d+)")[, 2])
    seed        <- as.integer(str_match(fname, "seed_(\\d+)")[, 2])
    t <- readRDS(path)
    nsv_keys <- grep("^n\\.sv_[1-9]", names(t), value = TRUE) # don't extract n.sv_0
    sv_list <- unlist(lapply(nsv_keys, function(key) {
        cd <- as.data.frame(colData(t[[key]][["SE"]]))
        sv_cols <- grep("^SV\\d+", colnames(cd), value = TRUE)
        
        # return named list: SV1_nsv1, SV2_nsv1, SV1_nsv2, ...
        setNames(
            lapply(sv_cols, function(sv) cd[[sv]]),
            paste0(sv_cols, "_", gsub("\\.", "", key))  # e.g. SV1_nsv1, SV2_nsv2
        )
    }), recursive = FALSE) 
    
    df <- as.data.frame(sv_list)
    
    cor_matrix <- cor(df)
    as.data.frame(as.table(cor_matrix)) %>%
        mutate(sample_size = sample_size, seed = seed)
  })
}

# SV value correlation curve
get_SV_value_correlation_curve <- function(files, path){
  map_dfr(files, function(path){
    fname       <- basename(path)
    sample_size <- as.integer(str_match(fname, "SampleSize_(\\d+)")[, 2])
    seed        <- as.integer(str_match(fname, "seed_(\\d+)")[, 2])
    t <- readRDS(path)
    nsv_keys <- grep("^n\\.sv_[0-9]", names(t), value = TRUE) # n.sv_0 is only for visualization of the suggested SV numbers

    map_dfr(seq_along(nsv_keys), function(key){
        cd <- as.data.frame(colData(t[[nsv_keys[key]]][["SE"]]))
        sv_cols <- grep("^SV\\d+", colnames(cd), value = TRUE)

        maxcor <- if(key==1) NA else { # because when key =1: n.sv = 0, there is no SVs
            if (key ==2) 0 else { # when key = 2: n.sv=1, there is no correlation score
            sv_mat <- cd[ ,sv_cols]
            last_sv <- sv_mat[, ncol(sv_mat)]
            prev_sv <- sv_mat[, -ncol(sv_mat)]
            max(abs(cor(last_sv, prev_sv)))
        }}
        data.frame(
            sample_size = factor(sample_size),
            seed = factor(seed),
            nsv = factor(nsv_keys[key]),
            maxcor = maxcor,
            BE = metadata(t[[key]][["SE"]])$BE,
            LEEK = metadata(t[[key]][["SE"]])$LEEK
        )|>
        pivot_longer(cols = c(BE, LEEK), 
             names_to = "sv_method", 
             values_to = "sv_suggested") |>
        mutate(sv_suggested = ifelse(
            as.integer(str_extract(sv_suggested, "\\d+")) > 6, # if the suggested number is larget than 6, change them to 6
            "n.sv_6",
            sv_suggested
        ))
    })
  })
}

# DEG ranking correlation AUC
get_DEG_ranking_cor_auc <- function(N, files, path){
  map_dfr(files, function(path){  # map_dfr loops over every file path in files, applies the function, and row-binds all results into one dataframe.
    fname       <- basename(path)
    sample_size <- as.integer(str_match(fname, "SampleSize_(\\d+)")[, 2])
    seed        <- as.integer(str_match(fname, "seed_(\\d+)")[, 2])

    t <- readRDS(path)
    nsv_keys <- grep("^n\\.sv_", names(t), value = TRUE) # including n.sv_0

    df <- as.data.frame(setNames(
        lapply(nsv_keys, function(key) rownames(t[[key]]$DEA$AvsB)),
        gsub("n\\.sv_", "nsv", nsv_keys)  # n.sv_0 → nsv0, n.sv_1 → nsv1, etc.
    ))


    nsv_cols <- colnames(df)
    pairs_full <- expand.grid(col1 = nsv_cols, col2 = nsv_cols, stringsAsFactors = FALSE) # create the pairs of every combination of columns (nsv0~nsv0, nsv0~nsv1, nsv0~nsv2, ...)

    map_dfr(1:nrow(pairs_full), function(i){
        col1 <- pairs_full$col1[i]
        col2 <- pairs_full$col2[i]

        prop <- sapply(1:N, function(a) {
            length(intersect(df[[col1]][1:a], df[[col2]][1:a])) / a
        })
        auc <- sum((head(prop, -1) + tail(prop, -1)) / 2)
        auc_normalized <- auc / N

        tibble(
            sample_size = sample_size,
            seed = seed,
            col1 = col1, 
            col2 = col2, 
            auc = auc_normalized
        )
    })
})
}

# p skew
get_pskew <- function(files, path){
  map_dfr(files, function(path){
    fname       <- basename(path)
    sample_size <- as.integer(str_match(fname, "SampleSize_(\\d+)")[, 2])
    seed        <- as.integer(str_match(fname, "seed_(\\d+)")[, 2])

    t <- readRDS(path)
    nsv_keys <- grep("^n\\.sv_", names(t), value = TRUE) # including n.sv_0

    df <- as.data.frame(setNames(
        lapply(nsv_keys, function(key) mean(t[[key]]$DEA$AvsB$PValue - 0.5)),
        gsub("n\\.sv_", "nsv", nsv_keys)  # n.sv_0 → nsv0, n.sv_1 → nsv1, etc.
    )) %>% mutate(sample_size = sample_size, seed = seed)
})
}

# MDS separation
get_MDS_separation <- function(files, path){
  map_dfr(files, function(path){
    fname       <- basename(path)
    sample_size <- as.integer(str_match(fname, "SampleSize_(\\d+)")[, 2])
    seed        <- as.integer(str_match(fname, "seed_(\\d+)")[, 2])

    t <- readRDS(path)
    nsv_keys <- grep("^n\\.sv_", names(t), value = TRUE) # including n.sv_0

    map_dfr(nsv_keys, function(key) {
        mds <- t[[key]][["MDSdata"]]
        coords <- cbind(x = mds$var.explained[1] * mds$x,  # add the weight to the two axis
                        y = mds$var.explained[2] * mds$y)
        group <- as.factor(t[[key]][["SE"]]$group)
        
        # Calculate Total Sum of Squares (SS_total)
        # This is the sum of squared distances from each point to the global centroid
        global_centroid <- colMeans(coords)
        ss_total <- sum(rowSums(scale(coords, center = global_centroid, scale = FALSE)^2))
        
        # Calculate Between-Group Sum of Squares (SS_between)
        # This is the sum of squared distances from group centroids to the global centroid,
        # weighted by group size
        group_centroids <- aggregate(coords ~ group, FUN = mean)
        group_counts    <- table(group)
        
        ss_between <- sum(sapply(levels(group), function(g) {
            diff <- group_centroids[group_centroids$group == g, c("x", "y")] - global_centroid
            group_counts[g] * sum(diff^2)
        }))
        
        data.frame(sample_size = sample_size, 
                   seed = seed, 
                   n_sv = key, 
                   pve = ss_between / ss_total,  # !!Proportion of Variance Explained
                   BE = metadata(t[[key]][["SE"]])$BE,
                   LEEK = metadata(t[[key]][["SE"]])$LEEK
                )|>
        pivot_longer(cols = c(BE, LEEK), 
             names_to = "sv_method", 
             values_to = "sv_suggested") |>
        mutate(sv_suggested = ifelse(
            as.integer(str_extract(sv_suggested, "\\d+")) > 6, # if the suggested number is larget than 6, change them to 6
            "n.sv_6",
            sv_suggested
        ))

    })
})
}

# SV contribution to correction
get_SV_contribution_to_correction <- function(files, path){
  map_dfr(files, function(path){
    fname       <- basename(path)
    sample_size <- as.integer(str_match(fname, "SampleSize_(\\d+)")[, 2])
    seed        <- as.integer(str_match(fname, "seed_(\\d+)")[, 2])

    t <- readRDS(path)
    nsv_keys <- grep("^n\\.sv_", names(t), value = TRUE)

    map_dfr(nsv_keys, function(key){
        sv_sumAbsBs <- metadata(t[[key]][["SE"]])$sumAbsBs
        if(is.null(sv_sumAbsBs)) return(NULL)  # there is no sumAbsBs for n.sv_0
        data.frame(
            sample_size = sample_size,
            seed        = seed,
            n_sv        = key,
            SV          = names(sv_sumAbsBs),
            SV_id       = paste0(key, "_", names(sv_sumAbsBs)),
            sumAbsBs    = as.numeric(sv_sumAbsBs)
        )
    })
  })
}

# SV contribution to logFC of the DEGs
get_SV_contribution_to_logFC_DEG <- function(files, path){
  map_dfr(files, function(path){
    fname       <- basename(path)
    sample_size <- as.integer(str_match(fname, "SampleSize_(\\d+)")[, 2])
    seed        <- as.integer(str_match(fname, "seed_(\\d+)")[, 2])
    
    t <- readRDS(path)
    nsv_keys <- grep("^n\\.sv_", names(t), value = TRUE)

    map_dfr(nsv_keys, function(key){
        dea <- t[[key]][["DEA"]]
        svs <- setdiff(names(dea), c("AvsB", "allSV"))
        if(length(svs) == 0) return(NULL)
        
        genes <- rownames(dea[["AvsB"]])[which(dea[["AvsB"]]$FDR<0.05)] # only look at DEGs
        df <- data.frame(
            logFC_AvsB = dea[["AvsB"]][genes, "logFC", drop=FALSE],
            sapply(svs, function(sv) dea[[sv]][genes, "logFC", drop=FALSE])
        )
        df <- abs(df)
        rownames(df) <- genes
        colnames(df) <- c("logFC_AvsB", svs)

        SV_contribute <- max(sapply(svs, function(sv) {
            mean(df[[sv]] - df$logFC_AvsB)
        }))
    data.frame(sample_size = sample_size, 
               seed = seed, 
               n_sv = key, 
               SV_contribute = SV_contribute,
               BE = metadata(t[[key]][["SE"]])$BE,
                   LEEK = metadata(t[[key]][["SE"]])$LEEK
                )|>
        pivot_longer(cols = c(BE, LEEK), 
             names_to = "sv_method", 
             values_to = "sv_suggested") |>
        mutate(sv_suggested = ifelse(
            as.integer(str_extract(sv_suggested, "\\d+")) > 6, # if the suggested number is larget than 6, change them to 6
            "n.sv_6",
            sv_suggested
        ))
    })
})
}

# separate true and false positive:

get_SV_contribution_to_logFC_DEG_sep_isDEG <- function(files, path){
  map_dfr(files, function(path){
    fname       <- basename(path)
    sample_size <- as.integer(str_match(fname, "SampleSize_(\\d+)")[, 2])
    seed        <- as.integer(str_match(fname, "seed_(\\d+)")[, 2])
    
    t <- readRDS(path)
    nsv_keys <- grep("^n\\.sv_", names(t), value = TRUE)

    map_dfr(nsv_keys, function(key){
        dea <- t[[key]][["DEA"]]
        svs <- setdiff(names(dea), "AvsB")
        if(length(svs) == 0) return(NULL)
        
        genes <- rownames(dea[["AvsB"]])[which(dea[["AvsB"]]$FDR<0.05)] # only look at DEGs
        df <- data.frame(
            logFC_AvsB = dea[["AvsB"]][genes, "logFC", drop = FALSE],
            sapply(svs, function(sv) dea[[sv]][genes, "logFC", drop = FALSE])
        )
        df <- abs(df)
        rownames(df) <- genes
        colnames(df) <- c("logFC_AvsB", svs)
        df$isDEG <- rowData(t[[key]][["SE"]][genes],)$isDEG

        map_dfr(c(0,1), function(isdeg){
            df_sub <- df[df$isDEG == isdeg, ]
            SV_contribute <- max(sapply(svs, function(sv) {
            mean(df_sub[[sv]] - df_sub$logFC_AvsB)
            }))
                        data.frame(sample_size = sample_size, 
               seed = seed, 
               n_sv = key, 
               SV_contribute = SV_contribute,
               isDEG = isdeg,
               BE = metadata(t[[key]][["SE"]])$BE,
               LEEK = metadata(t[[key]][["SE"]])$LEEK
            )|>
        pivot_longer(cols = c(BE, LEEK), 
             names_to = "sv_method", 
             values_to = "sv_suggested") |>
        mutate(sv_suggested = ifelse(
            as.integer(str_extract(sv_suggested, "\\d+")) > 6, # if the suggested number is larget than 6, change them to 6
            "n.sv_6",
            sv_suggested
        ))
        })
    })
})
}

# false discovery rate:
get_fdr <- function(files, path){
  map_dfr(files, function(path){
    fname       <- basename(path)
    sample_size <- as.integer(str_match(fname, "SampleSize_(\\d+)")[, 2])
    seed        <- as.integer(str_match(fname, "seed_(\\d+)")[, 2])
    
    t <- readRDS(path)
    nsv_keys <- grep("^n\\.sv_", names(t), value = TRUE)

    map_dfr(nsv_keys, function(key){
        dea <- t[[key]][["DEA"]]
        
        genes <- rownames(dea[["AvsB"]])[which(dea[["AvsB"]]$FDR<0.05)] # only look at DEGs
        isDEG <- rowData(t[[key]][["SE"]][genes],)$isDEG
        fdr <- sum(isDEG == 0)/length(genes)
        data.frame(sample_size = sample_size, 
               seed = seed, 
               n_sv = key, 
               fdr = fdr,
               BE = metadata(t[[key]][["SE"]])$BE,
               LEEK = metadata(t[[key]][["SE"]])$LEEK
            )|>
        pivot_longer(cols = c(BE, LEEK), 
             names_to = "sv_method", 
             values_to = "sv_suggested") |>
        mutate(sv_suggested = ifelse(
            as.integer(str_extract(sv_suggested, "\\d+")) > 6, # if the suggested number is larget than 6, change them to 6
            "n.sv_6",
            sv_suggested
        ))
    })
  })

}

get_fdr_curve_auc <- function(N, files, path){
  map_dfr(files, function(path){
    fname       <- basename(path)
    sample_size <- as.integer(str_match(fname, "SampleSize_(\\d+)")[, 2])
    seed        <- as.integer(str_match(fname, "seed_(\\d+)")[, 2])

    t <- readRDS(path)
    nsv_keys <- grep("^n\\.sv_", names(t), value = TRUE)
    map_dfr(nsv_keys, function(key){
        dea <- t[[key]][["DEA"]]
        
        genes <- rownames(dea[["AvsB"]])[1:N] # top N in DEA results
        isDEG <- rowData(t[[key]][["SE"]][genes],)$isDEG
        prop <- sapply(1:length(genes), function(a) {
            sum(isDEG[1:a] == 0) / a
        })
        auc <- sum((head(prop, -1) + tail(prop, -1)) / 2)
        auc_normalized <- auc / N

        data.frame(
            sample_size = sample_size, 
            seed = seed, 
            n_sv = key, 
            auc = auc_normalized,
            BE = metadata(t[[key]][["SE"]])$BE,
            LEEK = metadata(t[[key]][["SE"]])$LEEK
        )|>
        pivot_longer(cols = c(BE, LEEK), 
             names_to = "sv_method", 
             values_to = "sv_suggested") |>
        mutate(sv_suggested = ifelse(
            as.integer(str_extract(sv_suggested, "\\d+")) > 6, # if the suggested number is larget than 6, change them to 6
            "n.sv_6",
            sv_suggested
        ))
    })
  })

}

get_semi_F1_score <- function(files, path){
  map_dfr(files, function(path){
    fname       <- basename(path)
    sample_size <- as.integer(str_match(fname, "SampleSize_(\\d+)")[, 2])
    seed        <- as.integer(str_match(fname, "seed_(\\d+)")[, 2])

    t <- readRDS(path)
    nsv_keys <- grep("^n\\.sv_", names(t), value = TRUE)
    map_dfr(nsv_keys, function(key){
        dea <- t[[key]][["DEA"]]
        
        genes <- rownames(dea[["AvsB"]])[which(dea[["AvsB"]]$FDR<0.05)] # only look at DEGs
        isDEG <- rowData(t[[key]][["SE"]][genes],)$isDEG
        sensitivity <- sum(isDEG == 1) / sum(rowData(t[[key]][["SE"]])$isDEG) # 500 is the total true DEGs 
        fdr <- sum(isDEG == 0)/length(genes)

        F1 <- 2/(1/(1-fdr) + 1/sensitivity)

        data.frame(
            sample_size = sample_size, 
            seed = seed, 
            n_sv = key, 
            F1 = F1,
            BE = metadata(t[[key]][["SE"]])$BE,
            LEEK = metadata(t[[key]][["SE"]])$LEEK
        )|>
        pivot_longer(cols = c(BE, LEEK), 
             names_to = "sv_method", 
             values_to = "sv_suggested") |>
        mutate(sv_suggested = ifelse(
            as.integer(str_extract(sv_suggested, "\\d+")) > 6, # if the suggested number is larget than 6, change them to 6
            "n.sv_6",
            sv_suggested
        ))
    })
  })

}

# test if p_hist is stable
get_p_hist_stability <- function(files, path){
  map_dfr(files, function(path){
    fname <- basename(path)
    sample_size = as.integer(str_match(fname, "SampleSize_(\\d+)")[,2])
    seed <- as.integer(str_match(fname, "seed_(\\d+)")[,2])
    
    obj <- readRDS(path)
    nsv_keys <- names(obj)
    map_dfr(nsv_keys, function(key){
        pvals <- obj[[key]][["DEA"]]$AvsB$PValue

        tibble(
            sample_size = sample_size, 
            seed = factor(seed),
            n_sv = key,
            gene = rownames(obj[[key]][["DEA"]]$AvsB),
            PValue = pvals
        )
    })
})
}

get_significance_ranking_cor_auc <- function(N, files, path){
  map_dfr(files, function(path){  # map_dfr loops over every file path in files, applies the function, and row-binds all results into one dataframe.
    fname       <- basename(path)
    sample_size <- as.integer(str_match(fname, "SampleSize_(\\d+)")[, 2])
    seed        <- as.integer(str_match(fname, "seed_(\\d+)")[, 2])

    t <- readRDS(path)
    nsv_keys <- grep("^n\\.sv_", names(t), value = TRUE) 
    nsv_keys <- nsv_keys[nsv_keys != "n.sv_0"] # because n.sv_0 doesnt have allSV DEA results
      
    df <- lapply(nsv_keys, function(key) {
    list(
        AvsB  = rownames(t[[key]]$DEA$AvsB),
        allSV = rownames(t[[key]]$DEA$allSV)
        )
    }) 
    names(df) <- nsv_keys
    
    df <- unlist(df, recursive = FALSE)
    
    map_dfr(nsv_keys, function(key){
        prop <- sapply(1:N, function(a) {
            length(intersect(df[[paste0(key, ".AvsB")]][1:a], df[[paste0(key, ".allSV")]][1:a])) / a
        })
        auc <- sum((head(prop, -1) + tail(prop, -1)) / 2)
        auc_normalized <- auc / N

        tibble(
            sample_size = sample_size,
            seed = seed,
            n_sv = key,
            col1 = paste0(key, ".AvsB"), 
            col2 = paste0(key, ".allSV"), 
            auc = auc_normalized,
                           BE = metadata(t[[key]][["SE"]])$BE,
               LEEK = metadata(t[[key]][["SE"]])$LEEK
            )|>
        pivot_longer(cols = c(BE, LEEK), 
             names_to = "sv_method", 
             values_to = "sv_suggested") |>
        mutate(sv_suggested = ifelse(
            as.integer(str_extract(sv_suggested, "\\d+")) > 6, # if the suggested number is larget than 6, change them to 6
            "n.sv_6",
            sv_suggested
        ))
    })
  })
}
