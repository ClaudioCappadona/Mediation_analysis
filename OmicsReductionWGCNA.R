OmicsReductionWGCNA <- function (dataframe, plotType = c("scree", "score", "loading", 
    "scoreloading"), pcType, pcNum = 5, scale = TRUE, center = TRUE, 
    eigenLoading = "loading", rotate = "none", sizeVariable = NULL, 
    sizeVariableName = NULL, colorVariable = NULL, colorVariableName = NULL, 
    shapeVariable = NULL, shapeVariableName = NULL, scaleScoreloading = TRUE, 
    firstPC = 1, secondPC = 2, loadingsName = FALSE, loadingsCutvalue = NULL, 
    loadingsCutpercent = NULL, plot = FALSE, minModuleSize, ...) 
{
    #library(WGCNA)
    Args <- list(...)
    dataframe <- as.data.frame(dataframe)
    for (i in 1:ncol(dataframe)) {
        if (is.factor(dataframe[, i])) {
            dataframe[, i] <- factor_samesequence(dataframe[, 
                i])
        }
    }
    dataframe_numeric <- dataframe
    for (i in 1:ncol(dataframe)) {
        dataframe_numeric[, i] <- as.numeric(dataframe[, i])
    }
    if (scale == TRUE) {
        dataframe_numeric <- scale(dataframe_numeric, center = rep(0, 
            ncol(dataframe_numeric)))
    }
    powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
    sft = pickSoftThreshold(dataframe_numeric, powerVector = powers, 
        verbose = 5)
    R2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
    threshold <- 0.9
    for (i in 1:10) {
        if (length(table(R2 > threshold)) == 1) {
            if (names(table(R2 > threshold)) == "FALSE") {
                threshold <- threshold - 0.1
            }
        }
    }
    softPower_position <- c()
    for (i in 1:nrow(sft$fitIndices)) {
        if (R2[i] > threshold) {
            softPower_position <- c(softPower_position, sft$fitIndices[i, 
                1])
        }
    }
    softPower <- softPower_position[which.min(softPower_position)]
    adjacency <- adjacency(dataframe_numeric, power = softPower)
    TOM <- TOMsimilarity(adjacency)
    dissTOM <- 1 - TOM
    geneTree <- hclust(as.dist(dissTOM), method = "average")
    if (missing(minModuleSize)) {
        minModuleSize <- round(ncol(dataframe)/5)
    }
    
    if ("deepSplit" %in% names(Args)) {
        message(paste("using deepSplit:", Args$deepSplit))
        dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
        deepSplit = Args$deepSplit, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
    } else {
        message("using deepSplit default value: 2")
        dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
        deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
    }
    
    dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
        deepSplit = deepSplit, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
    dynamicColors <- labels2colors(dynamicMods)
    MEList <- moduleEigengenes(dataframe, colors = dynamicColors)
    MEs <- MEList$eigengenes
    MEDiss <- 1 - cor(MEs)
    if (ncol(MEDiss) > 2) {
        METree <- hclust(as.dist(MEDiss), method = "average")
    }
    if ("MEDissThres" %in% names(Args)) {
        message(paste("using MEDissThres:", Args$MEDissThres))
        MEDissThres <- Args$MEDissThres
    } else {
        message("using MEDissThres default value: 0.25")
        MEDissThres <- 0.25
    }
    
    if (length(unique(dynamicColors)) != 1) {
        merge <- mergeCloseModules(dataframe, dynamicColors, 
            cutHeight = MEDissThres, verbose = 3)
        mergedColors <- merge$colors
        mergedMEList <- moduleEigengenes(dataframe, colors = mergedColors)
        mergedMEs <- merge$newMEs
    }
    else {
        mergedColors <- rep("grey", ncol(dataframe))
        mergedMEList <- moduleEigengenes(dataframe, colors = mergedColors)
        mergedMEs <- MEs
    }
    moduleColors <- mergedColors
    colorOrder <- c("grey", standardColors(50))
    moduleLabels <- match(moduleColors, colorOrder) - 1
    MEs <- mergedMEs
    MEList <- mergedMEList
    colored_block_pca_object <- list()
    loadings_names_list <- list()
    for (i in 1:length(table(moduleColors))) {
        colored_block <- as.data.frame(dataframe[, moduleColors == 
            names(table(moduleColors))[i], drop = FALSE])
        if (pcNum > ncol(colored_block)) {
            pc_temp <- ncol(colored_block)
        }
        else {
            pc_temp <- pcNum
        }
        temp <- OmicsReductionPCA(colored_block, plotType = plotType, 
            pcType = pcType, pcNum = pc_temp, scale = scale, 
            center = center, eigenLoading = eigenLoading, rotate = rotate, 
            sizeVariable = sizeVariable, sizeVariableName = sizeVariableName, 
            colorVariable = colorVariable, colorVariableName = colorVariableName, 
            shapeVariable = shapeVariable, shapeVariableName = shapeVariableName, 
            scaleScoreloading = scaleScoreloading, firstPC = firstPC, 
            secondPC = secondPC, loadingsName = loadingsName, 
            loadingsCutvalue = loadingsCutvalue, loadingsCutpercent = loadingsCutpercent, 
            plot = FALSE)
        colored_block_pca_object[[i]] <- temp$object
        loadings_names_list[[i]] <- rownames(colored_block_pca_object[[i]]$loadings)
    }
    names(loadings_names_list) <- names(table(moduleColors))
    names(colored_block_pca_object) <- names(table(moduleColors))
    loadings_list <- list()
    scores_list <- list()
    variance_list <- list()
    variance_actual_list <- list()
    for (i in 1:length(colored_block_pca_object)) {
        loadings_list[[i]] <- colored_block_pca_object[[i]]$loadings
        scores_list[[i]] <- colored_block_pca_object[[i]]$scores
        variance_list[[i]] <- colored_block_pca_object[[i]]$variance
        variance_actual_list[[i]] <- colored_block_pca_object[[i]]$variance_actual
        if (i == 1) {
            scores_matrix <- as.data.frame(scores_list[[i]][, 
                1])
        }
        else {
            scores_matrix <- cbind(scores_matrix, scores_list[[i]][, 
                1])
        }
    }
    module_list <- list(loadings = loadings_list, scores_matrix = scale(scores_matrix), 
        scores_list = scores_list, colored_block_pca_object = colored_block_pca_object, 
        loadings_names_list = loadings_names_list, variance_list = variance_list, 
        variance_actual_list = variance_actual_list)
    plot_object <- list()
    scores <- scale(scores_matrix) %>% as.data.frame()
    colnames(scores) <- names(table(moduleColors))
    module_list$scores_matrix <- scores
    names(loadings_list) <- names(table(moduleColors))
    loadings_list_1PC <- list()
    for (i in 1:length(loadings_list)) {
        loadings_list_1PC[[i]] <- loadings_list[[i]][, 1, drop = FALSE]
    }
    names(loadings_list_1PC) <- names(table(moduleColors))
    variance <- MEList$varExplained
    variance <- as.numeric(MEList$varExplained[1, ])
    variance_firstPC_percentage <- unlist(lapply(variance_actual_list, 
        function(vec) {
            return(vec[1])
        }))/sum(unlist(variance_actual_list))
    variance_percentage <- unlist(lapply(variance_actual_list, 
        function(vec) {
            return(sum(vec))
        }))/sum(unlist(variance_actual_list))
    names(variance) <- colnames(scores)
    module_list$variance <- variance
    module_list$variance_percentage <- variance_percentage
    module_list$variance_firstPC_percentage <- variance_firstPC_percentage
    module_list$loadings <- loadings_list_1PC
    module_list$scores <- module_list$scores_matrix
    loading_lengths <- c()
    for (i in 1:length(module_list$loadings)) {
        loading_lengths[i] <- length(module_list$loadings[[i]][, 
            1])
    }
    module_list_loadings <- module_list$loadings[order(loading_lengths)]
    module_list_variance <- module_list$variance[order(loading_lengths)]
    if (length(module_list_variance) > 1) {
        for (i in 1:(length(module_list_loadings) - 1)) {
            if (length(module_list_loadings[[i]]) == length(module_list_loadings[[i + 
                1]])) {
                if (module_list_variance[i + 1] < module_list_variance[i]) {
                  loading_lengths[i] <- loading_lengths[i] + 
                    0.5
                }
            }
        }
    }
    module_list$loadings_names_list <- module_list$loadings_names_list[order(-loading_lengths), 
        drop = FALSE]
    module_list$loadings <- module_list$loadings[order(-loading_lengths), 
        drop = FALSE]
    module_list$variance <- module_list$variance[order(-loading_lengths), 
        drop = FALSE]
    module_list$variance_firstPC_percentage <- module_list$variance_firstPC_percentage[order(-loading_lengths), 
        drop = FALSE]
    module_list$variance_percentage <- module_list$variance_percentage[order(-loading_lengths), 
        drop = FALSE]
    module_list$scores <- module_list$scores[, order(-loading_lengths), 
        drop = FALSE]
    module_list$scores_matrix <- module_list$scores_matrix[, 
        order(-loading_lengths), drop = FALSE]
    module_list$scores_list <- module_list$scores_list[order(-loading_lengths), 
        drop = FALSE]
    module_list$colored_block_pca_object <- module_list$colored_block_pca_object[order(-loading_lengths), 
        drop = FALSE]
    scores <- module_list$scores
    variance <- module_list$variance
    variance_percentage <- module_list$variance_percentage
    variance_firstPC_percentage <- module_list$variance_firstPC_percentage
    scores_matrix <- module_list$scores_matrix
    scores_list <- module_list$scores_list
    colored_block_pca_object <- module_list$colored_block_pca_object
    loadings_list_1PC <- module_list$loadings
    loadings_names_list <- module_list$loadings_names_list
    if (length(loadings_names_list) >= 2) {
        loadings_names_1 <- loadings_names_list[[firstPC]]
        loadings_names_2 <- loadings_names_list[[secondPC]]
        loadings_names_1_fixed <- loadings_names_list[[firstPC]]
        loadings_names_2_fixed <- loadings_names_list[[secondPC]]
        loadings_1 <- loadings_list_1PC[[firstPC]]
        loadings_2 <- loadings_list_1PC[[secondPC]]
        loadings_1_fixed <- loadings_list_1PC[[firstPC]]
        loadings_2_fixed <- loadings_list_1PC[[secondPC]]
    }
    else {
        loadings_names_1 <- loadings_names_list[[1]]
        loadings_names_2 <- loadings_names_list[[1]]
        loadings_names_1_fixed <- loadings_names_list[[1]]
        loadings_names_2_fixed <- loadings_names_list[[1]]
        loadings_1 <- loadings_list_1PC[[1]]
        loadings_2 <- loadings_list_1PC[[1]]
        loadings_1_fixed <- loadings_list_1PC[[1]]
        loadings_2_fixed <- loadings_list_1PC[[1]]
        firstPC <- secondPC <- 1
    }
    if (plot == TRUE) {
        if ("scree" %in% plotType) {
            library(scales)
            library(dplyr)
            module_names <- as.data.frame(factor(colnames(scores), 
                levels = colnames(scores)))
            color_names <- module_names[, 1]
            module_names <- as.data.frame(paste0("module ", 1:nrow(module_names)))
            data <- cbind.data.frame(module_names, variance)
            rownames(data) <- module_names[, 1]
            colnames(data) <- c("module_names", "variance")
            scree_plot <- ggplot(data = data) + geom_bar(mapping = aes(x = module_names, 
                y = variance_percentage), stat = "identity", 
                fill = color_names, alpha = 0.1) + geom_bar(mapping = aes(x = module_names, 
                y = variance_firstPC_percentage), stat = "identity", 
                fill = color_names) + xlab("Modules (ordered by number of variables per cluster)") + 
                scale_y_continuous(name = "Variance Explained (darker color for PC1)", 
                  breaks = seq(0, 0.8, 0.1), labels = percent_format(accuracy = 5L)) + 
                theme_classic(base_size = 14)
            plot_object$scree_plot <- scree_plot
        }
        if ("score" %in% plotType) {
            score_plot <- ggplot() + geom_point(aes(x = scores[, 
                firstPC], y = scores[, secondPC], size = sizeVariable, 
                color = colorVariable, shape = shapeVariable)) + 
                geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, 
                linetype = 2) + scale_x_continuous(name = paste0("PC1 scores from module ", 
                firstPC, " of ", length(loadings_1_fixed[, 1]), 
                " variables (", round(100 * variance[firstPC], 
                  2), "%) "), limits = c(-max(abs(scores[, firstPC]), 
                abs(scores[, secondPC])), max(abs(scores[, firstPC]), 
                abs(scores[, secondPC])))) + scale_y_continuous(name = paste0("PC1 scores from module ", 
                secondPC, " of ", length(loadings_2_fixed[, 1]), 
                " variables (", round(100 * variance[secondPC], 
                  2), "%) "), limits = c(-max(abs(scores[, firstPC]), 
                abs(scores[, secondPC])), max(abs(scores[, firstPC]), 
                abs(scores[, secondPC])))) + theme_bw() + theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank())
            if (is.factor(colorVariable) | is.character(colorVariable)) {
                score_plot <- score_plot + scale_color_discrete(name = paste(colorVariableName)) + 
                  scale_size_continuous(name = paste(sizeVariableName)) + 
                  scale_shape_discrete(name = paste(shapeVariableName)) + 
                  theme_bw() + theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())
            }
            else if (is.numeric(colorVariable) | is.integer(colorVariable)) {
                score_plot <- score_plot + scale_color_gradientn(name = paste(colorVariableName), 
                  colours = wes_palette("Zissou1")) + scale_size_continuous(name = paste(sizeVariableName)) + 
                  scale_shape_discrete(name = paste(shapeVariableName)) + 
                  theme_bw() + theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())
            }
            else if (is.null(colorVariable)) {
                if (!is.null(sizeVariable)) {
                  score_plot <- score_plot + scale_size_continuous(name = paste(sizeVariableName)) + 
                    theme_bw() + theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank())
                }
                if (!is.null(shapeVariable)) {
                  score_plot <- score_plot + scale_shape_discrete(name = paste(shapeVariableName)) + 
                    theme_bw() + theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank())
                }
            }
            plot_object$score_plot <- score_plot
        }
        if ("loading" %in% plotType) {
            if (!is.null(loadingsCutvalue)) {
                loadings_1 <- loadings_1[abs(loadings_1_fixed) > 
                  loadingsCutvalue]
                loadings_2 <- loadings_2[abs(loadings_2_fixed) > 
                  loadingsCutvalue]
                loadings_names_1 <- loadings_names_1_fixed[abs(loadings_1_fixed) > 
                  loadingsCutvalue]
                loadings_names_2 <- loadings_names_2_fixed[abs(loadings_2_fixed) > 
                  loadingsCutvalue]
            }
            if (!is.null(loadingsCutpercent)) {
                loadings_1 <- loadings_1_fixed[abs(loadings_1_fixed) > 
                  max(abs(loadings_1_fixed)) * loadingsCutpercent]
                loadings_2 <- loadings_2_fixed[abs(loadings_2_fixed) > 
                  max(abs(loadings_2_fixed)) * loadingsCutpercent]
                loadings_names_1 <- loadings_names_1_fixed[abs(loadings_1_fixed) > 
                  max(abs(loadings_1_fixed)) * loadingsCutpercent]
                loadings_names_2 <- loadings_names_2_fixed[abs(loadings_2_fixed) > 
                  max(abs(loadings_2_fixed)) * loadingsCutpercent]
            }
            loadings_1 <- unlist(loadings_1)
            loadings_2 <- unlist(loadings_2)
            loading_plot <- ggplot() + geom_segment(aes(xend = loadings_1, 
                yend = 0), x = 0, y = 0, color = "black", arrow = arrow(length = unit(0.5, 
                "cm"))) + geom_segment(aes(xend = 0, yend = loadings_2), 
                x = 0, y = 0, color = "black", arrow = arrow(length = unit(0.5, 
                  "cm")))
            if (loadingsName == TRUE) {
                loading_plot <- loading_plot + geom_label_repel(aes(x = loadings_1, 
                  y = 0, label = loadings_names_1), size = 2, 
                  max.overlaps = Inf, vjust = "outward", color = "black", 
                  segment.color = "grey") + geom_label_repel(aes(x = 0, 
                  y = loadings_2, label = loadings_names_2), 
                  size = 2, max.overlaps = Inf, vjust = "outward", 
                  color = "black", segment.color = "grey")
            }
            loading_plot <- loading_plot + geom_vline(xintercept = 0, 
                linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + 
                scale_x_continuous(name = paste0("PC1 loadings from module ", 
                  firstPC, " of ", length(loadings_1_fixed[, 
                    1]), " variables (", round(100 * variance[firstPC], 
                    2), "%) "), limits = c(-max(abs(loadings_1), 
                  abs(loadings_2)), max(abs(loadings_1), abs(loadings_2)))) + 
                scale_y_continuous(name = paste0("PC1 loadings from module ", 
                  secondPC, " of ", length(loadings_2_fixed[, 
                    1]), " variables (", round(100 * variance[secondPC], 
                    2), "%) "), limits = c(-max(abs(loadings_1), 
                  abs(loadings_2)), max(abs(loadings_1), abs(loadings_2)))) + 
                theme_bw() + theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank())
            plot_object$loading_plot <- loading_plot
        }
        if ("scoreloading" %in% plotType) {
            if (!is.null(loadingsCutvalue)) {
                loadings_1 <- loadings_1_fixed[abs(loadings_1_fixed) > 
                  loadingsCutvalue]
                loadings_2 <- loadings_2_fixed[abs(loadings_2_fixed) > 
                  loadingsCutvalue]
                loadings_names_1 <- loadings_names_1_fixed[abs(loadings_1_fixed) > 
                  loadingsCutvalue]
                loadings_names_2 <- loadings_names_2_fixed[abs(loadings_2_fixed) > 
                  loadingsCutvalue]
            }
            if (!is.null(loadingsCutpercent)) {
                loadings_1 <- loadings_1_fixed[abs(loadings_1_fixed) > 
                  max(abs(loadings_1_fixed)) * loadingsCutpercent]
                loadings_2 <- loadings_2_fixed[abs(loadings_2_fixed) > 
                  max(abs(loadings_2_fixed)) * loadingsCutpercent]
                loadings_names_1 <- loadings_names_1_fixed[abs(loadings_1_fixed) > 
                  max(abs(loadings_1_fixed)) * loadingsCutpercent]
                loadings_names_2 <- loadings_names_2_fixed[abs(loadings_2_fixed) > 
                  max(abs(loadings_2_fixed)) * loadingsCutpercent]
            }
            loadings_1 <- unlist(loadings_1)
            loadings_2 <- unlist(loadings_2)
            if (scaleScoreloading == TRUE) {
                scoreloading_plot <- ggplot() + geom_point(aes(x = scores[, 
                  firstPC], y = scores[, secondPC], size = sizeVariable, 
                  shape = shapeVariable, color = colorVariable)) + 
                  geom_vline(xintercept = 0, linetype = 2) + 
                  geom_hline(yintercept = 0, linetype = 2) + 
                  geom_segment(aes(xend = loadings_1 * (max(abs(scores[, 
                    firstPC]), abs(scores[, secondPC]))/max(abs(loadings_1), 
                    abs(loadings_2))), yend = 0), x = 0, y = 0, 
                    color = "black", arrow = arrow(length = unit(0.5, 
                      "cm"))) + geom_segment(aes(yend = loadings_2 * 
                  (max(abs(scores[, firstPC]), abs(scores[, secondPC]))/max(abs(loadings_1), 
                    abs(loadings_2))), xend = 0), x = 0, y = 0, 
                  color = "black", arrow = arrow(length = unit(0.5, 
                    "cm")))
                if (loadingsName == TRUE) {
                  scoreloading_plot <- scoreloading_plot + geom_label_repel(aes(x = loadings_1 * 
                    (max(abs(scores[, firstPC]), abs(scores[, 
                      secondPC]))/max(abs(loadings_1), abs(loadings_2))), 
                    y = 0, label = loadings_names_1), size = 2, 
                    max.overlaps = Inf, vjust = "outward", color = "black", 
                    segment.color = "grey") + geom_label_repel(aes(x = 0, 
                    y = loadings_2 * (max(abs(scores[, firstPC]), 
                      abs(scores[, secondPC]))/max(abs(loadings_1), 
                      abs(loadings_2))), label = loadings_names_2), 
                    size = 2, max.overlaps = Inf, vjust = "outward", 
                    color = "black", segment.color = "grey")
                }
                scoreloading_plot <- scoreloading_plot + scale_y_continuous(name = paste0("PC1 scores from module ", 
                  firstPC, " of ", length(loadings_1_fixed[, 
                    1]), " variables (", round(100 * variance[firstPC], 
                    2), "%) "), limits = 1.01 * c(-max(abs(scores[, 
                  secondPC]), abs(scores[, firstPC])), max(abs(scores[, 
                  secondPC]), abs(scores[, firstPC]))), sec.axis = sec_axis(trans = ~./(max(abs(scores[, 
                  secondPC]), abs(scores[, firstPC]))/max(abs(loadings_2), 
                  abs(loadings_1))), name = paste("Loadings"), 
                  breaks = waiver())) + scale_x_continuous(name = paste0("PC1 scores from module ", 
                  secondPC, " of ", length(loadings_2_fixed[, 
                    1]), " variables (", round(100 * variance[secondPC], 
                    2), "%) "), limits = 1.01 * c(-max(abs(scores[, 
                  firstPC]), abs(scores[, secondPC])), max(abs(scores[, 
                  firstPC]), abs(scores[, secondPC]))), sec.axis = sec_axis(trans = ~./(max(abs(scores[, 
                  secondPC]), abs(scores[, firstPC]))/max(abs(loadings_2), 
                  abs(loadings_1))), name = paste("Loadings"), 
                  breaks = waiver())) + theme_bw() + theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())
                if (is.factor(colorVariable) | is.character(colorVariable)) {
                  scoreloading_plot <- scoreloading_plot + scale_color_discrete(name = paste(colorVariableName)) + 
                    scale_size_continuous(name = paste(sizeVariableName)) + 
                    scale_shape_discrete(name = paste(shapeVariableName)) + 
                    theme_bw() + theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank())
                }
                else if (is.numeric(colorVariable) | is.integer(colorVariable)) {
                  scoreloading_plot <- scoreloading_plot + scale_color_gradientn(name = paste(colorVariableName), 
                    colours = wes_palette("Zissou1")) + scale_size_continuous(name = paste(sizeVariableName)) + 
                    scale_shape_discrete(name = paste(shapeVariableName)) + 
                    theme_bw() + theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank())
                }
                else if (is.null(colorVariable)) {
                  if (!is.null(sizeVariable)) {
                    scoreloading_plot <- scoreloading_plot + 
                      scale_size_continuous(name = paste(sizeVariableName)) + 
                      theme_bw() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  }
                  if (!is.null(shapeVariable)) {
                    scoreloading_plot <- scoreloading_plot + 
                      scale_shape_discrete(name = paste(shapeVariableName)) + 
                      theme_bw() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  }
                }
                plot_object$scoreloading_plot <- scoreloading_plot
            }
            if (scaleScoreloading == FALSE) {
                scoreloading_plot <- ggplot() + geom_point(aes(x = scores[, 
                  firstPC], y = scores[, secondPC], size = sizeVariable, 
                  shape = shapeVariable, color = colorVariable)) + 
                  geom_vline(xintercept = 0, linetype = 2) + 
                  geom_hline(yintercept = 0, linetype = 2) + 
                  geom_segment(aes(xend = loadings_1, yend = 0), 
                    x = 0, y = 0, color = "black", arrow = arrow(length = unit(0.5, 
                      "cm"))) + geom_segment(aes(yend = loadings_2, 
                  xend = 0), x = 0, y = 0, color = "black", arrow = arrow(length = unit(0.5, 
                  "cm"))) + scale_y_continuous(paste0("PC1 scores from module ", 
                  firstPC, " of ", length(loadings_1_fixed[, 
                    1]), " variables (", round(100 * variance[firstPC], 
                    2), "%) "), limits = 1.01 * c(-max(abs(scores[, 
                  secondPC]), abs(scores[, firstPC])), max(abs(scores[, 
                  secondPC]), abs(scores[, firstPC])))) + scale_x_continuous(paste0("PC1 scores from module ", 
                  secondPC, " of ", length(loadings_2_fixed[, 
                    1]), " variables (", round(100 * variance[secondPC], 
                    2), "%) "), limits = 1.01 * c(-max(abs(scores[, 
                  secondPC]), abs(scores[, firstPC])), max(abs(scores[, 
                  secondPC]), abs(scores[, firstPC])))) + theme_bw() + 
                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                if (is.factor(colorVariable) | is.character(colorVariable)) {
                  scoreloading_plot <- scoreloading_plot + scale_color_discrete(name = paste(colorVariableName)) + 
                    scale_size_continuous(name = paste(sizeVariableName)) + 
                    scale_shape_discrete(name = paste(shapeVariableName)) + 
                    theme_bw() + theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank())
                }
                else if (is.numeric(colorVariable) | is.integer(colorVariable)) {
                  scoreloading_plot <- scoreloading_plot + scale_color_gradientn(name = paste(colorVariableName), 
                    colours = wes_palette("Zissou1")) + scale_size_continuous(name = paste(sizeVariableName)) + 
                    scale_shape_discrete(name = paste(shapeVariableName)) + 
                    theme_bw() + theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank())
                }
                else if (is.null(colorVariable)) {
                  if (!is.null(sizeVariable)) {
                    scoreloading_plot <- scoreloading_plot + 
                      scale_size_continuous(name = paste(sizeVariableName)) + 
                      theme_bw() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  }
                  if (!is.null(shapeVariable)) {
                    scoreloading_plot <- scoreloading_plot + 
                      scale_shape_discrete(name = paste(shapeVariableName)) + 
                      theme_bw() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  }
                }
                plot_object$scoreloading_plot <- scoreloading_plot
            }
        }
    }
    plot_object$object <- module_list
    new_module_names <- paste0("module", 1:length(plot_object$object$variance))
    colnames(plot_object$object$scores_matrix) <- new_module_names
    colnames(plot_object$object$scores) <- new_module_names
    names(plot_object$object$loadings) <- new_module_names
    names(plot_object$object$loadings_names_list) <- new_module_names
    names(plot_object$object$variance) <- new_module_names
    plot_object$method <- "WGCNA"
    return(plot_object)
}
