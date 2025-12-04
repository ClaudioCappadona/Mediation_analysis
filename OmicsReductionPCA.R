OmicsReductionPCA <- function (dataframe, plotType = c("scree", "score", "loading", 
    "scoreloading"), pcType, pcNum = 5, scale = TRUE, center = TRUE, 
    eigenLoading = "loading", rotate = "none", sizeVariable = NULL, 
    sizeVariableName = NULL, colorVariable = NULL, colorVariableName = NULL, 
    shapeVariable = NULL, shapeVariableName = NULL, scaleScoreloading = TRUE, 
    firstPC = 1, secondPC = 2, loadingsName = TRUE, loadingsCutvalue = NULL, 
    loadingsCutpercent = NULL, plot = TRUE) 
{
    library(wesanderson)
    frac_var <- function(x) {
        x^2/sum(x^2)
    }
    frac_var2 <- function(x) {
        x/sum(x)
    }
    library(ggplot2)
    library(ggrepel)
    library(tidyverse)
    if (!is.null(loadingsCutpercent)) {
        if (loadingsCutpercent < 0 | loadingsCutpercent > 1) {
            stop("percentage could only be from 0 to 1")
        }
    }
    if (!is.null(loadingsCutvalue)) {
        if (loadingsCutvalue < 0) {
            stop("This is an absolute value, should be bigger than 0")
        }
    }
    if (!is.null(loadingsCutpercent) & !is.null(loadingsCutvalue)) {
        loadingsCutvalue <- NULL
    }
    if (missing(pcType)) {
        pcType <- "prcomp"
    }
    if (missing(plotType)) {
        plotType <- c("scree", "score", "loading", "scoreloading")
    }
    dataframe <- as.data.frame(dataframe)
    dataframe_numeric <- dataframe
    for (i in 1:ncol(dataframe)) {
        dataframe_numeric[, i] <- as.numeric(dataframe[, i])
    }
    if (!pcType %in% c("prcomp", "principal")) {
        stop("This method is not implemented")
    }
    if (!all(plotType %in% c("scree", "score", "loading", "scoreloading"))) {
        stop("This method is not implemented")
    }
    if (is.null(colorVariableName)) {
        colorVariableName <- "colorVariable"
    }
    if (is.null(shapeVariableName)) {
        shapeVariableName <- "shapeVariable"
    }
    if (is.null(sizeVariableName)) {
        sizeVariableName <- "sizeVariable"
    }
    pca_object <- list()
    plot_object <- list()
    if (pcType == "prcomp") {
        pca <- prcomp(dataframe_numeric, scale. = scale, center = center)
        scores <- pca$x
        loadings <- pca$rotation
        scores <- pca$x %>% as.data.frame()
        loadings <- pca$rotation %>% as.data.frame()
        if (eigenLoading == "loading") {
            for (i in 1:ncol(loadings)) {
                loadings[, i] <- loadings[, i] * pca$sdev[i]
                scores[, i] <- scores[, i]/pca$sdev[i]
            }
        }
        pca_object$scores <- scores
        pca_object$loadings <- loadings
        variance <- pca$sdev[1:pcNum] %>% as_tibble() %>% frac_var()
        variance_actual <- pca$sdev[1:pcNum] %>% as.data.frame()
    }
    if (pcType == "principal") {
        if (!rotate %in% c("none", "varimax", "quartimax", "promax", 
            "oblimin", "simplimax", "cluster")) {
            stop("This rotate method nor supported")
        }
        library(psych)
        pca <- principal(dataframe_numeric, nfactors = pcNum, 
            rotate = rotate, scores = TRUE)
        scores <- pca$scores %>% as.data.frame()
        loadings <- pca$loadings
        class(loadings) <- "matrix"
        loadings <- as.data.frame(as.matrix(loadings))
        if (eigenLoading == "eigen") {
            for (i in 1:ncol(loadings)) {
                loadings[, i] <- loadings[, i]/sqrt(pca$values)[i]
                scores[, i] <- scores[, i] * sqrt(pca$values)[i]
            }
        }
        pca_object$scores <- scores
        pca_object$loadings <- loadings
        variance <- sqrt(pca$values[1:pcNum]) %>% as_tibble() %>% 
            frac_var()
        variance_actual <- pca$values[1:pcNum] %>% as.data.frame()
    }
    if (plot == TRUE) {
        if ("scree" %in% plotType) {
            library(scales)
            library(dplyr)
            if (pcType == "prcomp") {
                scree_plot <- pca$sdev[1:pcNum] %>% as_tibble() %>% 
                  frac_var() %>% mutate(Comp = factor(colnames(pca$x)[1:pcNum], 
                  levels = colnames(pca$x)[1:pcNum])) %>% dplyr::slice(1:pcNum) %>% 
                  ggplot(aes(x = Comp, y = value)) + geom_bar(stat = "identity", 
                  fill = "#4DC5F9") + geom_hline(yintercept = 0.03, 
                  linetype = 2) + xlab("Principal Components prcomp") + 
                  scale_y_continuous(name = "Variance Explained", 
                    breaks = seq(0, 0.8, 0.1), labels = percent_format(accuracy = 5L)) + 
                  theme_classic(base_size = 14)
                plot_object$scree_plot <- scree_plot
            }
            if (pcType == "principal") {
                scree_plot <- sqrt(pca$values[1:pcNum]) %>% as_tibble() %>% 
                  frac_var() %>% mutate(Comp = factor(colnames(pca$scores), 
                  levels = colnames(pca$scores))) %>% dplyr::slice(1:pcNum) %>% 
                  ggplot(aes(x = Comp, y = value)) + geom_bar(stat = "identity", 
                  fill = "#4DC5F9") + geom_hline(yintercept = 0.03, 
                  linetype = 2) + xlab("Principal Components principal") + 
                  scale_y_continuous(name = "Variance Explained", 
                    breaks = seq(0, 0.8, 0.1), labels = percent_format(accuracy = 5L)) + 
                  theme_classic(base_size = 14)
                plot_object$scree_plot <- scree_plot
            }
        }
        if ("score" %in% plotType) {
            for (i in 1:ncol(dataframe)) {
                if (is.factor(dataframe[, i])) {
                  dataframe[, i] <- factor_samesequence(dataframe[, 
                    i])
                }
            }
            score_plot <- ggplot(scores, aes(x = scores[, firstPC], 
                y = scores[, secondPC])) + geom_point(aes(size = sizeVariable, 
                color = colorVariable, shape = shapeVariable)) + 
                geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, 
                linetype = 2) + scale_x_continuous(name = paste("Scores PC", 
                firstPC, "(", round(100 * variance[firstPC, 1], 
                  2), "%) "), limits = c(-max(abs(scores[, firstPC]), 
                abs(scores[, secondPC])), max(abs(scores[, firstPC]), 
                abs(scores[, secondPC])))) + scale_y_continuous(name = paste("Scores PC", 
                secondPC, "(", round(100 * variance[secondPC, 
                  1], 2), "%) "), limits = c(-max(abs(scores[, 
                firstPC]), abs(scores[, secondPC])), max(abs(scores[, 
                firstPC]), abs(scores[, secondPC])))) + theme_bw() + 
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
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
                length_loading <- vector()
                for (i in 1:nrow(loadings)) {
                  length_loading[i] <- sqrt((loadings[i, firstPC])^2 + 
                    (loadings[i, secondPC])^2)
                }
                loadings <- loadings[length_loading > loadingsCutvalue, 
                  ]
            }
            if (!is.null(loadingsCutpercent)) {
                length_loading <- vector()
                for (i in 1:nrow(loadings)) {
                  length_loading[i] <- sqrt((loadings[i, firstPC])^2 + 
                    (loadings[i, secondPC])^2)
                }
                loadings <- loadings[length_loading > max(abs(length_loading)) * 
                  loadingsCutpercent, ]
            }
            if (loadingsName == TRUE) {
                loading_plot <- ggplot(loadings, aes(x = loadings[, 
                  firstPC], y = loadings[, secondPC])) + geom_segment(aes(xend = loadings[, 
                  firstPC], yend = loadings[, secondPC]), x = 0, 
                  y = 0, color = "black") + geom_label_repel(aes(x = loadings[, 
                  firstPC], y = loadings[, secondPC], label = rownames(loadings)), 
                  size = 2, max.overlaps = Inf, vjust = "outward", 
                  segment.color = "grey") + geom_vline(xintercept = 0, 
                  linetype = 2) + geom_hline(yintercept = 0, 
                  linetype = 2) + scale_x_continuous(name = paste("Loadings PC", 
                  firstPC, "(", round(100 * variance[firstPC, 
                    1], 2), "%) "), limits = c(-max(abs(loadings[, 
                  firstPC]), abs(loadings[, secondPC])), max(abs(loadings[, 
                  firstPC]), abs(loadings[, secondPC])))) + scale_y_continuous(name = paste("Loadings PC", 
                  secondPC, "(", round(100 * variance[secondPC, 
                    1], 2), "%) "), limits = c(-max(abs(loadings[, 
                  firstPC]), abs(loadings[, secondPC])), max(abs(loadings[, 
                  firstPC]), abs(loadings[, secondPC])))) + theme_bw() + 
                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                plot_object$loading_plot <- loading_plot
            }
            if (loadingsName == FALSE) {
                loading_plot <- ggplot(loadings, aes(x = loadings[, 
                  firstPC], y = loadings[, secondPC])) + geom_segment(aes(xend = loadings[, 
                  firstPC], yend = loadings[, secondPC]), x = 0, 
                  y = 0, color = "black") + geom_vline(xintercept = 0, 
                  linetype = 2) + geom_hline(yintercept = 0, 
                  linetype = 2) + scale_x_continuous(name = paste("Loadings PC", 
                  firstPC, "(", round(100 * variance[firstPC, 
                    1], 2), "%) "), limits = c(-max(abs(loadings[, 
                  firstPC]), abs(loadings[, secondPC])), max(abs(loadings[, 
                  firstPC]), abs(loadings[, secondPC])))) + scale_y_continuous(name = paste("Loadings PC", 
                  secondPC, "(", round(100 * variance[secondPC, 
                    1], 2), "%) "), limits = c(-max(abs(loadings[, 
                  firstPC]), abs(loadings[, secondPC])), max(abs(loadings[, 
                  firstPC]), abs(loadings[, secondPC])))) + theme_bw() + 
                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                plot_object$loading_plot <- loading_plot
            }
        }
        if ("scoreloading" %in% plotType) {
            if (!is.null(loadingsCutvalue)) {
                length_loading <- vector()
                for (i in 1:nrow(loadings)) {
                  length_loading[i] <- sqrt((loadings[i, firstPC])^2 + 
                    (loadings[i, secondPC])^2)
                }
                loadings <- loadings[length_loading > loadingsCutvalue, 
                  ]
            }
            if (!is.null(loadingsCutpercent)) {
                length_loading <- vector()
                for (i in 1:nrow(loadings)) {
                  length_loading[i] <- sqrt((loadings[i, firstPC])^2 + 
                    (loadings[i, secondPC])^2)
                }
                loadings <- loadings[length_loading > max(abs(length_loading)) * 
                  loadingsCutpercent, ]
            }
            for (i in 1:ncol(dataframe)) {
                if (is.factor(dataframe[, i])) {
                  dataframe[, i] <- factor_samesequence(dataframe[, 
                    i])
                }
            }
            if (scaleScoreloading == TRUE) {
                if (loadingsName == TRUE) {
                  scoreloading_plot <- ggplot(scores, aes(x = scores[, 
                    firstPC], y = scores[, secondPC])) + geom_point(aes(size = sizeVariable, 
                    shape = shapeVariable, color = colorVariable)) + 
                    geom_vline(xintercept = 0, linetype = 2) + 
                    geom_hline(yintercept = 0, linetype = 2) + 
                    geom_segment(data = loadings, aes(xend = loadings[, 
                      firstPC] * (max(abs(scores[, firstPC]), 
                      abs(scores[, secondPC]))/max(abs(loadings[, 
                      firstPC]), abs(loadings[, secondPC]))), 
                      yend = loadings[, secondPC] * (max(abs(scores[, 
                        firstPC]), abs(scores[, secondPC]))/max(abs(loadings[, 
                        firstPC]), abs(loadings[, secondPC])))), 
                      x = 0, y = 0, color = "black") + geom_label_repel(data = loadings, 
                    aes(x = loadings[, firstPC] * (max(abs(scores[, 
                      firstPC]), abs(scores[, secondPC]))/max(abs(loadings[, 
                      firstPC]), abs(loadings[, secondPC]))), 
                      y = loadings[, secondPC] * (max(abs(scores[, 
                        firstPC]), abs(scores[, secondPC]))/max(abs(loadings[, 
                        firstPC]), abs(loadings[, secondPC]))), 
                      label = rownames(loadings)), size = 2, 
                    max.overlaps = Inf, vjust = "inward", hjust = "inward", 
                    color = "black", segment.color = "grey") + 
                    scale_y_continuous(name = paste("Scores PC", 
                      secondPC, "(", round(100 * variance[secondPC, 
                        1], 2), "%) "), limits = 1.01 * c(-max(abs(scores[, 
                      secondPC]), abs(scores[, firstPC])), max(abs(scores[, 
                      secondPC]), abs(scores[, firstPC]))), sec.axis = sec_axis(trans = ~./(max(abs(scores[, 
                      secondPC]), abs(scores[, firstPC]))/max(abs(loadings[, 
                      secondPC]), abs(loadings[, firstPC]))), 
                      name = paste("Loadings PC", secondPC, "(", 
                        round(100 * variance[secondPC, 1], 2), 
                        "%) "), breaks = waiver())) + scale_x_continuous(name = paste("Scores PC", 
                    firstPC, "(", round(100 * variance[firstPC, 
                      1], 2), "%) "), limits = 1.01 * c(-max(abs(scores[, 
                    firstPC]), abs(scores[, secondPC])), max(abs(scores[, 
                    firstPC]), abs(scores[, secondPC]))), sec.axis = sec_axis(trans = ~./(max(abs(scores[, 
                    secondPC]), abs(scores[, firstPC]))/max(abs(loadings[, 
                    secondPC]), abs(loadings[, firstPC]))), name = paste("Loadings PC", 
                    firstPC, "(", round(100 * variance[firstPC, 
                      1], 2), "%) "), breaks = waiver())) + theme_bw() + 
                    theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  if (is.factor(colorVariable) | is.character(colorVariable)) {
                    scoreloading_plot <- scoreloading_plot + 
                      scale_color_discrete(name = paste(colorVariableName)) + 
                      scale_size_continuous(name = paste(sizeVariableName)) + 
                      scale_shape_discrete(name = paste(shapeVariableName)) + 
                      theme_bw() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  }
                  else if (is.numeric(colorVariable) | is.integer(colorVariable)) {
                    scoreloading_plot <- scoreloading_plot + 
                      scale_color_gradientn(name = paste(colorVariableName), 
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
                if (loadingsName == FALSE) {
                  scoreloading_plot <- ggplot(scores, aes(x = scores[, 
                    firstPC], y = scores[, secondPC])) + geom_point(aes(size = sizeVariable, 
                    shape = shapeVariable, color = colorVariable)) + 
                    geom_vline(xintercept = 0, linetype = 2) + 
                    geom_hline(yintercept = 0, linetype = 2) + 
                    geom_segment(data = loadings, aes(xend = loadings[, 
                      firstPC] * (max(abs(scores[, firstPC]), 
                      abs(scores[, secondPC]))/max(abs(loadings[, 
                      firstPC]), abs(loadings[, secondPC]))), 
                      yend = loadings[, secondPC] * (max(abs(scores[, 
                        firstPC]), abs(scores[, secondPC]))/max(abs(loadings[, 
                        firstPC]), abs(loadings[, secondPC])))), 
                      x = 0, y = 0, color = "black") + scale_y_continuous(name = paste("Scores PC", 
                    secondPC, "(", round(100 * variance[secondPC, 
                      1], 2), "%) "), limits = 1.01 * c(-max(abs(scores[, 
                    secondPC]), abs(scores[, firstPC])), max(abs(scores[, 
                    secondPC]), abs(scores[, firstPC]))), sec.axis = sec_axis(trans = ~./(max(abs(scores[, 
                    secondPC]), abs(scores[, firstPC]))/max(abs(loadings[, 
                    secondPC]), abs(loadings[, firstPC]))), name = paste("Loadings PC", 
                    secondPC, "(", round(100 * variance[secondPC, 
                      1], 2), "%) "), breaks = waiver())) + scale_x_continuous(name = paste("Scores PC", 
                    firstPC, "(", round(100 * variance[firstPC, 
                      1], 2), "%) "), limits = 1.01 * c(-max(abs(scores[, 
                    firstPC]), abs(scores[, secondPC])), max(abs(scores[, 
                    firstPC]), abs(scores[, secondPC]))), sec.axis = sec_axis(trans = ~./(max(abs(scores[, 
                    secondPC]), abs(scores[, firstPC]))/max(abs(loadings[, 
                    secondPC]), abs(loadings[, firstPC]))), name = paste("Loadings PC", 
                    firstPC, "(", round(100 * variance[firstPC, 
                      1], 2), "%) "), breaks = waiver())) + theme_bw() + 
                    theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  if (is.factor(colorVariable) | is.character(colorVariable)) {
                    scoreloading_plot <- scoreloading_plot + 
                      scale_color_discrete(name = paste(colorVariableName)) + 
                      scale_size_continuous(name = paste(sizeVariableName)) + 
                      scale_shape_discrete(name = paste(shapeVariableName)) + 
                      theme_bw() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  }
                  else if (is.numeric(colorVariable) | is.integer(colorVariable)) {
                    scoreloading_plot <- scoreloading_plot + 
                      scale_color_gradientn(name = paste(colorVariableName), 
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
            if (scaleScoreloading == FALSE) {
                if (loadingsName == TRUE) {
                  scoreloading_plot <- ggplot(scores, aes(x = scores[, 
                    firstPC], y = scores[, secondPC])) + geom_point(aes(size = sizeVariable, 
                    shape = shapeVariable, color = colorVariable)) + 
                    geom_vline(xintercept = 0, linetype = 2) + 
                    geom_hline(yintercept = 0, linetype = 2) + 
                    geom_segment(data = loadings, aes(xend = loadings[, 
                      firstPC], yend = loadings[, secondPC]), 
                      x = 0, y = 0, color = "black") + geom_label(data = loadings, 
                    aes(x = loadings[, firstPC], y = loadings[, 
                      secondPC], label = rownames(loadings)), 
                    size = 2, vjust = "inward", hjust = "inward", 
                    color = "black", segment.color = "grey") + 
                    scale_y_continuous(name = paste("Scores PC", 
                      secondPC, "(", round(100 * variance[firstPC, 
                        1], 2), "%) "), limits = 1.01 * c(-max(abs(scores[, 
                      secondPC]), abs(scores[, firstPC])), max(abs(scores[, 
                      secondPC]), abs(scores[, firstPC])))) + 
                    scale_x_continuous(name = paste("Scores PC", 
                      firstPC, "(", round(100 * variance[secondPC, 
                        1], 2), "%) "), limits = 1.01 * c(-max(abs(scores[, 
                      secondPC]), abs(scores[, firstPC])), max(abs(scores[, 
                      secondPC]), abs(scores[, firstPC])))) + 
                    theme_bw() + theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank())
                  if (is.factor(colorVariable) | is.character(colorVariable)) {
                    scoreloading_plot <- scoreloading_plot + 
                      scale_color_discrete(name = paste(colorVariableName)) + 
                      scale_size_continuous(name = paste(sizeVariableName)) + 
                      scale_shape_discrete(name = paste(shapeVariableName)) + 
                      theme_bw() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  }
                  else if (is.numeric(colorVariable) | is.integer(colorVariable)) {
                    scoreloading_plot <- scoreloading_plot + 
                      scale_color_gradientn(name = paste(colorVariableName), 
                        colours = wes_palette("Zissou1")) + scale_size_continuous(name = paste(sizeVariableName)) + 
                      scale_shape_discrete(name = paste(shapeVariableName)) + 
                      theme_bw() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  }
                  plot_object$scoreloading_plot <- scoreloading_plot
                }
                if (loadingsName == FALSE) {
                  scoreloading_plot <- ggplot(scores, aes(x = scores[, 
                    firstPC], y = scores[, secondPC])) + geom_point(aes(size = sizeVariable, 
                    shape = shapeVariable, color = colorVariable)) + 
                    geom_vline(xintercept = 0, linetype = 2) + 
                    geom_hline(yintercept = 0, linetype = 2) + 
                    geom_segment(data = loadings, aes(xend = loadings[, 
                      firstPC], yend = loadings[, secondPC]), 
                      x = 0, y = 0, color = "black") + scale_y_continuous(name = paste("Scores PC", 
                    secondPC, "(", round(100 * variance[firstPC, 
                      1], 2), "%) "), limits = 1.01 * c(-max(abs(scores[, 
                    secondPC]), abs(scores[, firstPC])), max(abs(scores[, 
                    secondPC]), abs(scores[, firstPC])))) + scale_x_continuous(name = paste("Scores PC", 
                    firstPC, "(", round(100 * variance[secondPC, 
                      1], 2), "%) "), limits = 1.01 * c(-max(abs(scores[, 
                    secondPC]), abs(scores[, firstPC])), max(abs(scores[, 
                    secondPC]), abs(scores[, firstPC])))) + theme_bw() + 
                    theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  if (is.factor(colorVariable) | is.character(colorVariable)) {
                    scoreloading_plot <- scoreloading_plot + 
                      scale_color_discrete(name = paste(colorVariableName)) + 
                      scale_size_continuous(name = paste(sizeVariableName)) + 
                      scale_shape_discrete(name = paste(shapeVariableName)) + 
                      theme_bw() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  }
                  else if (is.numeric(colorVariable) | is.integer(colorVariable)) {
                    scoreloading_plot <- scoreloading_plot + 
                      scale_color_gradientn(name = paste(colorVariableName), 
                        colours = wes_palette("Zissou1")) + scale_size_continuous(name = paste(sizeVariableName)) + 
                      scale_shape_discrete(name = paste(shapeVariableName)) + 
                      theme_bw() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
                  }
                  plot_object$scoreloading_plot <- scoreloading_plot
                }
            }
        }
    }
    plot_object$object <- pca_object
    variance <- as.numeric(variance[, 1])
    variance_actual <- as.numeric(variance_actual[, 1])
    plot_object$object$variance <- variance
    plot_object$object$variance_actual <- variance_actual
    plot_object$method <- "PCA"
    return(plot_object)
}