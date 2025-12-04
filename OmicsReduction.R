OmicsReduction <- function (dataframe, plottype = c("scree", "score", "loading", 
    "scoreloading"), pc_type, pc_num = 5, scale = TRUE, center = TRUE, 
    eigen_loading = "loading", rotate = "none", size_variable = NULL, 
    size_variable_name = NULL, color_variable = NULL, color_variable_name = NULL, 
    shape_variable = NULL, shape_variable_name = NULL, scale_scoreloading = TRUE, 
    first_PC = 1, second_PC = 2, loadings_name = FALSE, loadings_cutvalue = NULL, 
    loadings_cutpercent = NULL, plot = TRUE, option = c("PCA", 
        "WGCNA"), minModuleSize = 2, ...) 
{
    library(reshape2)
    if (missing(option)) {
        option = "PCA"
    }
    if (option == "PCA") {
        plot_object <- OmicsReductionPCA(dataframe, plottype = plottype, 
            pc_type = pc_type, pc_num = pc_num, scale = scale, 
            center = center, eigen_loading = eigen_loading, rotate = rotate, 
            size_variable = size_variable, size_variable_name = size_variable_name, 
            color_variable = color_variable, color_variable_name = color_variable_name, 
            shape_variable = shape_variable, shape_variable_name = shape_variable_name, 
            scale_scoreloading = scale_scoreloading, first_PC = first_PC, 
            second_PC = second_PC, loadings_name = loadings_name, 
            loadings_cutvalue = loadings_cutvalue, loadings_cutpercent = loadings_cutpercent, 
            plot = plot, ...)
    }
    if (option == "WGCNA") {
        plot_object <- OmicsReductionWGCNA(dataframe, plottype = plottype, 
            pc_type = pc_type, pc_num = pc_num, scale = scale, 
            center = center, eigen_loading = eigen_loading, rotate = rotate, 
            size_variable = size_variable, size_variable_name = size_variable_name, 
            color_variable = color_variable, color_variable_name = color_variable_name, 
            shape_variable = shape_variable, shape_variable_name = shape_variable_name, 
            scale_scoreloading = scale_scoreloading, first_PC = first_PC, 
            second_PC = second_PC, loadings_name = loadings_name, 
            loadings_cutvalue = loadings_cutvalue, loadings_cutpercent = loadings_cutpercent, 
            plot = plot, minModuleSize = minModuleSize, ...)
    }
    return(plot_object)
}
