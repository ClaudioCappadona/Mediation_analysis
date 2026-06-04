rm(list = ls())

sessionInfo()

library(CMAverse)
library(haven)
library(dplyr)
library(mediation)
library(ggplot2)
library(devtools)
library(WGCNA)
library(purrr)
library(parallel)
library(openxlsx)
library(argparse)


# +
options(stringsAsFactors = F)
parser <- ArgumentParser()


parser$add_argument("--max_cores", type = "character")
parser$add_argument("--outfolder", type="character")
parser$add_argument("--cohort", type="character")
parser$add_argument("--workdir", type="character")

Args <- parser$parse_args()
# -

max_cores <- as.numeric(Args$max_cores)
cohort <- Args$cohort
outfolder <- Args$outfolder
workdir <- Args$workdir

source(paste0(workdir,"OmicsReductionPCA.R"))
source(paste0(workdir,"OmicsReductionWGCNA.R"))
source(paste0(workdir,"OmicsReduction.R"))
source(paste0(workdir,"scale_by_iqr.R"))
source(paste0(workdir,"get_residuals.R"))

# +
#options(repr.matrix.max.cols = 1000, repr.matrix.max.rows = 50)
# -

# ### loading data

load("../CC_TriplotGUI_GenR_image")
colnames(exposuresDF_log_ID)[1] <- "Id_metabolomicsChild"

# ### Covariates: no transformation

# +
#str(covariatesDF_ID)

# +
#covariates
# -

# ### Exposures: scaled by IQR

exposuresDF_ID_IQR <- exposuresDF_ID %>% mutate(across(c(2:ncol(exposuresDF_ID)), scale_by_iqr))

# +
#exposures
# -

# ### Mediators: log transformed and scaled by IQR

omicsDF_ID_log_IQR <- omicsDF_ID %>% 
    mutate(across(c(2:ncol(omicsDF_ID)), log)) %>%
    mutate(across(c(2:ncol(omicsDF_ID)), scale_by_iqr))


# ### outcomes: residualized 

other_pheno_5y <- read_sav("./20251107_Additionalvariables_Chalmers.sav")

# including age and gender. Height is not present, will be added soon.
for_residualiz <- data_1 %>% dplyr::select(Id_metabolomicsChild, GENDER, agey5child)
for_residualiz <- merge(for_residualiz, other_pheno_5y %>% dplyr::select(Id_metabolomicsChild, length5child))

# +
outcomeDF_ID_resid <- merge(outcomeDF_ID, for_residualiz, by="Id_metabolomicsChild")
outcomeDF_log_ID_resid <- merge(for_residualiz, outcomeDF_log_ID, by.y= "outcomeDF_ID$Id_metabolomicsChild", by.x="Id_metabolomicsChild")
#outcomeDF_ID_resid$agey5child2 <- outcomeDF_ID_resid$agey5child * outcomeDF_ID_resid$agey5child

outcomeDF_ID_resid$agey5child_squared <- outcomeDF_ID_resid$agey5child * outcomeDF_ID_resid$agey5child
outcomeDF_ID_resid$INSUCH5_C_log <- log(outcomeDF_ID_resid$INSUCH5_C)
outcomeDF_ID_resid$TGCH5_C_log <- log(outcomeDF_ID_resid$TGCH5_C)

# +
# calculate systolic blood pressure residuals
outcomeDF_ID_resid <- get_residuals(data_original = outcomeDF_ID_resid, 
              id_var = "Id_metabolomicsChild", 
              predictors = c("GENDER", "agey5child", "length5child"), 
              response = "MeanSBP_ex1_5child",
              residual_name =  "MeanSBP_ex1_5child_res")

# calculate diastolic blood pressure residuals
outcomeDF_ID_resid <- get_residuals(data_original = outcomeDF_ID_resid, 
              id_var = "Id_metabolomicsChild", 
              predictors = c("GENDER", "agey5child", "length5child"), 
              response = "MeanDBP_ex1_5child",
              residual_name =  "MeanDBP_ex1_5child_res")

# calculate fat mass residuals
outcomeDF_ID_resid <- get_residuals(data_original = outcomeDF_ID_resid, 
              id_var = "Id_metabolomicsChild", 
              predictors = c("GENDER", "agey5child", "agey5child_squared","length5child"), 
              response = "avg_percent_fat",
              residual_name =  "avg_percent_fat_res")

# calculate Biological parameters residuals
outcomeDF_ID_resid <- get_residuals(data_original = outcomeDF_ID_resid, 
              id_var = "Id_metabolomicsChild", 
              predictors = c("GENDER", "agey5child"), 
              response = "CHOLCH5_C",
              residual_name =  "CHOLCH5_C_res")


outcomeDF_ID_resid <- get_residuals(data_original = outcomeDF_ID_resid, 
              id_var = "Id_metabolomicsChild", 
              predictors = c("GENDER", "agey5child"), 
              response = "HDLCH5_C",
              residual_name =  "HDLCH5_C_res")



outcomeDF_ID_resid <- get_residuals(data_original = outcomeDF_ID_resid, 
              id_var = "Id_metabolomicsChild", 
              predictors = c("GENDER", "agey5child"), 
              response = "HDLCH5_C",
              residual_name =  "HDLCH5_C_res")

outcomeDF_ID_resid <- get_residuals(data_original = outcomeDF_ID_resid, 
              id_var = "Id_metabolomicsChild", 
              predictors = c("GENDER", "agey5child"), 
              response = "TGCH5_C_log",
              residual_name =  "TGCH5_C_log_res")


outcomeDF_ID_resid <- get_residuals(data_original = outcomeDF_ID_resid, 
              id_var = "Id_metabolomicsChild", 
              predictors = c("GENDER", "agey5child"), 
              response = "INSUCH5_C_log",
              residual_name =  "INSUCH5_C_log_res")

# -

outcomes <- grep("*_res$",colnames(outcomeDF_ID_resid), value = T)
outcomes

intermediate_confounders <- c("DIAB_GRA", "HYPERTENSIE_QUEST", "PARITY", "BMI_0")

# ### computing WGCNA modules

ID_var <- "Id_metabolomicsChild"

# +
WGCNA_n_modules <- {}

#foreach(mod_size = c(2:12)) %dopar% {
for(mod_size in c(2:20)){
#log_ver <- "_log"
#mod_size <- 14
deepSplit <- 0
MEDissThres <- 0.4
    
    
enableWGCNAThreads(nThreads = max_cores)
    
reduced_Omics <<- OmicsReduction(dataframe = omicsDF_ID_log_IQR %>% dplyr::select(-ID_var),
                                   plottype = c("scree","score","loading","scoreloading"),
                                   pc_type = "principal", 
# Users can choose freely to use principal() instead of prcomp() to perform PCA.The result should be very similar, but principal() allows more functionalities
                                   first_PC = 1,
                                   second_PC= 2,
                                   option = "WGCNA",
# In this example we use WGCNA to perform data reduction, therefore specifying `option="WGCNA"`.
                                   scale = F,
                                   center = F,
                                   eigen_loading = "loading",
# If you want to show loading or scale it to eigevalues.
                                   rotate = "none",
# This argument is inherited from the principal() function
                                   #size_variable = auxilary2$FastingGlucose,
                                   #size_variable_name = "FastigGlucose",
                                   color_variable = covariatesDF[["GENDER"]],
                                   color_variable_name = "Gender",
                                   #shape_variable = auxilary2$PhysicalActivityIndex,
                                   #shape_variable_name = "PhysicalActivityIndex",
# The size,shape,color variables are used to differentiate points on the score plot and biplot 
                                   loadings_name = F,
# Show the names of loadings
                                   loadings_cutpercent = 0.2,
# The loadings bigger than 20% of all the loadings will show up
                                   minModuleSize = mod_size, 
                                   deepSplit = deepSplit,
                                   MEDissThres = MEDissThres
# The mimimum number of variable in each cluster of WGCNA
)
    
df <- data.frame(
    module_size = mod_size,
    n_modules = reduced_Omics$object$scores %>% as.data.frame() %>% ncol()
    )

WGCNA_n_modules <- rbind(WGCNA_n_modules, df)
}
# -

write.xlsx(x= WGCNA_n_modules, file= paste0(outfolder, cohort, "_WGCNA_n_modules.xlsx"))

pdf(file = paste0(outfolder, cohort, "_WGCNA_n_modules.pdf"), width = 7, height = 5)
    WGCNA_n_modules %>% plot(main = paste(cohort, "WGCNA n modules"))
dev.off()

write.table(x= WGCNA_n_modules[!duplicated(WGCNA_n_modules$n_modules),] %>% filter(n_modules < 10), 
            file= paste0(outfolder, cohort, "_WGCNA_selected_modules.tsv"), 
            quote = F, sep = "\t", row.names = F, col.names = T)

rm (mod_size)

save.image(paste0(outfolder, cohort, "_WGCNA_parameters.Rdata"))
