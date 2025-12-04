rm(list = ls())

sessionInfo()

# ### loading data and parameters

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

parser$add_argument("--mod_size", type = "character")
parser$add_argument("--max_cores", type = "character")
parser$add_argument("--outfolder", type="character")
parser$add_argument("--cohort", type="character")
parser$add_argument("--workdir", type="character")
parser$add_argument("--treatment", type="character")

Args <- parser$parse_args()

# +

outfolder <- Args$outfolder
cohort <- Args$cohort
workdir <- Args$workdir
treatment <- Args$treatment
# -

load(paste0(outfolder, cohort, "_WGCNA_parameters.Rdata"))

source(paste0(workdir,"OmicsReductionPCA.R"))
source(paste0(workdir,"OmicsReductionWGCNA.R"))
source(paste0(workdir,"OmicsReduction.R"))
source(paste0(workdir,"scale_by_iqr.R"))
source(paste0(workdir,"get_residuals.R"))

# +
#options(repr.matrix.max.cols = 1000, repr.matrix.max.rows = 50)

# +
options(stringsAsFactors = F)
parser <- ArgumentParser()

parser$add_argument("--mod_size", type = "character")
parser$add_argument("--max_cores", type = "character")
parser$add_argument("--outfolder", type="character")
parser$add_argument("--cohort", type="character")
parser$add_argument("--workdir", type="character")
parser$add_argument("--treatment", type="character")


Args <- parser$parse_args()
# -

max_cores <- as.numeric(Args$max_cores)
mod_size <- as.numeric(Args$mod_size)
cohort <- Args$cohort
outfolder <- Args$outfolder
workdir <- Args$workdir
treatment <- Args$treatment

# ### General info

print(paste0("processing ",cohort," cohort - module size ", mod_size))
print("")

# ### computing WGCNA modules

ID_var <- "Id_metabolomicsChild"

# +

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
                                   color_variable = covariatesDF_ID[["GENDER"]],
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

# -

metabol_cluster <- reduced_Omics$object$scores %>% as.data.frame()

metabol_cluster <- cbind(omicsDF_ID_log_IQR %>% dplyr::select(all_of(ID_var)),metabol_cluster) 

# ### CMAverse - intermediate confounders 

# +

#library(doParallel)
#registerDoParallel(cores = parallel::detectCores() - 1)
#cl <- 15
#registerDoParallel(cores = cl)

CMAverse_res <- {}


#CMAverse_analysis <- function(){

#combinations <- expand.grid(
#    treatment = c("minVSmax"),
#    exposure = colnames(exposuresDF)[1],
#    module = colnames(metabol_cluster),
#    outcome = c("MeanSBP_ex1_5child_res"),
#    stringsAsFactors = FALSE)

#CMAverse_analysis_res <- mclapply(1:nrow(combinations), function(i) {
    
#row <- combinations[i,]

#CMAverse_analysis_res <- apply(combinations,1, function(row) {
#    treatment = row$treatment
#    exposure = row$exposure
#    module = row$module
#    outcome = row$outcome
    #print(paste(exposure, module, outcome))

#CMAverse_analysis_res <- foreach(i = 1:nrow(combinations), .combine = 'rbind') %dopar% {
    #row <- combinations[i, ]
    #treatment = row$treatment
    #exposure = row$exposure
    #module = row$module
    #outcome = row$outcome
    #.packages = c("CMAverse")
    #.export = c("metabol_cluster", "exposuresDF")



#foreach(
#       treatment = c("minVSmax"),
#       exposure = colnames(exposuresDF)[1],
#       module = colnames(metabol_cluster),
#       outcome = c("MeanSBP_ex1_5child_res"),
#       .combine = 'rbind',
#       .packages = c("CMAverse"),
#       .export = c("metabol_cluster", "exposuresDF")
#       ) %do% {
    
    
#for (treatment in c("minVSmax", "quantile25vs75")){ #, "quantile10vs90")){
    
    #foreach(exposure = colnames(exposuresDF)[1]) %dopar% {
    for(exposure in colnames(exposuresDF)){
        
        #foreach(module = colnames(metabol_cluster)) %dopar% {
        for(module in colnames(metabol_cluster %>% dplyr::select(-ID_var))){
            
            #foreach(residualiz_outcome = c("MeanSBP_ex1_5child_res")) %dopar% {
            for(outcome in outcomes){
                
                analys_df <- covariatesDF_ID %>%
                merge(exposuresDF_ID_IQR %>% dplyr::select(all_of(ID_var), all_of(exposure))) %>%
                merge(metabol_cluster) %>% 
                merge(outcomeDF_ID_resid %>% dplyr::select(all_of(ID_var), all_of(outcome))) %>%
                merge(data_1 %>% dplyr::select(all_of(ID_var), all_of(intermediate_confounders)))
                
                analys_df_complete <- analys_df[complete.cases(analys_df),]
                       
                #analys_df_complete <- analys_df_complete[1:10,]
                
                for(L_matrix in c("conditional","not_conditional")){


if(treatment == "quantile25vs75"){
    a0 <- quantile(analys_df_complete %>% dplyr::select(exposure), 0.25, na.rm = TRUE)
    a1 <- quantile(analys_df_complete %>% dplyr::select(exposure), 0.75, na.rm = TRUE)
}
                
if(treatment == "quantile10vs90"){
    a0 <- quantile(analys_df_complete %>% dplyr::select(exposure), 0.10, na.rm = TRUE)
    a1 <- quantile(analys_df_complete %>% dplyr::select(exposure), 0.90, na.rm = TRUE)
}
        
if(treatment == "minVSmax"){
    a0 <- min(analys_df_complete %>% dplyr::select(exposure), na.rm = TRUE)
    a1 <- max(analys_df_complete %>% dplyr::select(exposure), na.rm = TRUE)
}    

if(conditional == "conditional"){
    L_matrix <- c(intermediate_confounders,  setdiff(colnames(metabol_cluster), module) %>% as.vector())
    L_matrix_reg <- list("logistic","logistic","ordinal","linear", rep("linear", length(setdiff(colnames(metabol_cluster), module)))) %>% flatten()
} else if(conditional == "not_conditional"){
    L_matrix <- c(intermediate_confounders)
    L_matrix_reg <-list("logistic","logistic","ordinal","linear")
}
                
n_sim <- 10          
#n_sim <- 101
max_iterations <- 5
iteration <- 1

res_gform <- cmest(
  data = analys_df_complete,
  model = "gformula",
  outcome = outcome,
  exposure = exposure,
  estimation = "imputation",
  mediator = module,
  basec = covariates,
  #basec = c("GENDER", "season", "EDUCM", "AGE_M_v2", "SMOKE_ALL", "ETHNMv2", "DELIVERY_MODE", "DIAB_GRA", "HYPERTENSIE_QUEST", "PARITY", "WEIGHT_0"),  
  #intermediate = c("DIAB_GRA", "HYPERTENSIE_QUEST", "PARITY", "WEIGHT_0"),  
  postc = L_matrix,
  mreg = list("linear"),
  postcreg = L_matrix_reg, #consider different distrib for BMI_0 since it's not linear
  yreg = "linear",
  astar = a0,
  a = a1,
  mval = list(0), # reference value for cluster mediators
  #mval = list(max(metabol_cluster)),
  EMint = FALSE, #set true for sensitivity analysis
  inference = "bootstrap",
  #boot = TRUE,
  nboot = n_sim,
  boot.ci.type = "bca"
)

mediation_summary <- summary(res_gform)

                
while(sum(mediation_summary$summarydf$P.val == 0) > 0 && 
   iteration < max_iterations) {
    message("Cannot find non-zero p-value, running iteration n.", iteration, "...\n")
    n_sim <- n_sim * 10
    iteration <- iteration + 1
    
      res_gform <- cmest(
      data = analys_df_complete,
      model = "gformula",
      outcome = outcome,
      exposure = exposure,
      estimation = "imputation",
      mediator = module,
      basec = covariates,
      postc = L_matrix, #change to BMI_0 or later if missing, see notes
      mreg = list("linear"),
      postcreg = L_matrix_reg, #consider different distrib for BMI_0 since it's not linear
      yreg = "linear",
      astar = a0,
      a = a1,
      mval = list(0), # reference value for cluster mediators
      EMint = FALSE, #set true for sensitivity analysis
      inference = "bootstrap",
      nboot = n_sim,
      boot.ci.type = "bca"
      )
    mediation_summary <- summary(res_gform)
  } 


df <- data.frame(
  Package = "CMAverse with intermediate",
  Conditional = L_matrix,
  Exposure = exposure,
  Treatment = treatment,
  Cluster = module,
  Outcome = outcome,
  Effect = row.names(mediation_summary$summarydf),
  Estimate = mediation_summary$summarydf$Estimate,
  Std.error = mediation_summary$summarydf$Std.error,
  CI.Lower = mediation_summary$summarydf$`95% CIL`,
  CI.Upper = mediation_summary$summarydf$`95% CIU`,
  p.value = mediation_summary$summarydf$P.val
)

CMAverse_res <- rbind(CMAverse_res, df)
    
   }#, mc.cores = 20)
  #}
 }
}
}


#CMAverse_df <<- append(CMAverse_df, df)

    #return(df)
                
            #}
    #return(CMAverse_df) 
        #}
    #}
#}
   
#}

#CMAverse_analysis_res <- CMAverse_analysis()
#CMAverse_analysis_df <- CMAverse_analysis_res  %>% do.call(rbind, .)

#CMA_df <- do.call(rbind, CMA_results)

#stopCluster(cl)                  
#registerDoSEQ()  


# +
#CMAverse_res %>% 
#    filter(p.value < 0.05) %>%
#    arrange(p.value) %>% 
#    filter(Effect=="rpnie") 

# +

write.xlsx(x = CMAverse_res, file = paste0(outfolder, cohort,"_",treatment, "_CMAverse_intermediate_outcomes_res_modsize",mod_size,".xlsx"))
write.table(x = CMAverse_res, file = paste0(outfolder, cohort"_",treatment, "_CMAverse_intermediate_outcomes_res_modsize",mod_size,".tsv"), 
                quote = F, sep = "\t", row.names = F, col.names = T)
# -


