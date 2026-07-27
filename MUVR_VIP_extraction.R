sessionInfo() #metabolomics

library(dplyr)
library(ggplot2)
library(ggrepel)
library(MUVR)
library(gtools)

# report function
report <- function(models) {
  nMod <- length(models)
  modelName <- type <- character(nMod)
  nVar <- q2 <- ber <- fitness <- matrix(NA, nrow = nMod, ncol = 3)
  colnames(nVar) <- paste0('nVar_', c('min', 'mid', 'max'))
  colnames(q2) <- paste0('q2_', c('min', 'mid', 'max'))
  colnames(ber) <- paste0('ber_', c('min', 'mid', 'max'))
  colnames(fitness) <- paste0('fitness_', c('min', 'mid', 'max'))
  i <- 0
  for (model in models) {
    i <- i + 1
    modelName[i] <- basename(model)
    cat('--------------------------------\n')
    cat('model:', modelName[i], '\n')
    
    if(endsWith(model,"rda")){
        load(file = model)
    } else if(endsWith(model,"rds")){
        mod <- readRDS(model)
    }
    
    nVar[i, ] <- mod$nVar
    if(any(class(mod) == 'Regression')) {
      type[i] <- 'R'
      q2[i,] <- fitness[i,] <- round(mod$fitMetric$Q2, 2)
      cat(paste0('Q2: ', paste(q2[i,], collapse = ' ')))
    }
    if(any(class(mod) == 'Classification')) {
      type[i] <- 'C'
      ber[i,] <- apply(mod$yClass, 2, function(x) MUVR:::getBER(mod$inData$Y, x))
      nClass <- length(unique(mod$inData$Y))
      fitness[i,] <- round((1 - ber[i,] / (1 - (1 / nClass))), 2)
      ber[i,] <- round(ber[i,], 2)
      cat(paste0('CR: ', paste(1 - round(mod$miss / length(mod$inData$Y), 2), collapse = ' ')))
      cat(paste0('\nBER: ', paste(ber[i,], collapse = ' ')))
      cat(paste0('\nfit: ', paste(fitness[i,], collapse = ' ')))
    }
    cat('\n')
  }
  resultDF <- data.frame(model = modelName, type, nVar, q2, ber, fitness)
  resultDF <- resultDF[order(resultDF$fitness_min, decreasing = TRUE),]
  return(resultDF)
}

# +
# define path to find MUVR models, and outpath to export the VIPs
path <- '/path/to/MUVR/mods/'

outpath <- "./"
# -

mods <-  list.files(path = path, pattern = "rda",recursive = T, full.names = TRUE)

mods

#different for each cohort
#mods <- grep("elapse|surf|O3",mods, value = T,invert = T)

report <- report(mods)

report <- report %>% 
select(model,q2_max,nVar_max) %>% mutate(model_name = recode(model,
                          "stringent_cohort_birth_73_pm10_preg.rda" = "PM10",
                          "stringent_cohort_birth_74_pm25_preg.rda" = "PM2.5",                       
                          "stringent_cohort_birth_132_pm25abs_preg.rda" = "Black carbon",                     
                          "stringent_cohort_birth_72_no2_preg.rda" = "NO2",
                          "stringent_cohort_birth_125_hum_preg.rda" = "Humidity",
                          "stringent_cohort_birth_122_tm_preg.rda" = "Temperature",
                          "stringent_cohort_birth_119_uvdvc_preg.rda" = "UV"))

report

path

model_VIPs <- list()
for(i in 1:nrow(report)){
    load(paste0(path,report$model[i])) # new object will be called "mod"
    
    model_VIPs[[report$model_name[i]]] <- MUVR::getVIP(mod,model = "max")
}

model_VIPs

saveRDS(model_VIPs, paste0(outpath, "VIPs_stringent_models_test.rds"))


