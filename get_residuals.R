# Function for linear regression with customizable predictors
get_residuals <- function(data_original,data_cc_response,id_var, predictors, response, residual_name) {
  
  #subset database on cc fot response to run the lm model
  #data_cc_response <-data_original[complete.cases(data_original[[response]]), ]
    
  data_cc_response <-data_original
  
  # Formulate the regression formula
  formula <- as.formula(paste(response, "~", paste(predictors, collapse = "+")))
  
  # Fit the linear regression model
  model <- lm(formula, data = data_cc_response, na.action = na.exclude)
  
  # Extract residuals from the model
  data_cc_response[[residual_name]] <-rstandard(model)
  
  residual_name<-subset(data_cc_response, select=c(id_var, residual_name))
  
  # Return the data frame with ONLY ID and residuals
  #return(residual_name)
  return(data_cc_response)
}
# keep a data frame with only ID and residuals