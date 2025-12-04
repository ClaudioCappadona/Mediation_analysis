scale_by_iqr <- function(x) {
  iqr_value <- IQR(x, na.rm = TRUE)  # Compute IQR
  median_value <- median(x, na.rm = TRUE)  # Compute median
  (x - median_value) / iqr_value  # Center by median and scale by IQR
}
        