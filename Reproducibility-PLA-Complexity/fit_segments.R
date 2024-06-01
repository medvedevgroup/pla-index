model <- function(x, c, a) {
  c * x^(-a)
}

fit_segments <- function(segments_fn, genome_len){
  lines <- readLines(segments_fn)
  
  data <- read.table(text=lines[1:length(lines)], header = FALSE, col.names = c("x", "y"))
  
  x <- data$x
  y <- data$y
  
  # Fit the model to the data
  fit <- nls(y ~ model(x, c, a), start = list(c = data$y[length(data$y)], a = 1))
  
  param_coef <- coef(fit)
  c <- param_coef["c"]
  a <- param_coef["a"]
  
  c_prime <- c / genome_len
  
  y_pred = predict(fit)
  
  rmse = sqrt(mean((y - y_pred)^2))
  
  summary_result <- summary(fit)
  
  p_values <- summary_result$coefficients[, "Pr(>|t|)"]
  # print(paste(p_values[1], p_values[2]))
  
  # Extract significant values based on a significance level (e.g., 0.05)
  significant_values <- summary_result$coefficients[p_values < 0.05, ]
  
  print(significant_values)
  
  print(paste("Fitted degree (alpha):", a,"Fitted constatnt (Beta): ",c_prime))
}

input_file <- "fit_segment_file.txt"
lines <- readLines(input_file)


segments_fn <- lines[1]
genome_len <- as.numeric(gsub('"', '',lines[2]))
fit_segments(segments_fn, genome_len)
