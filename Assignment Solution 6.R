####################################Assignmnet 6###########################################
#Student: Hao Wei (Uni-Jena), Chanda Albert Jacob Kasoma (Uni Jena), Raphael Kroes (Uni Jena)

####################################Question 1
library(readxl)  
library(dplyr)   
library(lubridate) 

#prt_data <- read_excel("PRT_data.xls", sheet = 2)
mon_data <- read_excel("Dataset_EA-MPD.xlsx", sheet = 4)
keep     <- c("date", "OIS_6M")
mon_data <- mon_data[keep]

##### monetary shock data
mon_data$month <- 0
mon_data$year  <- 0
mon_data$monthyear <- 0

mon_data$month <- format(mon_data$date, "%m")
mon_data$year  <- format(mon_data$date, "%Y")
mon_data$monthyear <- paste (mon_data$month, mon_data$year, sep = "")

tab       <- table(cut(mon_data$date, "month"))
mon_shock <- data.frame(monthyear = format(as.Date(names(tab)), "%m%Y",), shock_freq=as.vector(tab))



####################################Question 2
library(dplyr)
library(lubridate)

bond_data <- read.csv("IRLTLT01DEM156N.csv")
bond_data <- bond_data %>%
  rename(date = observation_date, government_bond = IRLTLT01DEM156N) %>%
  mutate(date = format(ymd(date), "%m%Y"))  


prod_data <- read.csv("Production_Index.csv")
prod_data <- prod_data %>%
  rename(date = TIME_PERIOD, production_index = OBS_VALUE) %>%
  mutate(date = format(ym(date), "%m%Y"))  

# merge
merged_data <- full_join(bond_data, prod_data, by = "date") %>%
  select(date, government_bond, production_index)


cpi_data <- read.csv("CPI.csv")
cpi_data <- cpi_data %>%
  rename(date = TIME_PERIOD, CPI = OBS_VALUE) %>%   
  mutate(date = format(ym(date), "%m%Y"))         


unemployment_data <- read.csv("unemployment_rate.csv")
unemployment_data <- unemployment_data %>%
  rename(date = TIME_PERIOD, unemployment_rate = OBS_VALUE) %>% 
  mutate(date = format(ym(date), "%m%Y"))                      

# Merge
macro_data <- merged_data %>%
  full_join(cpi_data, by = "date") %>%
  full_join(unemployment_data, by = "date")

macro_data <- macro_data %>%
  select(date, government_bond, production_index, CPI, unemployment_rate)

mon_shock <- mon_shock %>%
  rename(date = monthyear)
# merge
final_data <- full_join(macro_data, mon_shock, by = "date")
write.csv(final_data, "final_combined_data.csv", row.names = FALSE)


####################################Question 3
#Distributed Lag Model
library(dplyr)
library(lmtest)
library(sandwich)
library(ggplot2)
library(gridExtra) 

hmax <- 24

for (h in 1:hmax) {
  final_data[[paste0("OIS_6M_", h)]] <- dplyr::lag(final_data$shock_freq, h)
}

estimate_dl_model <- function(dep_var, shock_prefix, hmax, data) {
  formula <- as.formula(paste0(dep_var, " ~ ", 
                               paste0(shock_prefix, 1:hmax, collapse = " + ")))
  model <- lm(formula, data = data)
  
  coef <- coef(model)
  se <- sqrt(diag(vcovHC(model, type = "HC3")))
  results <- data.frame(
    lag = 1:hmax,
    coef = coef[grep(shock_prefix, names(coef))],
    lower = coef[grep(shock_prefix, names(coef))] - 1.64 * se[grep(shock_prefix, names(se))],
    upper = coef[grep(shock_prefix, names(coef))] + 1.64 * se[grep(shock_prefix, names(se))]
  )
  return(results)
}

results_unemp <- estimate_dl_model("unemployment_rate", "OIS_6M_", hmax, final_data)
results_prod <- estimate_dl_model("production_index", "OIS_6M_", hmax, final_data)
results_cpi <- estimate_dl_model("CPI", "OIS_6M_", hmax, final_data)
results_bond <- estimate_dl_model("government_bond", "OIS_6M_", hmax, final_data)

plot_response <- function(results, title, y_label) {
  ggplot(results, aes(x = lag)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "90% Conf. Interval"), alpha = 0.5) +
    geom_line(aes(y = coef, color = "Coefficient"), size = 1) +
    geom_line(aes(y = lower, color = "Lower Bound"), linetype = "dashed") +
    geom_line(aes(y = upper, color = "Upper Bound"), linetype = "dashed") +
    scale_color_manual(values = c("Coefficient" = "blue", "Lower Bound" = "blue", "Upper Bound" = "blue")) +
    scale_fill_manual(values = c("90% Conf. Interval" = "lightblue")) +
    labs(title = title, x = "Lag (in months)", y = y_label, color = "Legend", fill = "Legend") +
    theme_minimal()
}


p1 <- plot_response(results_unemp, "Response Coefficients Over Time, Unemployment Outcome", "Response Coefficient")
p2 <- plot_response(results_prod, "Response Coefficients Over Time, Production Outcome", "Response Coefficient")
p3 <- plot_response(results_cpi, "Response Coefficients Over Time, CPI Outcome", "Response Coefficient")
p4 <- plot_response(results_bond, "Response Coefficients Over Time, 10Y Bond Outcome", "Response Coefficient")


grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

#Local Linear Projections
library(dplyr)
library(lmtest)
library(sandwich)
library(ggplot2)
library(gridExtra) 

hmax <- 24

for (h in 0:hmax) {
  final_data[[paste0("unemployment_rate_lead_", h)]] <- dplyr::lead(final_data$unemployment_rate, h)
  final_data[[paste0("production_index_lead_", h)]] <- dplyr::lead(final_data$production_index, h)
  final_data[[paste0("CPI_lead_", h)]] <- dplyr::lead(final_data$CPI, h)
  final_data[[paste0("government_bond_lead_", h)]] <- dplyr::lead(final_data$government_bond, h)
}


estimate_llp_model <- function(dep_var_lead_prefix, shock_var, hmax, data) {
  results <- data.frame(lag = 0:hmax, coef = NA, lower = NA, upper = NA)
  
  for (h in 0:hmax) {

    lead_var <- paste0(dep_var_lead_prefix, h)
    formula <- as.formula(paste0(lead_var, " ~ ", shock_var))
    
    model <- lm(formula, data = data)

    coef <- coef(model)[shock_var]
    se <- sqrt(diag(vcovHC(model, type = "HC1")))[shock_var]

    results$coef[results$lag == h] <- coef
    results$lower[results$lag == h] <- coef - 1.64 * se
    results$upper[results$lag == h] <- coef + 1.64 * se
  }
  return(results)
}


results_llp_unemp <- estimate_llp_model("unemployment_rate_lead_", "shock_freq", hmax, final_data)
results_llp_prod <- estimate_llp_model("production_index_lead_", "shock_freq", hmax, final_data)
results_llp_cpi <- estimate_llp_model("CPI_lead_", "shock_freq", hmax, final_data)
results_llp_bond <- estimate_llp_model("government_bond_lead_", "shock_freq", hmax, final_data)

plot_response_llp <- function(results, title, y_label) {
  ggplot(results, aes(x = lag, y = coef)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "90% Conf. Interval"), alpha = 0.5) +
    geom_line(aes(color = "Coefficient"), size = 1) +
    geom_line(aes(y = lower, color = "Lower Bound"), linetype = "dashed") +
    geom_line(aes(y = upper, color = "Upper Bound"), linetype = "dashed") +
    scale_color_manual(values = c("Coefficient" = "darkgreen", "Lower Bound" = "darkgreen", "Upper Bound" = "darkgreen")) +
    scale_fill_manual(values = c("90% Conf. Interval" = "lightgreen")) +
    labs(title = title, x = "Horizon (in months)", y = y_label, color = "Legend", fill = "Legend") +
    theme_minimal()
}

p1_llp <- plot_response_llp(results_llp_unemp, "Local Linear Projection: Unemployment Rate", "Response Coefficient")
p2_llp <- plot_response_llp(results_llp_prod, "Local Linear Projection: Production Index", "Response Coefficient")
p3_llp <- plot_response_llp(results_llp_cpi, "Local Linear Projection: CPI", "Response Coefficient")
p4_llp <- plot_response_llp(results_llp_bond, "Local Linear Projection: 10Y Bond", "Response Coefficient")

grid.arrange(p1_llp, p2_llp, p3_llp, p4_llp, ncol = 2, nrow = 2)

####################################Question 4
#Distributed Lag Model
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(sandwich)
library(lmtest)
install.packages("patchwork")
library(patchwork)

# Parameters
horizons <- c(12, 36)  # Response horizons
hmax <- max(horizons)  # Maximum lag

# Generate lagged variables
final_data <- final_data %>%
  mutate(
    UNEMP_L1 = lag(unemployment_rate, 1),
    LOG_PROD_L1 = lag(production_index, 1),
    LOG_CPI_L1 = lag(CPI, 1),
    Y_BOND_L1 = lag(government_bond, 1)
  )

for (h in 1:hmax) {
  final_data[[paste0("OIS_6M_", h)]] <- dplyr::lag(final_data$shock_freq, h)
}

final_data_clean <- final_data %>%
  filter(complete.cases(across(starts_with("OIS_6M_"))))

# Distributed Lag Model function
estimate_dl_model_robust <- function(dep_var, shock_prefix, horizon, data) {
  formula <- as.formula(
    paste0(dep_var, " ~ ", paste0(shock_prefix, 1:horizon, collapse = " + "),
           " + UNEMP_L1 + LOG_PROD_L1 + LOG_CPI_L1 + Y_BOND_L1")
  )
  model <- lm(formula, data = data)
  coef <- coef(model)
  se <- sqrt(diag(vcovHC(model, type = "HC1")))
  
  results <- data.frame(
    lag = 1:horizon,
    coef = coef[grep(shock_prefix, names(coef))],
    lower = coef[grep(shock_prefix, names(coef))] - 1.64 * se[grep(shock_prefix, names(se))],
    upper = coef[grep(shock_prefix, names(coef))] + 1.64 * se[grep(shock_prefix, names(se))]
  )
  return(results)
}

# Create and combine plots
for (H in horizons) {
  cat("\nResponse Horizon: ", H, " months\n")
  
  results_unemp <- estimate_dl_model_robust("unemployment_rate", "OIS_6M_", H, final_data_clean)
  results_prod <- estimate_dl_model_robust("production_index", "OIS_6M_", H, final_data_clean)
  results_cpi <- estimate_dl_model_robust("CPI", "OIS_6M_", H, final_data_clean)
  results_bond <- estimate_dl_model_robust("government_bond", "OIS_6M_", H, final_data_clean)
  
  # Plot for unemployment
  plot_unemp <- ggplot(results_unemp, aes(x = lag, y = coef)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "navy", alpha = 0.2) +
    geom_line(color = "navy") +
    labs(title = "Unemployment Rate", x = "Lag (in months)", y = "Response Coefficient") +
    theme_minimal()
  
  # Plot for production index
  plot_prod <- ggplot(results_prod, aes(x = lag, y = coef)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "navy", alpha = 0.2) +
    geom_line(color = "navy") +
    labs(title = "Production Index", x = "Lag (in months)", y = "Response Coefficient") +
    theme_minimal()
  
  # Plot for CPI
  plot_cpi <- ggplot(results_cpi, aes(x = lag, y = coef)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "navy", alpha = 0.2) +
    geom_line(color = "navy") +
    labs(title = "CPI", x = "Lag (in months)", y = "Response Coefficient") +
    theme_minimal()
  
  # Plot for government bond
  plot_bond <- ggplot(results_bond, aes(x = lag, y = coef)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "navy", alpha = 0.2) +
    geom_line(color = "navy") +
    labs(title = "10Y Bond Yield", x = "Lag (in months)", y = "Response Coefficient") +
    theme_minimal()
  
  # Combine all four plots
  combined_plot <- (plot_unemp | plot_prod) / (plot_cpi | plot_bond) +
    plot_annotation(title = paste("Distributed Lag Model: Response Horizon =", H, "Months"))
  
  # Print and save the combined plot
  print(combined_plot)
  ggsave(paste0("DLM_Response_H", H, ".pdf"), combined_plot, width = 12, height = 8)
}


# Local Linear Projections

# Load required libraries
library(dplyr)
library(ggplot2)
library(sandwich)
library(lmtest)
library(patchwork)

# Set parameters
horizons <- c(12, 36)  # Response horizons
endogenous_vars <- c("unemployment_rate", "production_index", "CPI", "government_bond")

# Ensure all lagged variables are created
final_data <- final_data %>%
  mutate(
    UNEMP_L1 = lag(unemployment_rate, 1),
    LOG_PROD_L1 = lag(production_index, 1),
    LOG_CPI_L1 = lag(CPI, 1),
    Y_BOND_L1 = lag(government_bond, 1),
    shock_freq_L1 = lag(shock_freq, 1)  # Ensure this exists
  ) %>%
  na.omit()  # Drop rows with missing values

# Function to create lead variables
create_lead_vars <- function(data, var, horizon) {
  for (h in 0:horizon) {
    lead_col <- paste0(var, "_lead_", h)
    data[[lead_col]] <- dplyr::lead(data[[var]], h)
  }
  return(data)
}

# Function to estimate Local Linear Projections
estimate_llp_model_robust <- function(dep_var, shock_var, horizon, data) {
  results <- data.frame(lag = 0:horizon, coef = NA, lower = NA, upper = NA)
  
  # Loop over horizons to estimate the model
  for (h in 0:horizon) {
    lead_col <- paste0(dep_var, "_lead_", h)
    formula <- as.formula(paste0(lead_col, " ~ ", shock_var, " + ", 
                                 "shock_freq_L1 + UNEMP_L1 + LOG_PROD_L1 + LOG_CPI_L1 + Y_BOND_L1"))
    model <- lm(formula, data = data)
    se <- sqrt(diag(vcovHC(model, type = "HC1")))
    coef_val <- coef(model)[shock_var]
    results$coef[h + 1] <- coef_val
    results$lower[h + 1] <- coef_val - 1.64 * se[shock_var]
    results$upper[h + 1] <- coef_val + 1.64 * se[shock_var]
  }
  return(results)
}

# Create lead variables for each endogenous variable
for (var in endogenous_vars) {
  final_data <- create_lead_vars(final_data, var, max(horizons))
}

# Loop over horizons and plot results
for (H in horizons) {
  cat("\nResponse Horizon: ", H, " months\n")
  
  # Unemployment Rate
  results_unemp <- estimate_llp_model_robust("unemployment_rate", "shock_freq", H, final_data)
  plot_unemp <- ggplot(results_unemp, aes(x = lag, y = coef)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "darkgreen", alpha = 0.2) +
    geom_line(color = "darkgreen") +
    labs(title = "Unemployment Rate", x = "Horizon (in months)", y = "Response Coefficient") +
    theme_minimal()
  
  # Production Index
  results_prod <- estimate_llp_model_robust("production_index", "shock_freq", H, final_data)
  plot_prod <- ggplot(results_prod, aes(x = lag, y = coef)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "darkgreen", alpha = 0.2) +
    geom_line(color = "darkgreen") +
    labs(title = "Production Index", x = "Horizon (in months)", y = "Response Coefficient") +
    theme_minimal()
  
  # CPI
  results_cpi <- estimate_llp_model_robust("CPI", "shock_freq", H, final_data)
  plot_cpi <- ggplot(results_cpi, aes(x = lag, y = coef)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "darkgreen", alpha = 0.2) +
    geom_line(color = "darkgreen") +
    labs(title = "CPI", x = "Horizon (in months)", y = "Response Coefficient") +
    theme_minimal()
  
  # Government Bond Yield
  results_bond <- estimate_llp_model_robust("government_bond", "shock_freq", H, final_data)
  plot_bond <- ggplot(results_bond, aes(x = lag, y = coef)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "darkgreen", alpha = 0.2) +
    geom_line(color = "darkgreen") +
    labs(title = "10Y Bond Yield", x = "Horizon (in months)", y = "Response Coefficient") +
    theme_minimal()
  
  # Combine and display plots
  combined_plot <- (plot_unemp | plot_prod) / (plot_cpi | plot_bond) +
    plot_annotation(title = paste("Local Linear Projections: Response Horizon =", H, "Months"))
  
  print(combined_plot)
}
