# ============================================================
# PORTFOLIO OPTIMIZATION & BACKTESTING ANALYSIS (2004-2025)
# Retrospective Analysis: 2025 Market Leaders Backtested to 2004
# ============================================================
#
# CENTRAL QUESTION:
# What if you had known in 2004 which five companies would dominate 
# in 2025? This retrospective analysis tests portfolio optimization 
# strategies on MSFT, NVDA, AAPL, AMZN, BRK-A using 21 years of 
# actual data—a what-if scenario with perfect hindsight.
#
# WHY THIS MATTERS:
# Tests three key hypotheses: (1) Concentration benefit of mega-caps 
# vs S&P 500, (2) Added value of optimization vs equal weighting, 
# (3) Impact of transaction costs.
#
# DATA & METHODOLOGY:
# Period: Nov 2004 - Nov 2025 (21 years) | Training: 48-month 
# inception window | Rebalancing: Annual
# Stocks: MSFT, NVDA, AAPL, AMZN, BRK-A (top 5 market cap, Nov 2025)
# Optimization: CAPM + Global Minimum Variance (GMVP) + 
# Mean-Variance efficient frontier
# Backtesting: Rolling-window, no look-ahead bias
# Transaction costs: 10 bps per rebalancing
# Inflation: CPI-adjusted to Nov 2004 dollars (Nov 2004 - Nov 2024)
#
# KEY OUTPUTS:
# CAPM betas | Annual metrics (returns, volatility, Sharpe) | 
# Wealth evolution charts | Performance comparison
#
# DISCLAIMER:
# Retrospective analysis with perfect hindsight. Does NOT demonstrate 
# predictive ability. Real 2004 investors faced genuine uncertainty 
# about these companies' futures.
#
# ============================================================
rm(list = ls())
# ============================================================
# SECTION 1: ENVIRONMENT SETUP
# ============================================================

packages <- c(
  "quantmod",           # Financial data retrieval
  "dplyr",              # Data manipulation
  "tidyr",              # Data reshaping
  "zoo",                # Time series handling
  "ggplot2",            # Visualization
  "PerformanceAnalytics",# Financial metrics
  "lubridate"           # Date operations
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# SECTION 2: UTILITY FUNCTIONS
# ============================================================

#' Load Monthly Stock Returns
#'
#' @param ticker Character string for stock ticker
#' @param from Start date (YYYY-MM-DD format)
#' @param to End date (YYYY-MM-DD format)
#'
#' @return Data frame with OHLCV data and computed returns
#'
#' @details Retrieves historical pricing data and calculates
#' monthly percentage returns, excluding the first NA observation.

load_stock_monthly <- function(ticker, from = "1999-11-01", to = "2025-11-02") {
  
  data <- getSymbols(
    ticker,
    from = from,
    to = to,
    periodicity = "monthly",
    auto.assign = FALSE
  )
  
  df <- as.data.frame(data)
  colnames(df) <- c("open", "high", "low", "close", "volume", "adjusted_close")
  df$date <- as.Date(rownames(df))
  rownames(df) <- NULL
  
  df <- df %>%
    select(date, open, high, low, close, adjusted_close, volume) %>%
    mutate(
      returns = (adjusted_close / lag(adjusted_close) - 1) * 100,
      compound_growth = adjusted_close / adjusted_close[1]
    )
  
  return(df)
}

#' Geometric Mean of Returns
#'
#' @param x Numeric vector of returns (in %)
#'
#' @return Single numeric value representing geometric mean return (%)
#'
#' @details Computes compound annual growth rate (CAGR) from monthly
#' returns, properly handling NA values.

geometric_mean_returns <- function(x) {
  x_clean <- na.omit(x)
  n <- length(x_clean)
  ((prod(1 + x_clean / 100)) ^ (1 / n) - 1)
}
geometric_mean_rf <- function(rf_rates) {
  rf_clean <- na.omit(rf_rates)
  n <- length(rf_clean)
  if (n == 0) return(NA_real_)
  ((prod(1 + rf_clean / 100)) ^ (1 / n) - 1) * 100
}

#' Global Minimum Variance Portfolio (GMVP)
#'
#' @param window_data Matrix of returns
#'
#' @return List containing optimal weights and portfolio volatility
#'
#' @details Solves for the minimum-variance portfolio using the
#' closed-form solution: w = Σ⁻¹·1 / (1ᵀ·Σ⁻¹·1)

calculate_gmvp_weights <- function(window_data, min_weight = 0, max_weight = 1) {
  
  library(quadprog)
  
  sigma_mat <- cov(window_data)
  n_assets <- ncol(window_data)
  
  # Paramètres pour quadprog
  D <- 2 * sigma_mat
  d <- rep(0, n_assets)
  
  # Contrainte d'égalité: sum(w) = 1
  A_eq <- matrix(1, nrow = 1, ncol = n_assets)
  b_eq <- 1
  
  # Contraintes d'inégalité: min_weight <= w_i <= max_weight
  A_in <- rbind(
    diag(n_assets),                    # w_i >= min_weight
    -diag(n_assets)                    # w_i <= max_weight
  )
  b_in <- c(rep(min_weight, n_assets), rep(-max_weight, n_assets))
  
  # Résoudre le problème d'optimisation quadratique
  result <- quadprog::solve.QP(
    Dmat = D,
    dvec = d,
    Amat = t(rbind(A_eq, A_in)),
    bvec = c(b_eq, b_in),
    meq = 1
  )
  
  gmvp_weights <- result$solution
  gmvp_var <- as.numeric(t(gmvp_weights) %*% sigma_mat %*% gmvp_weights)
  gmvp_vol <- sqrt(gmvp_var)
  
  return(list(weights = gmvp_weights, volatility = gmvp_vol))
}

calculate_mv_portfolio <- function(window_data, 
                                   expected_returns,
                                   min_weight = 0, 
                                   max_weight = 1) {
  
  library(quadprog)
  
  sigma_mat <- cov(window_data)
  n_assets <- ncol(window_data)
  
  # Rendement cible = MOYENNE CAPM
  target_return <- mean(expected_returns)
  
  D <- 2 * sigma_mat
  d <- rep(0, n_assets)
  
  # DEUX contraintes d'égalité
  A_eq <- rbind(
    matrix(1, nrow = 1, ncol = n_assets),      # sum(w) = 1
    matrix(expected_returns, nrow = 1)         # E[R]^T w = target
  )
  
  b_eq <- c(1, target_return)
  
  # Contraintes d'inégalité (poids min/max)
  A_in <- rbind(
    diag(n_assets),
    -diag(n_assets)
  )
  
  b_in <- c(
    rep(min_weight, n_assets),
    rep(-max_weight, n_assets)
  )
  
  # Résoudre
  result <- quadprog::solve.QP(
    Dmat = D,
    dvec = d,
    Amat = t(rbind(A_eq, A_in)),
    bvec = c(b_eq, b_in),
    meq = 2  # 2 égalités
  )
  
  mv_weights <- result$solution
  mv_var <- as.numeric(t(mv_weights) %*% sigma_mat %*% mv_weights)
  mv_vol <- sqrt(mv_var)
  mv_ret <- as.numeric(t(mv_weights) %*% expected_returns)
  
  return(list(
    weights = mv_weights,
    volatility = mv_vol,
    expected_return = mv_ret,
    target = target_return
  ))
}


# ============================================================
# SECTION 3: DATA ACQUISITION & PREPARATION
# ============================================================

cat("\n--- DATA ACQUISITION ---\n")

# Load individual stocks
tickers <- c("NVDA", "MSFT", "AAPL", "AMZN", "BRK-A")

cat("Loading monthly stock data for:", paste(tickers, collapse = ", "), "\n")

stock_list <- lapply(tickers, load_stock_monthly)
names(stock_list) <- tickers

# Consolidate into single dataframe
stock_list_tagged <- lapply(seq_along(stock_list), function(i) {
  stock_list[[i]] %>% mutate(stock = names(stock_list)[i])
})

data_all <- bind_rows(stock_list_tagged)

# Load S&P 500 benchmark
cat("Loading S&P 500 (SPY) benchmark data...\n")

spy_data <- load_stock_monthly("SPY")
names(spy_data)[names(spy_data) == "adjusted_close"] <- "spy_close"
names(spy_data)[names(spy_data) == "returns"] <- "spy_returns"
names(spy_data)[names(spy_data) == "compound_growth"] <- "spy_compound_growth"

# Load risk-free rate (3-month US Treasury)
cat("Loading risk-free rate (3M Treasury)...\n")

rf_data <- getSymbols(
  "DGS3MO",
  src = "FRED",
  from = "1999-11-01",
  to = "2025-11-02",
  auto.assign = FALSE
) %>%
  as.data.frame()

rf_data$date <- as.Date(rownames(rf_data))
colnames(rf_data) <- c("rf_rate", "date")
rownames(rf_data) <- NULL

rf_data <- rf_data %>%
  mutate(
    rf_rate = as.numeric(as.character(rf_rate)),
    yearmonth = format(date, "%Y-%m")
  ) %>%
  group_by(yearmonth) %>%
  filter(date == min(date)) %>%
  ungroup() %>%
  select(date, rf_rate)

# Handle missing November 2025 data
last_date <- max(rf_data$date)
target_date <- as.Date("2025-11-01")

if (last_date < target_date) {
  new_row <- data.frame(date = target_date, rf_rate = NA)
  rf_data <- bind_rows(rf_data, new_row)
}

# Load Consumer Price Index for inflation adjustment
cat("Loading CPI data for inflation adjustment (through Nov 2024)...\n")

cpi_data <- getSymbols(
  "CPIAUCSL",
  src = "FRED",
  from = "2004-11-01",
  to = "2024-11-01",
  auto.assign = FALSE
) %>%
  as.data.frame()

cpi_data$date <- as.Date(rownames(cpi_data))
colnames(cpi_data) <- c("cpi", "date")
rownames(cpi_data) <- NULL

# ============================================================
# SECTION 4: DATA STRUCTURING
# ============================================================

cat("Structuring return matrices...\n")

# Create wide-format returns matrix
data_wide <- data_all %>%
  select(date, stock, returns) %>%
  pivot_wider(names_from = stock, values_from = returns) %>%
  arrange(date) %>%
  slice(-1)  # Remove first NA

returns_matrix_full <- data_wide %>%
  select(-date) %>%
  as.matrix()

# Subset for 150-month period analysis
n_150 <- 150
returns_matrix_150 <- returns_matrix_full[1:n_150, ]

# Align S&P 500 returns
spy_returns_full <- spy_data$spy_returns[-1]
spy_returns_150 <- spy_returns_full[1:n_150]

cat("✓ Data loaded and structured\n")
cat("  Full period: ", nrow(returns_matrix_full), " months\n")
cat("  Analysis begins: Nov 2004 | Analysis ends: Nov 2025\n\n")

# ============================================================
# SECTION 5: CAPM ANALYSIS
# ============================================================

cat("--- CAPM ANALYSIS ---\n")

# Calculate betas (systematic risk)
betas_full <- apply(returns_matrix_full, 2, function(x) {
  cov(x, spy_returns_full, use = "complete.obs") /
    var(spy_returns_full, na.rm = TRUE)
})

betas_150 <- apply(returns_matrix_150, 2, function(x) {
  cov(x, spy_returns_150, use = "complete.obs") /
    var(spy_returns_150, na.rm = TRUE)
})

# Calculate annualized S&P 500 expected return
mean_spy_monthly_full <- geometric_mean_returns(spy_returns_full)
mean_spy_annual_full <- ((1 + mean_spy_monthly_full) ^ 12 - 1) * 100

mean_spy_monthly_150 <- geometric_mean_returns(spy_returns_150)
mean_spy_annual_150 <- ((1 + mean_spy_monthly_150) ^ 12 - 1) * 100

# Risk-free rates (average over periods)
mean_rf_full <- geometric_mean_rf(rf_data$rf_rate)
mean_rf_150 <- geometric_mean_rf(rf_data$rf_rate[1:150])

# CAPM: E(Ri) = Rf + βi × (E(Rm) - Rf)
expected_returns_full <- mean_rf_full + betas_full * (mean_spy_annual_full - mean_rf_full)
expected_returns_150 <- mean_rf_150 + betas_150 * (mean_spy_annual_150 - mean_rf_150)

capm_summary <- data.frame(
  Stock = tickers,
  Beta_Full = round(betas_full, 3),
  Beta_150 = round(betas_150, 3),
  Expected_Return_Full_Pct = round(expected_returns_full, 2),
  Expected_Return_150mo_Pct = round(expected_returns_150, 2)
)

cat("\nCAPM Results (vs S&P 500 Benchmark):\n")
print(capm_summary)
cat("\n")

# ============================================================
# SECTION 6: ROLLING-WINDOW PORTFOLIO OPTIMIZATION
# ============================================================

cat("--- ROLLING-WINDOW PORTFOLIO OPTIMIZATION ---\n")

initial_periods <- 48  # 4-year inception period
periods_per_year <- 12
n_years <- 22
transaction_cost_rate <- 0.001  # 10 basis points

cat("Parameters:\n")
cat("  Inception period: ", initial_periods, " months (4 years)\n")
cat("  Rebalancing: Annual\n")
cat("  Transaction cost: ", transaction_cost_rate * 100, "%\n")
cat("  Analysis horizon: ", n_years, " years\n\n")

optimization_results <- list()

for (year in 1:n_years) {
  
  end_period <- initial_periods + (year * periods_per_year)
  
  if (end_period > nrow(returns_matrix_full)) {
    break
  }
  
  window_data <- returns_matrix_full[1:end_period, ]
  
  # ✅ Calculer les rendements CAPM pour cette fenêtre
  betas_window <- apply(window_data, 2, function(x) {
    cov(x, spy_returns_full[1:end_period], use = "complete.obs") /
      var(spy_returns_full[1:end_period], na.rm = TRUE)
  })
  
  mean_spy_window <- geometric_mean_returns(spy_returns_full[1:end_period])
  mean_spy_annual_window <- ((1 + mean_spy_window) ^ 12 - 1) * 100
  mean_rf_window <- geometric_mean_rf(rf_data$rf_rate[1:end_period])
  
  # ✅ Rendements CAPM pour cette fenêtre
  expected_returns_window <- mean_rf_window + 
    betas_window * 
    (mean_spy_annual_window - mean_rf_window)
  
  # Calculer les portefeuilles
  gmvp <- calculate_gmvp_weights(window_data)
  
  # ✅ Passer les rendements CAPM à MV !
  mv <- calculate_mv_portfolio(window_data, expected_returns_window)
  
  optimization_results[[year]] <- list(
    Year = year,
    Periods = end_period,
    GMVP_weights = gmvp$weights,
    GMVP_vol = gmvp$volatility,
    MV_weights = mv$weights,          # ← Différent de GMVP !
    MV_vol = mv$volatility,
    MV_target_return = mv$target
  )
}


cat("✓ Optimization complete (", length(optimization_results), " annual rebalances)\n\n")

# ============================================================
# SECTION 7: OUT-OF-SAMPLE BACKTESTING & PERFORMANCE
# ============================================================

cat("--- OUT-OF-SAMPLE BACKTESTING ---\n")

# Initialize wealth tracking
wealth_data <- data.frame(
  Year_Index = 0,
  Wealth_GMVP_Net = 1,
  Wealth_GMVP_Gross = 1,
  Wealth_MV_Net = 1,
  Wealth_MV_Gross = 1,
  Wealth_SP500 = 1,
  Wealth_GMVP_Net_Real = 1,
  Wealth_GMVP_Gross_Real = 1,
  Wealth_MV_Net_Real = 1,
  Wealth_MV_Gross_Real = 1,
  Wealth_SP500_Real = 1,
  stringsAsFactors = FALSE
)

# Cumulative wealth variables
wealth_gmvp_net <- 1
wealth_gmvp_gross <- 1
wealth_mv_net <- 1
wealth_mv_gross <- 1
wealth_sp500 <- 1

wealth_gmvp_net_real <- 1
wealth_gmvp_gross_real <- 1
wealth_mv_net_real <- 1
wealth_mv_gross_real <- 1
wealth_sp500_real <- 1

performance_summary <- list()

# Annual backtesting loop
for (year in 1:(length(optimization_results))) {
  
  initial_period <- initial_periods + year * periods_per_year
  end_period <- initial_periods + ((year + 1) * periods_per_year) - 1
  
  if (end_period > nrow(returns_matrix_full)) {
    break
  }
  # Extract monthly returns for the rebalancing year
  monthly_returns <- returns_matrix_full[initial_period:end_period, ]
  spy_returns_year <- spy_returns_full[initial_period:end_period]
  
  # --- TURNOVER & TRANSACTION COSTS ---
  if (year == 1) {
    # First year: assume equal-weight starting portfolio
    equal_weight <- rep(
      1 / length(optimization_results[[year]]$GMVP_weights),
      length(optimization_results[[year]]$GMVP_weights)
    )
    turnover_gmvp <- sum(abs(optimization_results[[year]]$GMVP_weights - equal_weight))
    turnover_mv <- sum(abs(optimization_results[[year]]$MV_weights - equal_weight))
  } else {
    # Subsequent years: compare to previous portfolio
    turnover_gmvp <- sum(abs(optimization_results[[year]]$GMVP_weights -
                               optimization_results[[year - 1]]$GMVP_weights))
    turnover_mv <- sum(abs(optimization_results[[year]]$MV_weights -
                             optimization_results[[year - 1]]$MV_weights))
  }
  
  transaction_costs_gmvp <- turnover_gmvp * transaction_cost_rate * 100
  transaction_costs_mv <- turnover_mv * transaction_cost_rate * 100
  
  # --- PORTFOLIO RETURNS ---
  
  portfolio_return_gmvp <- rowSums(sweep(monthly_returns, 2,
                                         optimization_results[[year]]$GMVP_weights, "*"))
  portfolio_return_mv <- rowSums(sweep(monthly_returns, 2,
                                       optimization_results[[year]]$MV_weights, "*"))
  
  # Annual compound returns (gross of costs)
  annual_return_gmvp_gross <- (prod(1 + portfolio_return_gmvp / 100, na.rm = TRUE) - 1) * 100
  annual_return_mv_gross <- (prod(1 + portfolio_return_mv / 100, na.rm = TRUE) - 1) * 100
  annual_return_sp500_gross <- (prod(1 + spy_returns_year / 100, na.rm = TRUE) - 1) * 100
  
  # Net of transaction costs
  annual_return_gmvp <- annual_return_gmvp_gross - transaction_costs_gmvp
  annual_return_mv <- annual_return_mv_gross - transaction_costs_mv
  
  # --- WEALTH ACCUMULATION ---
  
  wealth_gmvp_net <- wealth_gmvp_net * (1 + annual_return_gmvp / 100)
  wealth_gmvp_gross <- wealth_gmvp_gross * (1 + annual_return_gmvp_gross / 100)
  wealth_mv_net <- wealth_mv_net * (1 + annual_return_mv / 100)
  wealth_mv_gross <- wealth_mv_gross * (1 + annual_return_mv_gross / 100)
  wealth_sp500 <- wealth_sp500 * (1 + annual_return_sp500_gross / 100)
  
  # --- INFLATION ADJUSTMENT ---
  
  if (year <= 20) {
    # CPI data available through Nov 2024 (year 20)
    cpi_ratio <- cpi_data$cpi[year * periods_per_year] / cpi_data$cpi[1]
    
    wealth_gmvp_net_real <- wealth_gmvp_net / cpi_ratio
    wealth_gmvp_gross_real <- wealth_gmvp_gross / cpi_ratio
    wealth_mv_net_real <- wealth_mv_net / cpi_ratio
    wealth_mv_gross_real <- wealth_mv_gross / cpi_ratio
    wealth_sp500_real <- wealth_sp500 / cpi_ratio
  } else {
    # CPI data unavailable after Nov 2024 (year 21+)
    wealth_gmvp_net_real <- NA
    wealth_gmvp_gross_real <- NA
    wealth_mv_net_real <- NA
    wealth_mv_gross_real <- NA
    wealth_sp500_real <- NA
  }
  
  # --- RISK METRICS ---
  
  vol_gmvp <- sd(portfolio_return_gmvp, na.rm = TRUE) * sqrt(12)
  vol_mv <- sd(portfolio_return_mv, na.rm = TRUE) * sqrt(12)
  vol_sp500 <- sd(spy_returns_year, na.rm = TRUE) * sqrt(12)
  
  rf_subset <- rf_data$rf_rate[initial_period:end_period]
  rf_year <- geometric_mean_rf(rf_data$rf_rate[initial_period:end_period])
  
  sharpe_gmvp <- (annual_return_gmvp - rf_year) / vol_gmvp
  sharpe_mv <- (annual_return_mv - rf_year) / vol_mv
  sharpe_sp500 <- (annual_return_sp500_gross - rf_year) / vol_sp500
  
  # --- STORE RESULTS ---
  
  performance_summary[[year]] <- list(
    Year_Index = year,
    Year_Calendar = 2004 + year,
    Annual_Return_GMVP_Pct = round(annual_return_gmvp, 2),
    Annual_Return_MV_Pct = round(annual_return_mv, 2),
    Annual_Return_SP500_Pct = round(annual_return_sp500_gross, 2),
    Volatility_GMVP_Pct = round(vol_gmvp, 2),
    Volatility_MV_Pct = round(vol_mv, 2),
    Volatility_SP500_Pct = round(vol_sp500, 2),
    Sharpe_GMVP = round(sharpe_gmvp, 3),
    Sharpe_MV = round(sharpe_mv, 3),
    Sharpe_SP500 = round(sharpe_sp500, 3),
    Turnover_GMVP = round(turnover_gmvp, 3),
    Turnover_MV = round(turnover_mv, 3),
    Trans_Cost_GMVP_Bps = round(transaction_costs_gmvp, 4),
    Trans_Cost_MV_Bps = round(transaction_costs_mv, 4)
  )
  
  wealth_data <- rbind(wealth_data, data.frame(
    Year_Index = year,
    Wealth_GMVP_Net = wealth_gmvp_net,
    Wealth_GMVP_Gross = wealth_gmvp_gross,
    Wealth_MV_Net = wealth_mv_net,
    Wealth_MV_Gross = wealth_mv_gross,
    Wealth_SP500 = wealth_sp500,
    Wealth_GMVP_Net_Real = wealth_gmvp_net_real,
    Wealth_GMVP_Gross_Real = wealth_gmvp_gross_real,
    Wealth_MV_Net_Real = wealth_mv_net_real,
    Wealth_MV_Gross_Real = wealth_mv_gross_real,
    Wealth_SP500_Real = wealth_sp500_real,
    stringsAsFactors = FALSE
  ))
}

performance_df <- do.call(rbind, lapply(performance_summary, as.data.frame))
rownames(performance_df) <- NULL

cat("✓ Backtesting complete\n")
cat("Sample annual performance (first 5 years):\n")
print(head(performance_df, 5))
cat("\n")

# ============================================================
# SECTION 8: VISUALIZATION
# ============================================================


cat("--- GENERATING VISUALIZATIONS ---\n")

start_date <- ymd("2004-11-01")

wealth_long <- wealth_data %>%
  mutate(Date = start_date + years(Year_Index)) %>%
  pivot_longer(
    cols = contains("Wealth"),
    names_to = "Portfolio",
    values_to = "Wealth"
  ) %>%
  mutate(
    Type = case_when(
      grepl("_Net", Portfolio) ~ "Net",
      grepl("_Gross", Portfolio) ~ "Gross",
      grepl("SP500", Portfolio) ~ "Benchmark",
      TRUE ~ NA_character_
    ),
    Portfolio_Type = case_when(
      grepl("GMVP", Portfolio) ~ "GMVP",
      grepl("MV", Portfolio) ~ "MV",
      grepl("SP500", Portfolio) ~ "S&P 500",
      TRUE ~ NA_character_
    ),
    Real_Nominal = case_when(
      grepl("_Real", Portfolio) ~ "Real (CPI-Adjusted)",
      TRUE ~ "Nominal"
    )
  ) %>%
  filter(!is.na(Type)) %>%
  arrange(Date)

# Plot 1: Nominal Wealth Evolution
p1 <- ggplot(wealth_long %>% filter(Real_Nominal == "Nominal"),
             aes(x = Date, y = Wealth,
                 color = paste(Portfolio_Type, Type, sep = "_"),
                 linetype = Type)) +
  geom_line(size = 1.1) +
  scale_linetype_manual(
    values = c("Gross" = "solid", "Net" = "dashed", "Benchmark" = "dotted"),
    name = "Strategy Type"
  ) +
  scale_color_manual(
    values = c(
      "GMVP_Gross" = "#E41A1C", "GMVP_Net" = "#FF6B6B",
      "MV_Gross" = "#377EB8", "MV_Net" = "#6FB3E8",
      "S&P 500_Benchmark" = "#2CA02C"
    ),
    name = "Strategy",
    labels = c(
      "GMVP_Gross" = "GMVP Gross",
      "GMVP_Net" = "GMVP Net",
      "MV_Gross" = "MV Gross",
      "MV_Net" = "MV Net",
      "S&P 500_Benchmark" = "S&P 500"
    )
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
  labs(
    title = "Wealth Evolution: Optimized Portfolios vs S&P 500 Benchmark",
    subtitle = "Nominal Values | Nov 2004 - Nov 2025 | Initial Investment = $1",
    x = "Year",
    y = "Cumulative Wealth",
    caption = "GMVP = Global Minimum Variance | MV = Mean-Variance Efficient | Net = After 10bps transaction costs"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray90", size = 0.3)
  )

print(p1)

# Plot 2: Real Wealth Evolution (CPI-Adjusted)
p2 <- ggplot(wealth_long %>% filter(Real_Nominal == "Real (CPI-Adjusted)"),
             aes(x = Date, y = Wealth,
                 color = paste(Portfolio_Type, Type, sep = "_"),
                 linetype = Type)) +
  geom_line(size = 1.1) +
  scale_linetype_manual(
    values = c("Gross" = "solid", "Net" = "dashed", "Benchmark" = "dotted"),
    name = "Strategy Type"
  ) +
  scale_color_manual(
    values = c(
      "GMVP_Gross" = "#E41A1C", "GMVP_Net" = "#FF6B6B",
      "MV_Gross" = "#377EB8", "MV_Net" = "#6FB3E8",
      "S&P 500_Benchmark" = "#2CA02C"
    ),
    name = "Strategy",
    labels = c(
      "GMVP_Gross" = "GMVP Gross",
      "GMVP_Net" = "GMVP Net",
      "MV_Gross" = "MV Gross",
      "MV_Net" = "MV Net",
      "S&P 500_Benchmark" = "S&P 500"
    )
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
  labs(
    title = "Real Wealth Evolution: Optimized Portfolios vs S&P 500 Benchmark",
    subtitle = "CPI-Adjusted to Nov 2004 Dollars | Nov 2004 - Nov 2024",
    x = "Year",
    y = "Cumulative Wealth (Nov 2004 $)",
    caption = "Data: Nov 2004 - Nov 2024 | Deflated by Bureau of Labor Statistics CPI-U"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray90", size = 0.3)
  )

print(p2)

cat("✓ Visualizations saved (PNG format, 300 DPI)\n\n")

# ============================================================
# SECTION 9: COMPREHENSIVE SUMMARY STATISTICS
# ============================================================

cat("--- FINAL PERFORMANCE SUMMARY ---\n\n")

final_summary <- data.frame(
  Strategy = c("GMVP (Gross)", "GMVP (Net)", "MV (Gross)", "MV (Net)", "S&P 500"),
  Wealth_2025_Nominal = c(
    tail(wealth_data$Wealth_GMVP_Gross, 1),
    tail(wealth_data$Wealth_GMVP_Net, 1),
    tail(wealth_data$Wealth_MV_Gross, 1),
    tail(wealth_data$Wealth_MV_Net, 1),
    tail(wealth_data$Wealth_SP500, 1)
  ),
  Wealth_2024_Real_Nov2004Dollars = c(
    wealth_data$Wealth_GMVP_Gross_Real[21],
    wealth_data$Wealth_GMVP_Net_Real[21],
    wealth_data$Wealth_MV_Gross_Real[21],
    wealth_data$Wealth_MV_Net_Real[21],
    wealth_data$Wealth_SP500_Real[21]
  ),
  Annualized_Return_Nominal_Pct = NA,
  Annualized_Volatility_Pct = NA,
  Sharpe_Ratio = NA
)

# Annualized return (21-year period: Nov 2004 - Nov 2025)
final_summary$Annualized_Return_Nominal_Pct <-
  (final_summary$Wealth_2025_Nominal ^ (1 / 21) - 1) * 100

# Average metrics from annual performance
final_summary$Annualized_Volatility_Pct <- c(
  mean(performance_df$Volatility_GMVP_Pct),
  mean(performance_df$Volatility_GMVP_Pct),
  mean(performance_df$Volatility_MV_Pct),
  mean(performance_df$Volatility_MV_Pct),
  mean(performance_df$Volatility_SP500_Pct)
)

final_summary$Sharpe_Ratio <- c(
  mean(performance_df$Sharpe_GMVP),
  mean(performance_df$Sharpe_GMVP),
  mean(performance_df$Sharpe_MV),
  mean(performance_df$Sharpe_MV),
  mean(performance_df$Sharpe_SP500)
)

final_summary <- final_summary %>%
  mutate(
    Wealth_2025_Nominal = round(Wealth_2025_Nominal, 2),
    Wealth_2024_Real_Nov2004Dollars = round(Wealth_2024_Real_Nov2004Dollars, 2),
    Annualized_Return_Nominal_Pct = round(Annualized_Return_Nominal_Pct, 2),
    Annualized_Volatility_Pct = round(Annualized_Volatility_Pct, 2),
    Sharpe_Ratio = round(Sharpe_Ratio, 3)
  )

print(final_summary)
cat("\nInterpretation Guide:\n")
cat("  Wealth columns: Dollar value of $1 initial investment\n")
cat("  Real values: Deflated to Nov 2004 purchasing power\n")
cat("  Annualized metrics: Geometric mean returns & average volatility\n")
cat("  Sharpe Ratio: Risk-adjusted return (higher is better)\n\n")

# ============================================================
# SECTION 10: DATA EXPORT
# ============================================================
# Define output prefix
output_prefix <- "portfolio_2025_leaders_"

cat("--- EXPORTING RESULTS ---\n")

write.csv(capm_summary, paste0(output_prefix, "capm_summary.csv"), row.names = FALSE)
cat("✓ CAPM Analysis →", paste0(output_prefix, "capm_summary.csv\n"))

write.csv(performance_df, paste0(output_prefix, "annual_performance_detail.csv"), row.names = FALSE)
cat("✓ Annual Performance →", paste0(output_prefix, "annual_performance_detail.csv\n"))

write.csv(final_summary, paste0(output_prefix, "final_summary_statistics.csv"), row.names = FALSE)
cat("✓ Summary Statistics →", paste0(output_prefix, "final_summary_statistics.csv\n"))

write.csv(wealth_data, paste0(output_prefix, "wealth_evolution_timeseries.csv"), row.names = FALSE)
cat("✓ Wealth Time Series →", paste0(output_prefix, "wealth_evolution_timeseries.csv\n"))

gfsave(paste0(output_prefix, "wealth_nominal.png"), p1, width = 13, height = 7, dpi = 300)
ggsave(paste0(output_prefix, "wealth_real_cpi_adjusted.png"), p2, width = 13, height = 7, dpi = 300)

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("✓ ANALYSIS COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Output files:\n")
cat("  - Visualizations:", paste0(output_prefix, "wealth_nominal.png,"), paste0(output_prefix, "wealth_real_cpi_adjusted.png\n"))
cat("  - Data: 4 CSV files\n")
cat("Period: Nov 2004 - Nov 2025 (21 years)\n")
cat("Nominal results: Full period through Nov 2025\n")
cat("Real results: Through Nov 2024 (CPI data limitation)\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
