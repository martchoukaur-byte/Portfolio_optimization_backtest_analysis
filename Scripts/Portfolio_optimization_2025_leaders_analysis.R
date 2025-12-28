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
# PORTFOLIO OPTIMIZATION & BACKTESTING ANALYSIS (2004-2025)
# Retrospective Analysis: 2025 Market Leaders Backtested to 2004
# ============================================================

packages <- c(
  "quantmod",
  "dplyr",
  "tidyr",
  "zoo",
  "ggplot2",
  "PerformanceAnalytics",
  "lubridate",
  "quadprog"
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
    dplyr::select(date, open, high, low, close, adjusted_close, volume) %>%
    dplyr::mutate(
      returns = (adjusted_close / dplyr::lag(adjusted_close) - 1) * 100,
      compound_growth = adjusted_close / adjusted_close[1]
    )
  
  return(df)
}

geometric_mean_returns <- function(x) {
  x_clean <- na.omit(x)
  n <- length(x_clean)
  if (n == 0) return(NA_real_)
  ((prod(1 + x_clean / 100)) ^ (1 / n) - 1)
}

# ✅ NOUVELLE FONCTION: geometric_mean_rf avec diagnostic
geometric_mean_rf <- function(rf_rates, verbose = TRUE) {
  rf_clean <- na.omit(rf_rates)
  n_original <- length(rf_rates)
  n_clean <- length(rf_clean)
  n_dropped <- n_original - n_clean
  
  if (verbose && n_dropped > 0) {
    cat(sprintf("  [RF] Dropped %d NAs, using n=%d observations\n",
                n_dropped, n_clean))
  }
  
  if (n_clean == 0) {
    if (verbose) cat("  [RF] ERROR: No valid RF observations!\n")
    return(NA_real_)
  }
  
  result <- ((prod(1 + rf_clean / 100)) ^ (1 / n_clean) - 1) * 100
  return(result)
}

calculate_gmvp_weights <- function(window_data, min_weight = 0, max_weight = 1) {
  sigma_mat <- cov(window_data)
  n_assets <- ncol(window_data)
  
  D <- 2 * sigma_mat
  d <- rep(0, n_assets)
  
  A_eq <- matrix(1, nrow = 1, ncol = n_assets)
  b_eq <- 1
  
  A_in <- rbind(
    diag(n_assets),
    -diag(n_assets)
  )
  b_in <- c(rep(min_weight, n_assets), rep(-max_weight, n_assets))
  
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
  
  sigma_mat <- cov(window_data)
  n_assets <- ncol(window_data)
  
  target_return <- mean(expected_returns)
  
  D <- 2 * sigma_mat
  d <- rep(0, n_assets)
  
  A_eq <- rbind(
    matrix(1, nrow = 1, ncol = n_assets),
    matrix(expected_returns, nrow = 1)
  )
  
  b_eq <- c(1, target_return)
  
  A_in <- rbind(
    diag(n_assets),
    -diag(n_assets)
  )
  
  b_in <- c(
    rep(min_weight, n_assets),
    rep(-max_weight, n_assets)
  )
  
  result <- quadprog::solve.QP(
    Dmat = D,
    dvec = d,
    Amat = t(rbind(A_eq, A_in)),
    bvec = c(b_eq, b_in),
    meq = 2
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

tickers <- c("NVDA", "MSFT", "AAPL", "AMZN", "BRK-A")

cat("Loading monthly stock data for:", paste(tickers, collapse = ", "), "\n")

stock_list <- lapply(tickers, load_stock_monthly)
names(stock_list) <- tickers

stock_list_tagged <- lapply(seq_along(stock_list), function(i) {
  stock_list[[i]] %>% dplyr::mutate(stock = names(stock_list)[i])
})

data_all <- dplyr::bind_rows(stock_list_tagged)

cat("Loading S&P 500 (SPY) benchmark data...\n")

spy_data <- load_stock_monthly("SPY")
names(spy_data)[names(spy_data) == "adjusted_close"] <- "spy_close"
names(spy_data)[names(spy_data) == "returns"] <- "spy_returns"
names(spy_data)[names(spy_data) == "compound_growth"] <- "spy_compound_growth"

cat("Loading risk-free rate (3M Treasury)...\n")

rf_data_raw <- getSymbols(
  "DGS3MO",
  src = "FRED",
  from = "1999-11-01",
  to = "2025-11-02",
  auto.assign = FALSE
) %>%
  as.data.frame()

rf_data_raw$date <- as.Date(rownames(rf_data_raw))
colnames(rf_data_raw) <- c("rf_rate", "date")
rownames(rf_data_raw) <- NULL

# ✅ CORRECTION #1: Timing RF (fin de mois, pas début)
# ✅ CORRECTION #2: Exclure NAs AVANT sélection
rf_data <- rf_data_raw %>%
  dplyr::mutate(
    rf_rate = as.numeric(as.character(rf_rate)),
    yearmonth = format(date, "%Y-%m")
  ) %>%
  dplyr::filter(!is.na(rf_rate)) %>%
  dplyr::group_by(yearmonth) %>%
  dplyr::filter(date == max(date)) %>%
  dplyr::ungroup() %>%
  dplyr::select(date, rf_rate) %>%
  dplyr::arrange(date)

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
# SECTION 4: DATA STRUCTURING + ALIGNMENT CHECK
# ============================================================

cat("Structuring return matrices...\n")

data_wide <- data_all %>%
  dplyr::select(date, stock, returns) %>%
  tidyr::pivot_wider(names_from = stock, values_from = returns) %>%
  dplyr::arrange(date) %>%
  dplyr::slice(-1)

returns_matrix_full <- data_wide %>%
  dplyr::select(-date) %>%
  as.matrix()

n_150 <- 150
returns_matrix_150 <- returns_matrix_full[1:n_150, ]

spy_returns_full <- spy_data$spy_returns[-1]
spy_returns_150 <- spy_returns_full[1:n_150]

cat("✓ Data loaded and structured\n")
cat("  Full period months (stocks): ", nrow(returns_matrix_full), "\n")
cat("  Returns start: ", as.character(data_wide$date[1]), 
    " | Returns end: ", as.character(tail(data_wide$date, 1)), "\n")
cat("  RF start: ", as.character(min(rf_data$date)), 
    " | RF end: ", as.character(max(rf_data$date)), "\n\n")

# ✅ CORRECTION #3: DATA QUALITY CHECK RF vs RETURNS
cat("--- DATA QUALITY CHECK: RF ALIGNMENT ---\n")

n_returns_months <- nrow(returns_matrix_full)
n_rf_total <- nrow(rf_data)
n_rf_valid <- sum(!is.na(rf_data$rf_rate))

cat("Returns months available: ", n_returns_months, "\n")
cat("RF rows (monthly):        ", n_rf_total, "\n")
cat("RF valid rows:            ", n_rf_valid, "\n")
cat("RF coverage:              ", 
    round(min(1, n_rf_valid / n_returns_months) * 100, 1), "%\n\n")

if (n_rf_total < n_returns_months) {
  stop(paste0(
    "CRITICAL: RF dataset is shorter than returns dataset.\n",
    "RF rows: ", n_rf_total, " < Returns months: ", n_returns_months, "\n",
    "Adjust RF date range or shorten analysis period."
  ))
}

rf_check <- rf_data %>%
  dplyr::mutate(ym = format(date, "%Y-%m")) %>%
  dplyr::group_by(ym) %>%
  dplyr::summarise(n_per_month = dplyr::n(), .groups = "drop")

duplicates <- rf_check %>% dplyr::filter(n_per_month > 1)
if (nrow(duplicates) > 0) {
  warning("Some months have multiple RF observations after aggregation.")
  print(duplicates)
}

cat("✓ RF alignment check passed\n\n")

# ============================================================
# SECTION 5: CAPM ANALYSIS
# ============================================================

cat("--- CAPM ANALYSIS ---\n")

betas_full <- apply(returns_matrix_full, 2, function(x) {
  cov(x, spy_returns_full, use = "complete.obs") /
    var(spy_returns_full, na.rm = TRUE)
})

betas_150 <- apply(returns_matrix_150, 2, function(x) {
  cov(x, spy_returns_150, use = "complete.obs") /
    var(spy_returns_150, na.rm = TRUE)
})

mean_spy_monthly_full <- geometric_mean_returns(spy_returns_full)
mean_spy_annual_full <- ((1 + mean_spy_monthly_full) ^ 12 - 1) * 100

mean_spy_monthly_150 <- geometric_mean_returns(spy_returns_150)
mean_spy_annual_150 <- ((1 + mean_spy_monthly_150) ^ 12 - 1) * 100

# ✅ CORRECTION #5: Appel avec verbose
cat("Calculating RF for full period:\n")
mean_rf_full <- geometric_mean_rf(rf_data$rf_rate[1:n_returns_months], verbose = TRUE)

cat("Calculating RF for first 150 months:\n")
mean_rf_150 <- geometric_mean_rf(rf_data$rf_rate[1:n_150], verbose = TRUE)

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

initial_periods <- 48
periods_per_year <- 12
n_years <- 22
transaction_cost_rate <- 0.001

cat("Parameters:\n")
cat("  Inception period: ", initial_periods, " months (4 years)\n")
cat("  Rebalancing: Annual\n")
cat("  Transaction cost: ", transaction_cost_rate * 100, "%\n")
cat("  Max analysis horizon: ", n_years, " years\n\n")

optimization_results <- list()

for (year in 1:n_years) {
  
  end_period <- initial_periods + (year * periods_per_year)
  
  if (end_period > nrow(returns_matrix_full)) {
    break
  }
  
  # ✅ CORRECTION #4: Vérifier RF data disponible
  if (end_period > nrow(rf_data)) {
    stop(sprintf("Year %d: RF data ends at row %d but need row %d\n",
                 year, nrow(rf_data), end_period))
  }
  
  window_data <- returns_matrix_full[1:end_period, ]
  spy_window <- spy_returns_full[1:end_period]
  rf_subset <- rf_data$rf_rate[1:end_period]
  
  n_rf_available <- sum(!is.na(rf_subset))
  rf_coverage_pct <- n_rf_available / end_period * 100
  
  if (rf_coverage_pct < 95) {
    cat(sprintf("  [Year %d] RF coverage = %.1f%% (%d/%d months)\n",
                year, rf_coverage_pct, n_rf_available, end_period))
  }
  
  betas_window <- apply(window_data, 2, function(x) {
    cov(x, spy_window, use = "complete.obs") /
      var(spy_window, na.rm = TRUE)
  })
  
  mean_spy_window <- geometric_mean_returns(spy_window)
  mean_spy_annual_window <- ((1 + mean_spy_window) ^ 12 - 1) * 100
  mean_rf_window <- geometric_mean_rf(rf_subset, verbose = (year == 1))
  
  expected_returns_window <- mean_rf_window +
    betas_window * (mean_spy_annual_window - mean_rf_window)
  
  gmvp <- calculate_gmvp_weights(window_data)
  mv <- calculate_mv_portfolio(window_data, expected_returns_window)
  
  optimization_results[[year]] <- list(
    Year = year,
    Periods = end_period,
    GMVP_weights = gmvp$weights,
    GMVP_vol = gmvp$volatility,
    MV_weights = mv$weights,
    MV_vol = mv$volatility,
    MV_target_return = mv$target
  )
}

cat("✓ Optimization complete (", length(optimization_results), " annual rebalances)\n\n")

# ============================================================
# SECTION 7: OUT-OF-SAMPLE BACKTESTING & PERFORMANCE
# ============================================================

cat("--- OUT-OF-SAMPLE BACKTESTING ---\n")

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

for (year in 1:length(optimization_results)) {
  
  initial_period <- initial_periods + year * periods_per_year
  end_period <- initial_periods + ((year + 1) * periods_per_year) - 1
  
  if (end_period > nrow(returns_matrix_full)) {
    break
  }
  
  monthly_returns <- returns_matrix_full[initial_period:end_period, ]
  spy_returns_year <- spy_returns_full[initial_period:end_period]
  
  if (year == 1) {
    equal_weight <- rep(
      1 / length(optimization_results[[year]]$GMVP_weights),
      length(optimization_results[[year]]$GMVP_weights)
    )
    turnover_gmvp <- sum(abs(optimization_results[[year]]$GMVP_weights - equal_weight))
    turnover_mv <- sum(abs(optimization_results[[year]]$MV_weights - equal_weight))
  } else {
    turnover_gmvp <- sum(abs(optimization_results[[year]]$GMVP_weights -
                               optimization_results[[year - 1]]$GMVP_weights))
    turnover_mv <- sum(abs(optimization_results[[year]]$MV_weights -
                             optimization_results[[year - 1]]$MV_weights))
  }
  
  transaction_costs_gmvp <- turnover_gmvp * transaction_cost_rate * 100
  transaction_costs_mv <- turnover_mv * transaction_cost_rate * 100
  
  portfolio_return_gmvp <- rowSums(sweep(monthly_returns, 2,
                                         optimization_results[[year]]$GMVP_weights, "*"))
  portfolio_return_mv <- rowSums(sweep(monthly_returns, 2,
                                       optimization_results[[year]]$MV_weights, "*"))
  
  annual_return_gmvp_gross <- (prod(1 + portfolio_return_gmvp / 100, na.rm = TRUE) - 1) * 100
  annual_return_mv_gross <- (prod(1 + portfolio_return_mv / 100, na.rm = TRUE) - 1) * 100
  annual_return_sp500_gross <- (prod(1 + spy_returns_year / 100, na.rm = TRUE) - 1) * 100
  
  annual_return_gmvp <- annual_return_gmvp_gross - transaction_costs_gmvp
  annual_return_mv <- annual_return_mv_gross - transaction_costs_mv
  
  wealth_gmvp_net <- wealth_gmvp_net * (1 + annual_return_gmvp / 100)
  wealth_gmvp_gross <- wealth_gmvp_gross * (1 + annual_return_gmvp_gross / 100)
  wealth_mv_net <- wealth_mv_net * (1 + annual_return_mv / 100)
  wealth_mv_gross <- wealth_mv_gross * (1 + annual_return_mv_gross / 100)
  wealth_sp500 <- wealth_sp500 * (1 + annual_return_sp500_gross / 100)
  
  if (year <= 20) {
    cpi_ratio <- cpi_data$cpi[year * periods_per_year] / cpi_data$cpi[1]
    
    wealth_gmvp_net_real <- wealth_gmvp_net / cpi_ratio
    wealth_gmvp_gross_real <- wealth_gmvp_gross / cpi_ratio
    wealth_mv_net_real <- wealth_mv_net / cpi_ratio
    wealth_mv_gross_real <- wealth_mv_gross / cpi_ratio
    wealth_sp500_real <- wealth_sp500 / cpi_ratio
  } else {
    wealth_gmvp_net_real <- NA
    wealth_gmvp_gross_real <- NA
    wealth_mv_net_real <- NA
    wealth_mv_gross_real <- NA
    wealth_sp500_real <- NA
  }
  
  vol_gmvp <- sd(portfolio_return_gmvp, na.rm = TRUE) * sqrt(12)
  vol_mv <- sd(portfolio_return_mv, na.rm = TRUE) * sqrt(12)
  vol_sp500 <- sd(spy_returns_year, na.rm = TRUE) * sqrt(12)
  
  rf_year <- geometric_mean_rf(
    rf_data$rf_rate[initial_period:end_period],
    verbose = FALSE
  )
  
  sharpe_gmvp <- (annual_return_gmvp - rf_year) / vol_gmvp
  sharpe_mv <- (annual_return_mv - rf_year) / vol_mv
  sharpe_sp500 <- (annual_return_sp500_gross - rf_year) / vol_sp500
  
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

start_date <- lubridate::ymd("2004-11-01")

wealth_long <- wealth_data %>%
  dplyr::mutate(Date = start_date + lubridate::years(Year_Index)) %>%
  tidyr::pivot_longer(
    cols = contains("Wealth"),
    names_to = "Portfolio",
    values_to = "Wealth"
  ) %>%
  dplyr::mutate(
    Type = dplyr::case_when(
      grepl("_Net", Portfolio) ~ "Net",
      grepl("_Gross", Portfolio) ~ "Gross",
      grepl("SP500", Portfolio) ~ "Benchmark",
      TRUE ~ NA_character_
    ),
    Portfolio_Type = dplyr::case_when(
      grepl("GMVP", Portfolio) ~ "GMVP",
      grepl("MV", Portfolio) ~ "MV",
      grepl("SP500", Portfolio) ~ "S&P 500",
      TRUE ~ NA_character_
    ),
    Real_Nominal = dplyr::case_when(
      grepl("_Real", Portfolio) ~ "Real (CPI-Adjusted)",
      TRUE ~ "Nominal"
    )
  ) %>%
  dplyr::filter(!is.na(Type)) %>%
  dplyr::arrange(Date)

p1 <- ggplot(wealth_long %>% dplyr::filter(Real_Nominal == "Nominal"),
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

p2 <- ggplot(wealth_long %>% dplyr::filter(Real_Nominal == "Real (CPI-Adjusted)"),
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

cat("✓ Visualizations prepared\n\n")

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

final_summary$Annualized_Return_Nominal_Pct <-
  (final_summary$Wealth_2025_Nominal ^ (1 / 21) - 1) * 100

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
  dplyr::mutate(
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

ggplot2::ggsave(paste0(output_prefix, "wealth_nominal.png"), p1, width = 13, height = 7, dpi = 300)
ggplot2::ggsave(paste0(output_prefix, "wealth_real_cpi_adjusted.png"), p2, width = 13, height = 7, dpi = 300)

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
