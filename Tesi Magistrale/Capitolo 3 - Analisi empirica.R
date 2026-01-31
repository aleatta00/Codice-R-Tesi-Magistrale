###############################################################################
# CODICE PER CAPITOLO 3: ANALISI EMPIRICA - DAX & S&P500
#
# Il presente codice implementa l'analisi empirica discussa nel Capitolo 3,
# applicando i metodi di inferenza bayesiana sviluppati nel Capitolo 2
# a dati reali di mercati finanziari.
#
# Contenuto:
#   - Caricamento e statistiche descrittive dei dati
#   - Trasformazione in rendimenti logaritmici
#   - Analisi esplorativa (grafici, correlazioni)
#   - Trasformazione a pseudo-osservazioni
#   - Inferenza bayesiana per copule Clayton, Gumbel, Frank
#   - Model comparison (DIC, WAIC)
#   - Diagnostica MCMC
###############################################################################

# -------------------------------------------------------------
# Librerie
# -------------------------------------------------------------
suppressPackageStartupMessages({
  library(copula)
  library(ggplot2)
  library(gridExtra)
  library(moments)
  library(coda)
  library(readxl)
  library(dplyr)
})

set.seed(42)

###############################################################################
# SEZIONE 3.1: CARICAMENTO E DESCRIZIONE DATI
###############################################################################

# -------------------------------------------------------------
# Caricamento dataset
# -------------------------------------------------------------
data <- read_excel("Dataset.xlsx")
data <- na.omit(data)
n <- nrow(data)

cat("=== Dataset Caricato ===\n")
cat("Numero osservazioni:", n, "\n")
cat("Periodo:", as.character(min(data$Date)), "-", 
    as.character(max(data$Date)), "\n\n")

# -------------------------------------------------------------
# Statistiche descrittive prezzi
# -------------------------------------------------------------
cat("=== Statistiche Descrittive Prezzi ===\n\n")

summary_prices <- data.frame(
  Statistica = c("Media", "Mediana", "Dev.Std", "Min", "Max", 
                 "Asimmetria", "Curtosi"),
  DAX = c(
    mean(data$DAX),
    median(data$DAX),
    sd(data$DAX),
    min(data$DAX),
    max(data$DAX),
    skewness(data$DAX),
    kurtosis(data$DAX)
  ),
  SP500 = c(
    mean(data$SP500),
    median(data$SP500),
    sd(data$SP500),
    min(data$SP500),
    max(data$SP500),
    skewness(data$SP500),
    kurtosis(data$SP500)
  )
)
summary_prices[, c("DAX", "SP500")] <-
  round(summary_prices[, c("DAX", "SP500")], 2)
print(summary_prices, 2)

# -------------------------------------------------------------
# Visualizzazione serie temporali prezzi
# -------------------------------------------------------------

plot(data$Date, data$DAX, type = "l", col = "blue", lwd = 1.5,
     xlab = "Data", ylab = "DAX", 
     main = "Serie Temporale DAX (2005-2024)")
grid()

plot(data$Date, data$SP500, type = "l", col = "red", lwd = 1.5,
     xlab = "Data", ylab = "S&P 500", 
     main = "Serie Temporale S&P 500 (2005-2024)")
grid()


###############################################################################
# SEZIONE 3.2: TRASFORMAZIONE IN RENDIMENTI LOGARITMICI
###############################################################################

# -------------------------------------------------------------
# Calcolo rendimenti logaritmici
# -------------------------------------------------------------
data$r_DAX <- c(NA, diff(log(data$DAX)))
data$r_SP500 <- c(NA, diff(log(data$SP500)))
data <- na.omit(data)

# Conversione in percentuale per tabelle
returns_pct <- data.frame(
  r_DAX = data$r_DAX * 100,
  r_SP500 = data$r_SP500 * 100
)

# -------------------------------------------------------------
# Statistiche descrittive rendimenti
# -------------------------------------------------------------
cat("\n=== Statistiche Descrittive Rendimenti (%) ===\n\n")

summary_returns <- data.frame(
  Statistica = c("Media", "Mediana", "Dev.Std", "Min", "Max", 
                 "Asimmetria", "Curtosi (eccesso)"),
  DAX = c(
    mean(returns_pct$r_DAX),
    median(returns_pct$r_DAX),
    sd(returns_pct$r_DAX),
    min(returns_pct$r_DAX),
    max(returns_pct$r_DAX),
    skewness(returns_pct$r_DAX),
    kurtosis(returns_pct$r_DAX) - 3
  ),
  SP500 = c(
    mean(returns_pct$r_SP500),
    median(returns_pct$r_SP500),
    sd(returns_pct$r_SP500),
    min(returns_pct$r_SP500),
    max(returns_pct$r_SP500),
    skewness(returns_pct$r_SP500),
    kurtosis(returns_pct$r_SP500) - 3
  )
)
summary_returns[, c("DAX", "SP500")] <-
  round(summary_returns[, c("DAX", "SP500")], 4)

print(summary_returns)

# -------------------------------------------------------------
# Visualizzazione serie temporali rendimenti
# -------------------------------------------------------------

plot(data$Date, returns_pct$r_DAX, type = "l", col = "blue", lwd = 0.5,
     xlab = "Data", ylab = "Rendimento DAX (%)", 
     main = "Rendimenti Logaritmici DAX")
abline(h = 0, col = "gray", lty = 2)
grid()

plot(data$Date, returns_pct$r_SP500, type = "l", col = "red", lwd = 0.5,
     xlab = "Data", ylab = "Rendimento S&P 500 (%)",
     main = "Rendimenti Logaritmici S&P 500")
abline(h = 0, col = "gray", lty = 2)
grid()

# -------------------------------------------------------------
# Istogrammi con overlay normale
# -------------------------------------------------------------

# DAX
hist(returns_pct$r_DAX, breaks = 50, freq = FALSE, 
     col = rgb(0, 0, 1, 0.5), border = "white",
     xlab = "Rendimento DAX (%)", ylab = "Densità",
     main = "Distribuzione Rendimenti DAX")

x_seq <- seq(min(returns_pct$r_DAX), max(returns_pct$r_DAX), length.out = 200)
lines(x_seq, dnorm(x_seq, mean = mean(returns_pct$r_DAX), 
                   sd = sd(returns_pct$r_DAX)),
      col = "red", lwd = 2)
legend("topleft", legend = c("Empirica", "Normale"), 
       fill = c(rgb(0, 0, 1, 0.5), NA), 
       border = c("white", NA),
       lty = c(NA, 1), col = c(NA, "red"), lwd = c(NA, 2),
       bty = "n")
grid()

# S&P500
hist(returns_pct$r_SP500, breaks = 50, freq = FALSE,
     col = rgb(1, 0, 0, 0.5), border = "white",
     xlab = "Rendimento S&P 500 (%)", ylab = "Densità",
     main = "Distribuzione Rendimenti S&P 500")

x_seq <- seq(min(returns_pct$r_SP500), max(returns_pct$r_SP500), length.out = 200)
lines(x_seq, dnorm(x_seq, mean = mean(returns_pct$r_SP500), 
                   sd = sd(returns_pct$r_SP500)),
      col = "darkred", lwd = 2)
legend("topleft", legend = c("Empirica", "Normale"), 
       fill = c(rgb(1, 0, 0, 0.5), NA),
       border = c("white", NA),
       lty = c(NA, 1), col = c(NA, "darkred"), lwd = c(NA, 2),
       bty = "n")
grid()

# -------------------------------------------------------------
# Q-Q plots
# -------------------------------------------------------------

# DAX
qqnorm(returns_pct$r_DAX, pch = 16, col = rgb(0, 0, 1, 0.4),
       xlab = "Quantili Teorici (Normale)", 
       ylab = "Quantili Empirici DAX (%)",
       main = "Q-Q Plot DAX vs Normale")
qqline(returns_pct$r_DAX, col = "red", lwd = 2)
grid()

# S&P500
qqnorm(returns_pct$r_SP500, pch = 16, col = rgb(1, 0, 0, 0.4),
       xlab = "Quantili Teorici (Normale)", 
       ylab = "Quantili Empirici S&P 500 (%)",
       main = "Q-Q Plot S&P 500 vs Normale")
qqline(returns_pct$r_SP500, col = "darkred", lwd = 2)
grid()

# -------------------------------------------------------------
# Scatterplot rendimenti
# -------------------------------------------------------------

plot(returns_pct$r_DAX, returns_pct$r_SP500, 
     pch = 16, col = rgb(0, 0, 1, 0.3), cex = 0.8,
     xlab = "Rendimento DAX (%)", ylab = "Rendimento S&P 500 (%)",
     main = "Scatterplot Rendimenti DAX vs S&P 500")
abline(lm(data$r_SP500 ~ data$r_DAX), col = "red", lwd = 2)
grid()

# -------------------------------------------------------------
# Misure di correlazione
# -------------------------------------------------------------
cat("\n=== Misure di Correlazione ===\n\n")

rho_pearson <- cor(data$r_DAX, data$r_SP500, method = "pearson")
tau_kendall <- cor(data$r_DAX, data$r_SP500, method = "kendall")
rho_spearman <- cor(data$r_DAX, data$r_SP500, method = "spearman")

cat("Pearson ρ:   ", round(rho_pearson, 4), "\n")
cat("Kendall τ:   ", round(tau_kendall, 4), "\n")
cat("Spearman ρ_S:", round(rho_spearman, 4), "\n\n")


###############################################################################
# SEZIONE 3.3: TRASFORMAZIONE A PSEUDO-OSSERVAZIONI
###############################################################################

# -------------------------------------------------------------
# Probability integral transform empirico
# -------------------------------------------------------------
u <- rank(data$r_DAX) / (nrow(data) + 1)
v <- rank(data$r_SP500) / (nrow(data) + 1)

cat("=== Pseudo-osservazioni Create ===\n")
cat("Dimensione:", length(u), "\n")
cat("Intervallo u: [", round(min(u), 4), ",", round(max(u), 4), "]\n")
cat("Intervallo v: [", round(min(v), 4), ",", round(max(v), 4), "]\n\n")

# -------------------------------------------------------------
# Visualizzazione pseudo-osservazioni
# -------------------------------------------------------------

plot(u, v, pch = 16, col = rgb(0, 0, 1, 0.3), cex = 0.8,
     xlab = "u (DAX)", ylab = "v (S&P 500)",
     main = "Pseudo-osservazioni",
     xlim = c(0, 1), ylim = c(0, 1))
grid()


###############################################################################
# FUNZIONI AUSILIARIE PER INFERENZA BAYESIANA
###############################################################################

# -------------------------------------------------------------
# Utility: trasformazioni
# -------------------------------------------------------------
clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

tau_from_z <- function(z) tanh(z)
z_from_tau <- function(tau) atanh(tau)

# -------------------------------------------------------------
# Supporto tau per famiglia di copula
# -------------------------------------------------------------
tau_support <- function(family = c("clayton", "gumbel", "frank")) {
  family <- match.arg(family)
  if (family %in% c("clayton", "gumbel")) return(c(0, 1))
  return(c(-1, 1))
}

# -------------------------------------------------------------
# Prior uniforme su tau
# -------------------------------------------------------------
log_prior_tau <- function(tau, family = c("clayton", "gumbel", "frank")) {
  family <- match.arg(family)
  if (!is.finite(tau)) return(-Inf)
  
  s <- tau_support(family)
  lo <- s[1]; hi <- s[2]
  
  if (tau <= lo || tau >= hi) return(-Inf)
  
  if (family %in% c("clayton", "gumbel")) return(0)
  return(-log(2))
}

# -------------------------------------------------------------
# Creazione oggetto copula
# -------------------------------------------------------------
make_copula <- function(theta, family = c("clayton", "gumbel", "frank")) {
  family <- match.arg(family)
  switch(
    family,
    clayton = claytonCopula(param = theta, dim = 2),
    gumbel  = gumbelCopula(param  = theta, dim = 2),
    frank   = frankCopula(param   = theta, dim = 2)
  )
}

# -------------------------------------------------------------
# Trasformazione tau -> theta via iTau
# -------------------------------------------------------------
theta_from_tau <- function(tau, family = c("clayton", "gumbel", "frank")) {
  family <- match.arg(family)
  
  cop0 <- switch(
    family,
    clayton = claytonCopula(param = 1, dim = 2),
    gumbel  = gumbelCopula(param  = 2, dim = 2),
    frank   = frankCopula(param   = 5, dim = 2)
  )
  
  th <- tryCatch(iTau(cop0, tau), error = function(e) NA_real_)
  th
}

# -------------------------------------------------------------
# Log-verosimiglianza copula
# -------------------------------------------------------------
loglik_copula <- function(u, v, theta, family = c("clayton", "gumbel", "frank")) {
  family <- match.arg(family)
  if (!is.finite(theta)) return(-Inf)
  
  cop <- make_copula(theta, family)
  x <- cbind(u, v)
  
  ll <- tryCatch(sum(dCopula(x, copula = cop, log = TRUE)),
                 error = function(e) -Inf)
  
  if (!is.finite(ll)) -Inf else ll
}

# -------------------------------------------------------------
# Log-posterior in spazio z = atanh(tau)
# -------------------------------------------------------------
logpost_z <- function(z, u, v, family = c("clayton", "gumbel", "frank")) {
  family <- match.arg(family)
  
  tau <- tau_from_z(z)
  
  lp <- log_prior_tau(tau, family)
  if (!is.finite(lp)) return(-Inf)
  
  theta <- theta_from_tau(tau, family)
  if (!is.finite(theta)) return(-Inf)
  
  ll <- loglik_copula(u, v, theta, family)
  if (!is.finite(ll)) return(-Inf)
  
  log_jac <- log(1 - tau^2)
  
  ll + lp + log_jac
}

# -------------------------------------------------------------
# Algoritmo Metropolis-Hastings Random Walk
# -------------------------------------------------------------
mh_rw_1d <- function(logpost_fun, init,
                     n_iter = 50000, burn_in = 5000,
                     sigma = 0.6,
                     adapt = TRUE, target = 0.25, n_adapt = 3000,
                     sigma_min = 0.05, sigma_max = 2.0,
                     verbose_every = 10000) {
  
  n_tot <- n_iter + burn_in
  chain <- numeric(n_tot)
  chain[1] <- init
  
  lp_curr <- logpost_fun(init)
  if (!is.finite(lp_curr)) stop("Valore iniziale: log-posterior non finito.")
  
  acc <- 0L
  
  for (t in 2:n_tot) {
    prop <- rnorm(1, mean = chain[t - 1], sd = sigma)
    lp_prop <- logpost_fun(prop)
    
    log_alpha <- lp_prop - lp_curr
    
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      chain[t] <- prop
      lp_curr <- lp_prop
      acc <- acc + 1L
    } else {
      chain[t] <- chain[t - 1]
    }
    
    # Adattamento sigma durante burn-in
    if (adapt && t <= min(n_adapt, burn_in)) {
      ar <- acc / t
      if (ar > target) sigma <- sigma * 1.02 else sigma <- sigma * 0.98
      sigma <- min(max(sigma, sigma_min), sigma_max)
    }
    
    if (!is.null(verbose_every) && verbose_every > 0 && t %% verbose_every == 0) {
      cat("  Iter", t, "/", n_tot,
          "- Acc:", round(acc / t, 3),
          "- sigma:", round(sigma, 4), "\n")
    }
  }
  
  post_z <- chain[(burn_in + 1):n_tot]
  
  list(
    chain_z = chain,
    post_z = post_z,
    acc_rate = acc / n_tot,
    sigma_last = sigma
  )
}

# -------------------------------------------------------------
# Wrapper MCMC per copule
# -------------------------------------------------------------
mcmc_copula_tau <- function(u, v,
                            family = c("clayton", "gumbel", "frank"),
                            n_iter = 50000, burn_in = 5000,
                            sigma = 0.6, target = 0.25,
                            seed = NULL,
                            verbose_every = 10000) {
  
  family <- match.arg(family)
  if (!is.null(seed)) set.seed(seed)
  
  # Inizializzazione da Kendall empirico
  tau_init <- suppressWarnings(cor(u, v, method = "kendall"))
  if (!is.finite(tau_init)) tau_init <- 0
  
  s <- tau_support(family)
  if (family %in% c("clayton", "gumbel")) {
    tau_init <- clamp(tau_init, 0.05, 0.95)
  } else {
    tau_init <- clamp(tau_init, -0.95, 0.95)
  }
  
  z_init <- z_from_tau(tau_init)
  lp <- function(z) logpost_z(z, u = u, v = v, family = family)
  
  cat("\n=== MCMC", toupper(family), "===\n")
  cat("  tau init =", round(tau_init, 4), "| z init =", round(z_init, 4), "\n")
  
  fit <- mh_rw_1d(
    logpost_fun = lp,
    init = z_init,
    n_iter = n_iter,
    burn_in = burn_in,
    sigma = sigma,
    adapt = TRUE,
    target = target,
    n_adapt = min(3000, burn_in),
    verbose_every = verbose_every
  )
  
  tau_post <- tau_from_z(fit$post_z)
  theta_post <- vapply(tau_post, theta_from_tau, numeric(1), family = family)
  
  m_tau <- mcmc(tau_post)
  
  list(
    family = family,
    tau_init = tau_init,
    tau_samples = tau_post,
    theta_samples = theta_post,
    chain_z = fit$chain_z,
    acc_rate = fit$acc_rate,
    sigma_last = fit$sigma_last,
    ess_tau = as.numeric(effectiveSize(m_tau)),
    hpd_tau = HPDinterval(m_tau, prob = 0.95)
  )
}

# -------------------------------------------------------------
# Tail dependence coefficients
# -------------------------------------------------------------
tail_dependence <- function(theta_samples, family = c("clayton", "gumbel", "frank")) {
  family <- match.arg(family)
  
  if (family == "clayton") {
    lamL <- 2^(-1 / theta_samples)
    return(list(lambda_L = lamL, lambda_U = rep(0, length(theta_samples))))
  }
  
  if (family == "gumbel") {
    lamU <- 2 - 2^(1 / theta_samples)
    return(list(lambda_L = rep(0, length(theta_samples)), lambda_U = lamU))
  }
  
  # Frank: no tail dependence
  return(list(lambda_L = rep(0, length(theta_samples)), 
              lambda_U = rep(0, length(theta_samples))))
}


###############################################################################
# SEZIONE 3.4: INFERENZA BAYESIANA PER COPULE ARCHIMEDEE
###############################################################################

# -------------------------------------------------------------
# Esecuzione MCMC per le tre famiglie
# -------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("INFERENZA BAYESIANA SU COPULE ARCHIMEDEE\n")
cat(strrep("=", 70), "\n")

set.seed(42)
fit_clayton <- mcmc_copula_tau(u, v, family = "clayton",
                               n_iter = 50000, burn_in = 5000,
                               sigma = 0.6, target = 0.25,
                               seed = 42)

set.seed(43)
fit_gumbel  <- mcmc_copula_tau(u, v, family = "gumbel",
                               n_iter = 50000, burn_in = 5000,
                               sigma = 0.6, target = 0.25,
                               seed = 43)

set.seed(44)
fit_frank   <- mcmc_copula_tau(u, v, family = "frank",
                               n_iter = 50000, burn_in = 5000,
                               sigma = 0.6, target = 0.25,
                               seed = 44)

# -------------------------------------------------------------
# Statistiche riassuntive
# -------------------------------------------------------------

summarize_fit <- function(fit) {
  tau_mean <- mean(fit$tau_samples)
  tau_med  <- median(fit$tau_samples)
  tau_ci   <- quantile(fit$tau_samples, c(0.025, 0.975))
  
  th_mean  <- mean(fit$theta_samples)
  th_ci    <- quantile(fit$theta_samples, c(0.025, 0.975))
  
  td <- tail_dependence(fit$theta_samples, fit$family)
  
  lamL_mean <- mean(td$lambda_L)
  lamU_mean <- mean(td$lambda_U)
  
  lamL_ci <- quantile(td$lambda_L, c(0.025, 0.975))
  lamU_ci <- quantile(td$lambda_U, c(0.025, 0.975))
  
  list(
    tau_mean = tau_mean,
    tau_median = tau_med,
    tau_ci_low = unname(tau_ci[1]),
    tau_ci_high = unname(tau_ci[2]),
    theta_mean = th_mean,
    theta_ci_low = unname(th_ci[1]),
    theta_ci_high = unname(th_ci[2]),
    lambda_L_mean = lamL_mean,
    lambda_L_ci_low = unname(lamL_ci[1]),
    lambda_L_ci_high = unname(lamL_ci[2]),
    lambda_U_mean = lamU_mean,
    lambda_U_ci_low = unname(lamU_ci[1]),
    lambda_U_ci_high = unname(lamU_ci[2]),
    ESS = fit$ess_tau,
    Acc_Rate = fit$acc_rate
  )
}

S_clayton <- summarize_fit(fit_clayton)
S_gumbel  <- summarize_fit(fit_gumbel)
S_frank   <- summarize_fit(fit_frank)

# -------------------------------------------------------------
# Tabella comparativa
# -------------------------------------------------------------
cat("\n=== Confronto tra Famiglie di Copule ===\n\n")

comparison_table <- data.frame(
  Famiglia   = c("Clayton", "Gumbel", "Frank"),
  tau_mean   = c(S_clayton$tau_mean,   S_gumbel$tau_mean,   S_frank$tau_mean),
  tau_median = c(S_clayton$tau_median, S_gumbel$tau_median, S_frank$tau_median),
  tau_ci_low = c(S_clayton$tau_ci_low, S_gumbel$tau_ci_low, S_frank$tau_ci_low),
  tau_ci_high= c(S_clayton$tau_ci_high,S_gumbel$tau_ci_high,S_frank$tau_ci_high),
  theta_mean = c(S_clayton$theta_mean, S_gumbel$theta_mean, S_frank$theta_mean),
  lambda_L   = c(S_clayton$lambda_L_mean, S_gumbel$lambda_L_mean, S_frank$lambda_L_mean),
  lambda_U   = c(S_clayton$lambda_U_mean, S_gumbel$lambda_U_mean, S_frank$lambda_U_mean),
  ESS        = c(S_clayton$ESS, S_gumbel$ESS, S_frank$ESS),
  Acc_Rate   = c(S_clayton$Acc_Rate, S_gumbel$Acc_Rate, S_frank$Acc_Rate)
)

comparison_table <- data.frame(
  Famiglia   = c("Clayton", "Gumbel", "Frank"),
  tau_mean   = c(S_clayton$tau_mean,   S_gumbel$tau_mean,   S_frank$tau_mean),
  tau_median = c(S_clayton$tau_median, S_gumbel$tau_median, S_frank$tau_median),
  tau_ci_low = c(S_clayton$tau_ci_low, S_gumbel$tau_ci_low, S_frank$tau_ci_low),
  tau_ci_high= c(S_clayton$tau_ci_high,S_gumbel$tau_ci_high,S_frank$tau_ci_high),
  theta_mean = c(S_clayton$theta_mean, S_gumbel$theta_mean, S_frank$theta_mean),
  lambda_L   = c(S_clayton$lambda_L_mean, S_gumbel$lambda_L_mean, S_frank$lambda_L_mean),
  lambda_U   = c(S_clayton$lambda_U_mean, S_gumbel$lambda_U_mean, S_frank$lambda_U_mean),
  ESS        = c(S_clayton$ESS, S_gumbel$ESS, S_frank$ESS),
  Acc_Rate   = c(S_clayton$Acc_Rate, S_gumbel$Acc_Rate, S_frank$Acc_Rate),
  check.names = FALSE
)


num_cols <- sapply(comparison_table, is.numeric)
comparison_table[num_cols] <- lapply(comparison_table[num_cols], round, 4)

print(comparison_table)




cat("\nHPD 95% su tau:\n")
cat("Clayton:", round(as.numeric(fit_clayton$hpd_tau[1]), 4), 
    round(as.numeric(fit_clayton$hpd_tau[2]), 4), "\n")
cat("Gumbel :", round(as.numeric(fit_gumbel$hpd_tau[1]), 4),  
    round(as.numeric(fit_gumbel$hpd_tau[2]), 4),  "\n")
cat("Frank  :", round(as.numeric(fit_frank$hpd_tau[1]), 4),   
    round(as.numeric(fit_frank$hpd_tau[2]), 4),   "\n\n")


###############################################################################
# SEZIONE 3.5: MODEL COMPARISON (DIC e WAIC)
###############################################################################

cat(strrep("=", 70), "\n")
cat("MODEL COMPARISON: DIC e WAIC\n")
cat(strrep("=", 70), "\n\n")

# -------------------------------------------------------------
# Funzione per calcolare DIC
# -----------------------------------------------------------



log_mean_exp <- function(x) {
  # x: vettore di log-likelihood (uno per campione) per una osservazione
  m <- max(x)
  if (!is.finite(m)) return(-Inf)
  m + log(mean(exp(x - m)))
}
compute_WAIC_fast <- function(u, v, tau_samples, family, S_use = 2000, seed = 1) {
  n <- length(u)
  X <- cbind(u, v)
  
  if (!is.null(S_use) && length(tau_samples) > S_use) {
    set.seed(seed)
    tau_samples <- tau_samples[sample.int(length(tau_samples), S_use)]
  }
  
  S <- length(tau_samples)
  loglik_matrix <- matrix(NA_real_, nrow = S, ncol = n)
  
  n_failed <- 0
  for (s in seq_len(S)) {
    theta <- theta_from_tau(tau_samples[s], family)
    if (!is.finite(theta)) stop("WAIC: theta non finito")
    
    cop <- make_copula(theta, family)
    ll <- tryCatch(dCopula(X, copula = cop, log = TRUE),
                   error = function(e) rep(NA_real_, n))
    if (any(!is.finite(ll))) n_failed <- n_failed + 1
    loglik_matrix[s, ] <- ll
  }
  
  if (any(!is.finite(loglik_matrix))) {
    stop("WAIC: loglik non finiti. Risolvi instabilità numerica invece di filtrare per-colonna.")
  }
  
  lppd <- sum(apply(loglik_matrix, 2, log_mean_exp))
  pW   <- sum(apply(loglik_matrix, 2, var))
  WAIC <- -2 * (lppd - pW)
  
  list(WAIC = WAIC, p_WAIC = pW, lppd = lppd, S_used = S, n_failed = n_failed)
}
compute_DIC_fast <- function(u, v, tau_samples, family, S_use = 5000, seed = 1) {
  X <- cbind(u, v)
  
  if (!is.null(S_use) && length(tau_samples) > S_use) {
    set.seed(seed)
    idx <- sample.int(length(tau_samples), S_use)
    tau_samples <- tau_samples[idx]
  }
  
  # loglik per draw
  loglik_samples <- vapply(tau_samples, function(tau) {
    theta <- theta_from_tau(tau, family)
    if (!is.finite(theta)) return(NA_real_)
    cop <- make_copula(theta, family)
    tryCatch(sum(dCopula(X, copula = cop, log = TRUE)), error = function(e) NA_real_)
  }, numeric(1))
  
  loglik_samples <- loglik_samples[is.finite(loglik_samples)]
  if (length(loglik_samples) < 10) stop("Troppi loglik non finiti: controlla theta_from_tau/make_copula.")
  
  # Dbar
  Dbar <- -2 * mean(loglik_samples)
  
  # theta_bar: meglio media in scala theta
  theta_samples <- vapply(tau_samples, function(tau) theta_from_tau(tau, family), numeric(1))
  theta_bar <- mean(theta_samples[is.finite(theta_samples)])
  
  cop_bar <- make_copula(theta_bar, family)
  loglik_bar <- tryCatch(sum(dCopula(X, copula = cop_bar, log = TRUE)), error = function(e) NA_real_)
  if (!is.finite(loglik_bar)) stop("loglik_bar non finito: controlla theta_bar e make_copula.")
  
  Dhat <- -2 * loglik_bar
  
  p_DIC <- Dbar - Dhat
  DIC <- Dbar + p_DIC  # = 2*Dbar - Dhat
  
  list(DIC = DIC, p_DIC = p_DIC, Dbar = Dbar, Dhat = Dhat, S_used = length(loglik_samples))
}

cat("Calcolo DIC e WAIC (FAST)...\n\n")

DIC_clayton  <- compute_DIC_fast(u, v, fit_clayton$tau_samples, "clayton", S_use = 10000, seed = 1)
WAIC_clayton <- compute_WAIC_fast(u, v, fit_clayton$tau_samples, "clayton", S_use = 10000, seed = 1)

DIC_gumbel   <- compute_DIC_fast(u, v, fit_gumbel$tau_samples, "gumbel",  S_use = 10000, seed = 1)
WAIC_gumbel  <- compute_WAIC_fast(u, v, fit_gumbel$tau_samples, "gumbel",  S_use = 10000, seed = 1)

DIC_frank    <- compute_DIC_fast(u, v, fit_frank$tau_samples, "frank",    S_use = 10000, seed = 1)
WAIC_frank   <- compute_WAIC_fast(u, v, fit_frank$tau_samples, "frank",    S_use = 10000, seed = 1)

cat("=== Model Comparison ===\n\n")

model_comparison <- data.frame(
  Famiglia = c("Clayton", "Gumbel", "Frank"),
  DIC = c(DIC_clayton$DIC, DIC_gumbel$DIC, DIC_frank$DIC),
  p_DIC = c(DIC_clayton$p_DIC, DIC_gumbel$p_DIC, DIC_frank$p_DIC),
  WAIC = c(WAIC_clayton$WAIC, WAIC_gumbel$WAIC, WAIC_frank$WAIC),
  p_WAIC = c(WAIC_clayton$p_WAIC, WAIC_gumbel$p_WAIC, WAIC_frank$p_WAIC)
)
num_cols <- sapply(model_comparison, is.numeric)
model_comparison[num_cols] <- lapply(model_comparison[num_cols], round, 2)
print(model_comparison)

#print(round(model_comparison, 2))

# Identifica il modello migliore
best_DIC <- which.min(model_comparison$DIC)
best_WAIC <- which.min(model_comparison$WAIC)

cat("\nModello migliore (DIC):", model_comparison$Famiglia[best_DIC], "\n")
cat("Modello migliore (WAIC):", model_comparison$Famiglia[best_WAIC], "\n\n")

# ==========================================================
# DIAGNOSTICA MCMC (Trace + Densità + ACF)
# ==========================================================

plot_mcmc_basic <- function(tau_chain, stats, family_name = "", col_line = "steelblue", lag_max = 50) {
  tau_mean <- stats$tau_mean
  ci_low   <- stats$tau_ci_low
  ci_high  <- stats$tau_ci_high
  
  # 1) Trace plot
  plot(tau_chain, type = "l", col = col_line,
       xlab = "Iterazione (post burn-in)", ylab = expression(tau),
       main = paste("Trace plot -", family_name))
  abline(h = tau_mean, col = "red", lwd = 2, lty = 2)
  grid()
  
  # 2) Densità a posteriori
  d <- density(tau_chain)
  plot(d, col = col_line, lwd = 2,
       main = paste("Densità a posteriori di", expression(tau), "-", family_name),
       xlab = expression(tau), ylab = "Densità")
  abline(v = tau_mean, col = "red", lwd = 2, lty = 2)
  abline(v = c(ci_low, ci_high), col = "red", lwd = 1.5, lty = 3)
  grid()
  
  # 3) ACF
  acf(tau_chain, main = paste("ACF -", family_name), lag.max = lag_max)
}

# Esegui per tutte le famiglie
plot_mcmc_basic(fit_clayton$tau_samples, S_clayton, "Clayton", col_line = "blue",  lag_max = 50)
plot_mcmc_basic(fit_gumbel$tau_samples,  S_gumbel,  "Gumbel",  col_line = "red",   lag_max = 50)
plot_mcmc_basic(fit_frank$tau_samples,   S_frank,   "Frank",   col_line = "green", lag_max = 50)


# -------------------------------------------------------------
# Geweke diagnostics
# -------------------------------------------------------------

cat("=== Geweke Diagnostics ===\n\n")

geweke_clayton <- geweke.diag(mcmc(fit_clayton$tau_samples))
geweke_gumbel  <- geweke.diag(mcmc(fit_gumbel$tau_samples))
geweke_frank   <- geweke.diag(mcmc(fit_frank$tau_samples))

cat("Clayton - Z-score:", round(geweke_clayton$z, 4), "\n")
cat("Gumbel  - Z-score:", round(geweke_gumbel$z, 4), "\n")
cat("Frank   - Z-score:", round(geweke_frank$z, 4), "\n")
cat("\n(Valori assoluti < 2 indicano convergenza)\n\n")

# -------------------------------------------------------------
# Heidelberger-Welch diagnostics
# -------------------------------------------------------------

cat("=== Heidelberger-Welch Diagnostics ===\n\n")

hw_clayton <- heidel.diag(mcmc(fit_clayton$tau_samples))
hw_gumbel  <- heidel.diag(mcmc(fit_gumbel$tau_samples))
hw_frank   <- heidel.diag(mcmc(fit_frank$tau_samples))

cat("Clayton:\n")
print(hw_clayton)
cat("\nGumbel:\n")
print(hw_gumbel)
cat("\nFrank:\n")
print(hw_frank)
cat("\n")


###############################################################################
# SEZIONE 3.6: ANALISI PER SOTTO-PERIODI (VERSIONE CORRETTA)
###############################################################################

cat("\n", strrep("=", 70), "\n")
cat("ANALISI PER SOTTO-PERIODI\n")
cat(strrep("=", 70), "\n\n")

# -------------------------------------------------------------
# Definizione dei sotto-periodi
# -------------------------------------------------------------

# IMPORTANTE: Assicurati che data$Date sia in formato Date
if (!inherits(data$Date, "Date")) {
  data$Date <- as.Date(data$Date)
}

# Definisci i cut-points per i sotto-periodi
periods <- list(
  list(name = "Pre-crisi (2005-2007)", 
       start = as.Date("2005-01-01"), 
       end = as.Date("2007-12-31")),
  
  list(name = "Crisi GFC (2008-2009)", 
       start = as.Date("2008-01-01"), 
       end = as.Date("2009-12-31")),
  
  list(name = "Ripresa e crisi del debito (2010-2012)", 
       start = as.Date("2010-01-01"), 
       end = as.Date("2012-12-31")),
  
  list(name = "Bull market (2013-2019)", 
       start = as.Date("2013-01-01"), 
       end = as.Date("2019-12-31")),
  
  list(name = "COVID-19 (2020)", 
       start = as.Date("2020-01-01"), 
       end = as.Date("2020-12-31")),
  
  list(name = "Post-COVID (2021-2024)", 
       start = as.Date("2021-01-01"), 
       end = as.Date("2024-12-31"))
)

cat("Sotto-periodi definiti:\n")
for (i in seq_along(periods)) {
  cat(sprintf("%d. %s: %s - %s\n", 
              i, 
              periods[[i]]$name,
              periods[[i]]$start,
              periods[[i]]$end))
}
cat("\n")

# Verifica i dati
cat("Verifica dataset:\n")
cat(sprintf("  Classe Date: %s\n", class(data$Date)))
cat(sprintf("  Range date: %s - %s\n", min(data$Date), max(data$Date)))
cat(sprintf("  N totale osservazioni: %d\n\n", nrow(data)))

# -------------------------------------------------------------
# Funzione per stimare copula su un sotto-periodo (CORRETTA)
# -------------------------------------------------------------

estimate_subperiod <- function(data_full, period_info, 
                               family = "clayton",
                               n_iter = 30000, 
                               burn_in = 3000,
                               sigma = 0.6,
                               seed = NULL) {
  
  # Converti esplicitamente le date
  start_date <- as.Date(period_info$start)
  end_date <- as.Date(period_info$end)
  
  # Estrai dati del sotto-periodo
  mask <- data_full$Date >= start_date & data_full$Date <= end_date
  
  data_sub <- data_full[mask, ]
  n_obs <- nrow(data_sub)
  
  cat(sprintf("\n--- %s ---\n", period_info$name))
  cat(sprintf("Periodo: %s - %s\n", start_date, end_date))
  cat(sprintf("Osservazioni trovate: %d\n", n_obs))
  
  if (n_obs < 50) {
    warning(sprintf("Periodo '%s': solo %d osservazioni. Saltato.", 
                    period_info$name, n_obs))
    return(NULL)
  }
  
  cat(sprintf("Range effettivo: %s - %s\n", min(data_sub$Date), max(data_sub$Date)))
  
  # Calcola pseudo-osservazioni per il sotto-periodo
  u_sub <- rank(data_sub$r_DAX) / (n_obs + 1)
  v_sub <- rank(data_sub$r_SP500) / (n_obs + 1)
  
  # Kendall empirico
  tau_empirical <- cor(data_sub$r_DAX, data_sub$r_SP500, method = "kendall")
  cat(sprintf("Tau di Kendall empirico: %.4f\n", tau_empirical))
  
  # Stima bayesiana
  cat("Esecuzione MCMC...\n")
  fit <- mcmc_copula_tau(
    u = u_sub, 
    v = v_sub,
    family = family,
    n_iter = n_iter,
    burn_in = burn_in,
    sigma = sigma,
    target = 0.25,
    seed = seed,
    verbose_every = NULL
  )
  
  # Calcola statistiche
  tau_mean <- mean(fit$tau_samples)
  tau_median <- median(fit$tau_samples)
  tau_ci <- quantile(fit$tau_samples, c(0.025, 0.975))
  
  theta_mean <- mean(fit$theta_samples)
  theta_ci <- quantile(fit$theta_samples, c(0.025, 0.975))
  
  # Tail dependence
  td <- tail_dependence(fit$theta_samples, family)
  
  lambda_L_mean <- mean(td$lambda_L)
  lambda_L_ci <- quantile(td$lambda_L, c(0.025, 0.975))
  
  lambda_U_mean <- mean(td$lambda_U)
  lambda_U_ci <- quantile(td$lambda_U, c(0.025, 0.975))
  
  cat(sprintf("Risultati:\n"))
  cat(sprintf("  E[tau]: %.3f [%.3f, %.3f]\n", 
              tau_mean, tau_ci[1], tau_ci[2]))
  cat(sprintf("  E[theta]: %.3f [%.3f, %.3f]\n", 
              theta_mean, theta_ci[1], theta_ci[2]))
  cat(sprintf("  E[lambda_L]: %.3f [%.3f, %.3f]\n", 
              lambda_L_mean, lambda_L_ci[1], lambda_L_ci[2]))
  cat(sprintf("  E[lambda_U]: %.3f [%.3f, %.3f]\n", 
              lambda_U_mean, lambda_U_ci[1], lambda_U_ci[2]))
  cat(sprintf("  ESS: %.0f | Acc rate: %.3f\n", fit$ess_tau, fit$acc_rate))
  
  list(
    period_name = period_info$name,
    n_obs = n_obs,
    start_date = start_date,
    end_date = end_date,
    tau_empirical = tau_empirical,
    tau_mean = tau_mean,
    tau_median = tau_median,
    tau_ci_low = as.numeric(tau_ci[1]),
    tau_ci_high = as.numeric(tau_ci[2]),
    theta_mean = theta_mean,
    theta_ci_low = as.numeric(theta_ci[1]),
    theta_ci_high = as.numeric(theta_ci[2]),
    lambda_L_mean = lambda_L_mean,
    lambda_L_ci_low = as.numeric(lambda_L_ci[1]),
    lambda_L_ci_high = as.numeric(lambda_L_ci[2]),
    lambda_U_mean = lambda_U_mean,
    lambda_U_ci_low = as.numeric(lambda_U_ci[1]),
    lambda_U_ci_high = as.numeric(lambda_U_ci[2]),
    ess = fit$ess_tau,
    acc_rate = fit$acc_rate,
    tau_samples = fit$tau_samples,
    theta_samples = fit$theta_samples
  )
}

# -------------------------------------------------------------
# Stima Clayton per tutti i sotto-periodi
# -------------------------------------------------------------

cat("\n", strrep("-", 70), "\n")
cat("STIMA COPULA DI CLAYTON PER SOTTO-PERIODI\n")
cat(strrep("-", 70), "\n")

set.seed(123)
results_clayton <- list()

for (i in seq_along(periods)) {
  res <- estimate_subperiod(
    data_full = data,
    period_info = periods[[i]],
    family = "clayton",
    n_iter = 30000,
    burn_in = 3000,
    sigma = 0.6,
    seed = 100 + i
  )
  
  if (!is.null(res)) {
    results_clayton[[i]] <- res
  }
}

# Rimuovi eventuali NULL
results_clayton <- results_clayton[!sapply(results_clayton, is.null)]

cat(sprintf("\n\nPeriodi analizzati con successo: %d\n", length(results_clayton)))

# Verifica che abbiamo risultati
if (length(results_clayton) == 0) {
  cat("ERRORE: Nessun sotto-periodo analizzato con successo!\n")
  cat("Verifica i dati e le date.\n")
} else {
  
  # -------------------------------------------------------------
  # Tabella riassuntiva Clayton
  # -------------------------------------------------------------
  
  cat("\n", strrep("=", 70), "\n")
  cat("TABELLA RIASSUNTIVA: COPULA DI CLAYTON\n")
  cat(strrep("=", 70), "\n\n")
  
  # Costruisci la tabella con conversione esplicita a numeric
  table_clayton <- data.frame(
    Periodo = sapply(results_clayton, function(x) x$period_name),
    n = as.numeric(sapply(results_clayton, function(x) x$n_obs)),
    tau_hat = as.numeric(sapply(results_clayton, function(x) x$tau_mean)),
    tau_ci_low = as.numeric(sapply(results_clayton, function(x) x$tau_ci_low)),
    tau_ci_high = as.numeric(sapply(results_clayton, function(x) x$tau_ci_high)),
    lambda_L_hat = as.numeric(sapply(results_clayton, function(x) x$lambda_L_mean)),
    lambda_L_ci_low = as.numeric(sapply(results_clayton, function(x) x$lambda_L_ci_low)),
    lambda_L_ci_high = as.numeric(sapply(results_clayton, function(x) x$lambda_L_ci_high)),
    stringsAsFactors = FALSE
  )
  
  # Formatta intervalli di credibilità
  table_clayton$tau_IC <- sprintf("[%.3f, %.3f]", 
                                  table_clayton$tau_ci_low, 
                                  table_clayton$tau_ci_high)
  
  table_clayton$lambda_L_IC <- sprintf("[%.3f, %.3f]", 
                                       table_clayton$lambda_L_ci_low, 
                                       table_clayton$lambda_L_ci_high)
  
  # Tabella finale per output LaTeX
  table_final <- data.frame(
    Periodo = table_clayton$Periodo,
    n = table_clayton$n,
    tau_hat = round(table_clayton$tau_hat, 3),
    `IC_95_tau` = table_clayton$tau_IC,
    lambda_L_hat = round(table_clayton$lambda_L_hat, 3),
    `IC_95_lambda_L` = table_clayton$lambda_L_IC,
    stringsAsFactors = FALSE
  )
  
  # Stampa tabella
  print(table_final, row.names = FALSE)
  
  # Esporta per LaTeX (opzionale)
  cat("\n\nFormato LaTeX:\n")
  for (i in 1:nrow(table_final)) {
    cat(sprintf("%s & %d & %.3f & %s & %.3f & %s \\\\\n",
                table_final$Periodo[i],
                table_final$n[i],
                table_final$tau_hat[i],
                table_final$IC_95_tau[i],
                table_final$lambda_L_hat[i],
                table_final$IC_95_lambda_L[i]))
  }
  
  # -------------------------------------------------------------
  # Stima Gumbel per periodi di crescita
  # -------------------------------------------------------------
  
  cat("\n", strrep("-", 70), "\n")
  cat("STIMA COPULA DI GUMBEL PER PERIODI DI CRESCITA\n")
  cat(strrep("-", 70), "\n")
  
  # Identifica i periodi di crescita
  growth_periods_idx <- c(4, 6)  # Bull market e Post-COVID
  
  results_gumbel <- list()
  
  for (idx in growth_periods_idx) {
    res <- estimate_subperiod(
      data_full = data,
      period_info = periods[[idx]],
      family = "gumbel",
      n_iter = 30000,
      burn_in = 3000,
      sigma = 0.6,
      seed = 200 + idx
    )
    
    if (!is.null(res)) {
      results_gumbel[[length(results_gumbel) + 1]] <- res
    }
  }
  
  # Tabella riassuntiva Gumbel
  if (length(results_gumbel) > 0) {
    cat("\n", strrep("=", 70), "\n")
    cat("TABELLA RIASSUNTIVA: COPULA DI GUMBEL (Periodi di crescita)\n")
    cat(strrep("=", 70), "\n\n")
    
    table_gumbel <- data.frame(
      Periodo = sapply(results_gumbel, function(x) x$period_name),
      n = as.numeric(sapply(results_gumbel, function(x) x$n_obs)),
      tau_hat = as.numeric(sapply(results_gumbel, function(x) x$tau_mean)),
      tau_ci_low = as.numeric(sapply(results_gumbel, function(x) x$tau_ci_low)),
      tau_ci_high = as.numeric(sapply(results_gumbel, function(x) x$tau_ci_high)),
      lambda_U_hat = as.numeric(sapply(results_gumbel, function(x) x$lambda_U_mean)),
      lambda_U_ci_low = as.numeric(sapply(results_gumbel, function(x) x$lambda_U_ci_low)),
      lambda_U_ci_high = as.numeric(sapply(results_gumbel, function(x) x$lambda_U_ci_high)),
      stringsAsFactors = FALSE
    )
    
    table_gumbel$tau_IC <- sprintf("[%.3f, %.3f]", 
                                   table_gumbel$tau_ci_low, 
                                   table_gumbel$tau_ci_high)
    
    table_gumbel$lambda_U_IC <- sprintf("[%.3f, %.3f]", 
                                        table_gumbel$lambda_U_ci_low, 
                                        table_gumbel$lambda_U_ci_high)
    
    table_gumbel_final <- data.frame(
      Periodo = table_gumbel$Periodo,
      n = table_gumbel$n,
      tau_hat = round(table_gumbel$tau_hat, 3),
      `IC_95_tau` = table_gumbel$tau_IC,
      lambda_U_hat = round(table_gumbel$lambda_U_hat, 3),
      `IC_95_lambda_U` = table_gumbel$lambda_U_IC,
      stringsAsFactors = FALSE
    )
    
    print(table_gumbel_final, row.names = FALSE)
    
    # Formato LaTeX
    cat("\n\nFormato LaTeX:\n")
    for (i in 1:nrow(table_gumbel_final)) {
      cat(sprintf("%s & %d & %.3f & %s & %.3f & %s \\\\\n",
                  table_gumbel_final$Periodo[i],
                  table_gumbel_final$n[i],
                  table_gumbel_final$tau_hat[i],
                  table_gumbel_final$IC_95_tau[i],
                  table_gumbel_final$lambda_U_hat[i],
                  table_gumbel_final$IC_95_lambda_U[i]))
    }
  }
  
  # -------------------------------------------------------------
  # Visualizzazione: Evoluzione temporale
  # -------------------------------------------------------------
  
  cat("\n", strrep("-", 70), "\n")
  cat("VISUALIZZAZIONE: Evoluzione temporale\n")
  cat(strrep("-", 70), "\n\n")
  
  # Prepara dati per il grafico
  plot_data <- data.frame(
    periodo = sapply(results_clayton, function(x) x$period_name),
    start_date = as.Date(sapply(results_clayton, function(x) as.character(x$start_date))),
    end_date = as.Date(sapply(results_clayton, function(x) as.character(x$end_date))),
    tau = as.numeric(sapply(results_clayton, function(x) x$tau_mean)),
    tau_low = as.numeric(sapply(results_clayton, function(x) x$tau_ci_low)),
    tau_high = as.numeric(sapply(results_clayton, function(x) x$tau_ci_high)),
    lambda_L = as.numeric(sapply(results_clayton, function(x) x$lambda_L_mean)),
    lambda_L_low = as.numeric(sapply(results_clayton, function(x) x$lambda_L_ci_low)),
    lambda_L_high = as.numeric(sapply(results_clayton, function(x) x$lambda_L_ci_high)),
    stringsAsFactors = FALSE
  )
  
  # Calcola punto medio per ogni periodo
  plot_data$midpoint <- plot_data$start_date + 
    as.numeric(difftime(plot_data$end_date, plot_data$start_date, units = "days")) / 2
  
  # Verifica che ci siano dati validi
  if (nrow(plot_data) > 0 && all(is.finite(plot_data$tau))) {
    
    # Grafico 1: Evoluzione di tau
    par(mfrow = c(2, 1), mar = c(5, 4, 4, 2))
    
    plot(plot_data$midpoint, plot_data$tau, 
         type = "b", pch = 19, col = "blue", lwd = 2, cex = 1.5,
         ylim = range(c(plot_data$tau_low, plot_data$tau_high)),
         xlab = "Periodo", 
         ylab = expression(tau~"(Kendall)"),
         main = "Evoluzione temporale della dipendenza (Clayton)")
    
    # Intervalli di credibilità
    arrows(plot_data$midpoint, plot_data$tau_low,
           plot_data$midpoint, plot_data$tau_high,
           angle = 90, code = 3, length = 0.05, col = "blue", lwd = 1.5)
    
    # Media intero periodo
    abline(h = mean(fit_clayton$tau_samples), col = "red", lty = 2, lwd = 2)
    
    grid()
    
    legend("topleft", 
           legend = c("Stima per sotto-periodo", "IC 95%", "Media intero periodo"),
           col = c("blue", "blue", "red"),
           lty = c(1, 1, 2),
           pch = c(19, NA, NA),
           lwd = c(2, 1.5, 2),
           bty = "n")
    
    # Grafico 2: Evoluzione di lambda_L
    plot(plot_data$midpoint, plot_data$lambda_L, 
         type = "b", pch = 19, col = "darkred", lwd = 2, cex = 1.5,
         ylim = range(c(plot_data$lambda_L_low, plot_data$lambda_L_high)),
         xlab = "Periodo", 
         ylab = expression(lambda[L]),
         main = "Evoluzione della dipendenza nella coda inferiore (Clayton)")
    
    arrows(plot_data$midpoint, plot_data$lambda_L_low,
           plot_data$midpoint, plot_data$lambda_L_high,
           angle = 90, code = 3, length = 0.05, col = "darkred", lwd = 1.5)
    
    # Media intero periodo
    td_full <- tail_dependence(fit_clayton$theta_samples, "clayton")
    abline(h = mean(td_full$lambda_L), col = "red", lty = 2, lwd = 2)
    
    grid()
    
    legend("topleft", 
           legend = c("Stima per sotto-periodo", "IC 95%", "Media intero periodo"),
           col = c("darkred", "darkred", "red"),
           lty = c(1, 1, 2),
           pch = c(19, NA, NA),
           lwd = c(2, 1.5, 2),
           bty = "n")
    
    par(mfrow = c(1, 1))
    
  } else {
    cat("Impossibile creare grafici: dati non validi\n")
  }
  
  # -------------------------------------------------------------
  # Statistiche riassuntive
  # -------------------------------------------------------------
  
  cat("\n", strrep("=", 70), "\n")
  cat("STATISTICHE RIASSUNTIVE: VARIAZIONE TEMPORALE\n")
  cat(strrep("=", 70), "\n\n")
  
  tau_values <- as.numeric(sapply(results_clayton, function(x) x$tau_mean))
  lambda_L_values <- as.numeric(sapply(results_clayton, function(x) x$lambda_L_mean))
  
  cat("Tau di Kendall:\n")
  cat(sprintf("  Min:    %.3f (%s)\n", 
              min(tau_values), 
              results_clayton[[which.min(tau_values)]]$period_name))
  cat(sprintf("  Max:    %.3f (%s)\n", 
              max(tau_values), 
              results_clayton[[which.max(tau_values)]]$period_name))
  cat(sprintf("  Range:  %.3f\n", max(tau_values) - min(tau_values)))
  cat(sprintf("  CV:     %.2f%%\n", 100 * sd(tau_values) / mean(tau_values)))
  
  cat("\nLambda_L (dipendenza coda inferiore):\n")
  cat(sprintf("  Min:    %.3f (%s)\n", 
              min(lambda_L_values), 
              results_clayton[[which.min(lambda_L_values)]]$period_name))
  cat(sprintf("  Max:    %.3f (%s)\n", 
              max(lambda_L_values), 
              results_clayton[[which.max(lambda_L_values)]]$period_name))
  cat(sprintf("  Range:  %.3f\n", max(lambda_L_values) - min(lambda_L_values)))
  cat(sprintf("  CV:     %.2f%%\n", 100 * sd(lambda_L_values) / mean(lambda_L_values)))
  
  # Confronto crisi vs. non-crisi
  if (length(results_clayton) >= 4) {
    crisis_idx <- c(2, 5)  # GFC e COVID
    available_crisis <- crisis_idx[crisis_idx <= length(results_clayton)]
    normal_idx <- setdiff(1:length(results_clayton), available_crisis)
    
    if (length(available_crisis) > 0 && length(normal_idx) > 0) {
      tau_crisis <- mean(sapply(results_clayton[available_crisis], function(x) x$tau_mean))
      tau_normal <- mean(sapply(results_clayton[normal_idx], function(x) x$tau_mean))
      
      lambda_crisis <- mean(sapply(results_clayton[available_crisis], function(x) x$lambda_L_mean))
      lambda_normal <- mean(sapply(results_clayton[normal_idx], function(x) x$lambda_L_mean))
      
      cat("\nConfronto periodi di crisi vs. periodi normali:\n")
      cat(sprintf("  Tau medio (crisi):     %.3f\n", tau_crisis))
      cat(sprintf("  Tau medio (normale):   %.3f\n", tau_normal))
      cat(sprintf("  Incremento relativo:   %.1f%%\n", 
                  100 * (tau_crisis - tau_normal) / tau_normal))
      
      cat(sprintf("\n  Lambda_L medio (crisi):   %.3f\n", lambda_crisis))
      cat(sprintf("  Lambda_L medio (normale): %.3f\n", lambda_normal))
      cat(sprintf("  Incremento relativo:      %.1f%%\n", 
                  100 * (lambda_crisis - lambda_normal) / lambda_normal))
    }
  }
}

cat("\n", strrep("=", 70), "\n")
cat("ANALISI PER SOTTO-PERIODI COMPLETATA\n")
cat(strrep("=", 70), "\n\n")


###############################################################################
# SEZIONE 3.6 BIS: ANALISI PER SOTTO-PERIODI - GUMBEL (TUTTI I PERIODI)
###############################################################################

cat("\n", strrep("-", 70), "\n")
cat("STIMA COPULA DI GUMBEL PER SOTTO-PERIODI\n")
cat(strrep("-", 70), "\n")

# -------------------------------------------------------------
# Stima Gumbel per tutti i sotto-periodi (stessa identica logica di Clayton)
# -------------------------------------------------------------

set.seed(321)
results_gumbel_all <- list()

for (i in seq_along(periods)) {
  res <- estimate_subperiod(
    data_full = data,
    period_info = periods[[i]],
    family = "gumbel",
    n_iter = 30000,
    burn_in = 3000,
    sigma = 0.6,
    seed = 300 + i
  )
  
  if (!is.null(res)) {
    results_gumbel_all[[i]] <- res
  }
}

# Rimuovi eventuali NULL
results_gumbel_all <- results_gumbel_all[!sapply(results_gumbel_all, is.null)]

cat(sprintf("\n\nPeriodi analizzati con successo (Gumbel): %d\n", length(results_gumbel_all)))

if (length(results_gumbel_all) == 0) {
  cat("ERRORE: Nessun sotto-periodo Gumbel analizzato con successo!\n")
} else {
  
  # -------------------------------------------------------------
  # Tabella riassuntiva Gumbel (tutti i periodi)
  # -------------------------------------------------------------
  
  cat("\n", strrep("=", 70), "\n")
  cat("TABELLA RIASSUNTIVA: COPULA DI GUMBEL (SOTTO-PERIODI)\n")
  cat(strrep("=", 70), "\n\n")
  
  table_gumbel_all <- data.frame(
    Periodo = sapply(results_gumbel_all, function(x) x$period_name),
    n = as.numeric(sapply(results_gumbel_all, function(x) x$n_obs)),
    tau_hat = as.numeric(sapply(results_gumbel_all, function(x) x$tau_mean)),
    tau_ci_low = as.numeric(sapply(results_gumbel_all, function(x) x$tau_ci_low)),
    tau_ci_high = as.numeric(sapply(results_gumbel_all, function(x) x$tau_ci_high)),
    lambda_U_hat = as.numeric(sapply(results_gumbel_all, function(x) x$lambda_U_mean)),
    lambda_U_ci_low = as.numeric(sapply(results_gumbel_all, function(x) x$lambda_U_ci_low)),
    lambda_U_ci_high = as.numeric(sapply(results_gumbel_all, function(x) x$lambda_U_ci_high)),
    stringsAsFactors = FALSE
  )
  
  table_gumbel_all$tau_IC <- sprintf("[%.3f, %.3f]",
                                     table_gumbel_all$tau_ci_low,
                                     table_gumbel_all$tau_ci_high)
  
  table_gumbel_all$lambda_U_IC <- sprintf("[%.3f, %.3f]",
                                          table_gumbel_all$lambda_U_ci_low,
                                          table_gumbel_all$lambda_U_ci_high)
  
  table_gumbel_all_final <- data.frame(
    Periodo = table_gumbel_all$Periodo,
    n = table_gumbel_all$n,
    tau_hat = round(table_gumbel_all$tau_hat, 3),
    `IC_95_tau` = table_gumbel_all$tau_IC,
    lambda_U_hat = round(table_gumbel_all$lambda_U_hat, 3),
    `IC_95_lambda_U` = table_gumbel_all$lambda_U_IC,
    stringsAsFactors = FALSE
  )
  
  print(table_gumbel_all_final, row.names = FALSE)
  
  cat("\n\nFormato LaTeX:\n")
  for (i in 1:nrow(table_gumbel_all_final)) {
    cat(sprintf("%s & %d & %.3f & %s & %.3f & %s \\\\\n",
                table_gumbel_all_final$Periodo[i],
                table_gumbel_all_final$n[i],
                table_gumbel_all_final$tau_hat[i],
                table_gumbel_all_final$IC_95_tau[i],
                table_gumbel_all_final$lambda_U_hat[i],
                table_gumbel_all_final$IC_95_lambda_U[i]))
  }
  
  # -------------------------------------------------------------
  # Visualizzazione: Evoluzione temporale (Gumbel)
  # -------------------------------------------------------------
  
  cat("\n", strrep("-", 70), "\n")
  cat("VISUALIZZAZIONE: Evoluzione temporale (Gumbel)\n")
  cat(strrep("-", 70), "\n\n")
  
  plot_data_gumbel <- data.frame(
    periodo = sapply(results_gumbel_all, function(x) x$period_name),
    start_date = as.Date(sapply(results_gumbel_all, function(x) as.character(x$start_date))),
    end_date = as.Date(sapply(results_gumbel_all, function(x) as.character(x$end_date))),
    tau = as.numeric(sapply(results_gumbel_all, function(x) x$tau_mean)),
    tau_low = as.numeric(sapply(results_gumbel_all, function(x) x$tau_ci_low)),
    tau_high = as.numeric(sapply(results_gumbel_all, function(x) x$tau_ci_high)),
    lambda_U = as.numeric(sapply(results_gumbel_all, function(x) x$lambda_U_mean)),
    lambda_U_low = as.numeric(sapply(results_gumbel_all, function(x) x$lambda_U_ci_low)),
    lambda_U_high = as.numeric(sapply(results_gumbel_all, function(x) x$lambda_U_ci_high)),
    stringsAsFactors = FALSE
  )
  
  plot_data_gumbel$midpoint <- plot_data_gumbel$start_date +
    as.numeric(difftime(plot_data_gumbel$end_date, plot_data_gumbel$start_date, units = "days")) / 2
  
  if (nrow(plot_data_gumbel) > 0 && all(is.finite(plot_data_gumbel$tau))) {
    
    par(mfrow = c(2, 1), mar = c(5, 4, 4, 2))
    
    # Grafico 1: tau
    plot(plot_data_gumbel$midpoint, plot_data_gumbel$tau,
         type = "b", pch = 19, col = "darkgreen", lwd = 2, cex = 1.5,
         ylim = range(c(plot_data_gumbel$tau_low, plot_data_gumbel$tau_high)),
         xlab = "Periodo",
         ylab = expression(tau~"(Kendall)"),
         main = "Evoluzione temporale della dipendenza (Gumbel)")
    
    arrows(plot_data_gumbel$midpoint, plot_data_gumbel$tau_low,
           plot_data_gumbel$midpoint, plot_data_gumbel$tau_high,
           angle = 90, code = 3, length = 0.05, col = "darkgreen", lwd = 1.5)
    
    # Media intero periodo (fit_gumbel dall'intero campione, già stimato in Sez. 3.4)
    abline(h = mean(fit_gumbel$tau_samples), col = "red", lty = 2, lwd = 2)
    
    grid()
    
    legend("topleft",
           legend = c("Stima per sotto-periodo", "IC 95%", "Media intero periodo"),
           col = c("darkgreen", "darkgreen", "red"),
           lty = c(1, 1, 2),
           pch = c(19, NA, NA),
           lwd = c(2, 1.5, 2),
           bty = "n")
    
    # Grafico 2: lambda_U
    plot(plot_data_gumbel$midpoint, plot_data_gumbel$lambda_U,
         type = "b", pch = 19, col = "darkorange", lwd = 2, cex = 1.5,
         ylim = range(c(plot_data_gumbel$lambda_U_low, plot_data_gumbel$lambda_U_high)),
         xlab = "Periodo",
         ylab = expression(lambda[U]),
         main = "Evoluzione della dipendenza nella coda superiore (Gumbel)")
    
    arrows(plot_data_gumbel$midpoint, plot_data_gumbel$lambda_U_low,
           plot_data_gumbel$midpoint, plot_data_gumbel$lambda_U_high,
           angle = 90, code = 3, length = 0.05, col = "darkorange", lwd = 1.5)
    
    # Media intero periodo per lambda_U
    td_full_g <- tail_dependence(fit_gumbel$theta_samples, "gumbel")
    abline(h = mean(td_full_g$lambda_U), col = "red", lty = 2, lwd = 2)
    
    grid()
    
    legend("topleft",
           legend = c("Stima per sotto-periodo", "IC 95%", "Media intero periodo"),
           col = c("darkorange", "darkorange", "red"),
           lty = c(1, 1, 2),
           pch = c(19, NA, NA),
           lwd = c(2, 1.5, 2),
           bty = "n")
    
    par(mfrow = c(1, 1))
    
  } else {
    cat("Impossibile creare grafici Gumbel: dati non validi\n")
  }
  
  # -------------------------------------------------------------
  # Statistiche riassuntive (Gumbel)
  # -------------------------------------------------------------
  
  cat("\n", strrep("=", 70), "\n")
  cat("STATISTICHE RIASSUNTIVE: VARIAZIONE TEMPORALE (Gumbel)\n")
  cat(strrep("=", 70), "\n\n")
  
  tau_values_g <- as.numeric(sapply(results_gumbel_all, function(x) x$tau_mean))
  lambda_U_values_g <- as.numeric(sapply(results_gumbel_all, function(x) x$lambda_U_mean))
  
  cat("Tau di Kendall (Gumbel):\n")
  cat(sprintf("  Min:    %.3f (%s)\n",
              min(tau_values_g),
              results_gumbel_all[[which.min(tau_values_g)]]$period_name))
  cat(sprintf("  Max:    %.3f (%s)\n",
              max(tau_values_g),
              results_gumbel_all[[which.max(tau_values_g)]]$period_name))
  cat(sprintf("  Range:  %.3f\n", max(tau_values_g) - min(tau_values_g)))
  cat(sprintf("  CV:     %.2f%%\n", 100 * sd(tau_values_g) / mean(tau_values_g)))
  
  cat("\nLambda_U (dipendenza coda superiore, Gumbel):\n")
  cat(sprintf("  Min:    %.3f (%s)\n",
              min(lambda_U_values_g),
              results_gumbel_all[[which.min(lambda_U_values_g)]]$period_name))
  cat(sprintf("  Max:    %.3f (%s)\n",
              max(lambda_U_values_g),
              results_gumbel_all[[which.max(lambda_U_values_g)]]$period_name))
  cat(sprintf("  Range:  %.3f\n", max(lambda_U_values_g) - min(lambda_U_values_g)))
  cat(sprintf("  CV:     %.2f%%\n", 100 * sd(lambda_U_values_g) / mean(lambda_U_values_g)))
}

cat("\n", strrep("-", 70), "\n")
cat("ANALISI GUMBEL PER SOTTO-PERIODI COMPLETATA\n")
cat(strrep("-", 70), "\n\n")
