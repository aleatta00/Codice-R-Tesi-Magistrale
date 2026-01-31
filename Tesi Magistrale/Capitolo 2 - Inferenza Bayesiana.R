###############################################################################
# CODICE PER CAPITOLO 2: INFERENZA BAYESIANA SU COPULE ARCHIMEDEE
#
# Il presente codice implementa i metodi di inferenza bayesiana discussi
# nel Capitolo 2, con particolare attenzione alla parametrizzazione basata
# sul coefficiente di Kendall tau.
#
# Contenuto:
#   - Funzioni ausiliarie per trasformazioni e calcoli
#   - Verifica equivalenza parametrizzazioni (Sez. 2.6.1)
#   - Confronto efficienza MCMC (Sez. 2.6.2)
#   - Analisi sensibilità prior (Sez. 2.6.3)
#   - Studio Monte Carlo (Sez. 2.7)
#   - Confronto Frequentista vs Bayesiano
###############################################################################

# -------------------------------------------------------------
# Librerie
# -------------------------------------------------------------
suppressPackageStartupMessages({
  library(copula)
  library(ggplot2)
  library(gridExtra)
  library(coda)
  library(parallel)
  library(knitr)
  library(dplyr)
})

set.seed(42)
options(digits = 4)

###############################################################################
# FUNZIONI AUSILIARIE
###############################################################################

# -------------------------------------------------------------
# Trasformazioni Clayton: tau <-> theta
# -------------------------------------------------------------

tau_to_theta_clayton <- function(tau) {
  2 * tau / (1 - tau)
}

theta_to_tau_clayton <- function(theta) {
  theta / (theta + 2)
}

# Trasformazioni per spazio non vincolato
tau_from_z <- function(z) tanh(z)
z_from_tau <- function(tau) atanh(tau)

# -------------------------------------------------------------
# Log-verosimiglianza copula di Clayton
# -------------------------------------------------------------

loglik_clayton <- function(theta, u, v) {
  if (theta <= 0) return(-Inf)
  cop <- claytonCopula(theta, dim = 2)
  tryCatch(
    sum(log(dCopula(cbind(u, v), cop))),
    error = function(e) -Inf
  )
}

# Versione alternativa con calcolo diretto (più stabile)
log_clayton_density <- function(u, v, theta) {
  if (!is.finite(theta) || theta <= 0) return(-Inf)
  
  term1 <- log1p(theta)
  term2 <- -(1 + theta) * (log(u) + log(v))
  
  lu <- -theta * log(u)
  lv <- -theta * log(v)
  m  <- max(lu, lv)
  s  <- exp(lu - m) + exp(lv - m) - exp(-m)
  
  if (!is.finite(s) || s <= 0) return(-Inf)
  
  term1 + term2 - (2 + 1/theta) * (log(s) + m)
}

log_lik_tau <- function(tau, u, v) {
  if (!is.finite(tau) || tau <= -1 || tau >= 1) return(-Inf)
  theta <- tau_to_theta_clayton(tau)
  sum(mapply(log_clayton_density, u, v, MoreArgs = list(theta = theta)))
}

# -------------------------------------------------------------
# Distribuzioni a priori
# -------------------------------------------------------------

# Prior uniforme su tau in (-1, 1)
log_prior_uniform_tau <- function(tau) {
  if (tau <= -1 || tau >= 1) return(-Inf)
  -log(2)
}

# Prior Beta(a, b) riscalata su (-1, 1)
log_prior_beta_tau <- function(tau, a = 2, b = 2) {
  if (tau <= -1 || tau >= 1) return(-Inf)
  t <- (tau + 1) / 2
  dbeta(t, a, b, log = TRUE) - log(2)
}

# Prior indotta su theta da Unif(tau): pi(theta) = 1/(theta+2)^2
log_prior_induced_theta <- function(theta) {
  if (theta <= 0) return(-Inf)
  -2 * log(theta + 2)
}

# Prior di Jeffreys su theta: pi(theta) propto 1/theta
log_prior_jeffreys_theta <- function(theta) {
  if (theta <= 0) return(-Inf)
  -log(theta)
}

# Prior di Jeffreys indotta su tau (Clayton)
log_prior_jeffreys_induced_tau <- function(tau, eps = 0.01) {
  if (tau <= eps || tau >= 1 - eps) return(-Inf)
  -log(tau) - log(1 - tau)
}

# -------------------------------------------------------------
# Funzioni per calcoli numerici su griglia
# -------------------------------------------------------------

# Normalizza densità su griglia uniforme
normalize_on_grid <- function(x, dt) {
  x / (sum(x) * dt)
}

# Calcola quantili da griglia discreta
quantile_from_grid <- function(grid, density, probs) {
  cum_dens <- cumsum(density * c(diff(grid)[1], diff(grid)))
  cum_dens <- cum_dens / max(cum_dens)
  sapply(probs, function(p) grid[which.min(abs(cum_dens - p))])
}

# -------------------------------------------------------------
# Algoritmo MCMC: Random Walk Metropolis-Hastings
# -------------------------------------------------------------

mcmc_rw_mh <- function(log_target, init, n_iter = 50000, n_burn = 5000,
                       sigma_prop = 0.5, adapt = TRUE, target_acc = 0.234) {
  
  n_total <- n_iter + n_burn
  chain <- numeric(n_total)
  chain[1] <- init
  n_accept <- 0
  
  for (t in 2:n_total) {
    proposal <- rnorm(1, chain[t-1], sigma_prop)
    
    log_alpha <- log_target(proposal) - log_target(chain[t-1])
    
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      chain[t] <- proposal
      n_accept <- n_accept + 1
    } else {
      chain[t] <- chain[t-1]
    }
    
    # Adattamento sigma durante primi 3000 iter di burn-in
    if (adapt && t <= min(3000, n_burn)) {
      acc_rate_current <- n_accept / t
      if (acc_rate_current > target_acc) {
        sigma_prop <- sigma_prop * 1.05
      } else {
        sigma_prop <- sigma_prop * 0.95
      }
      sigma_prop <- max(0.01, min(sigma_prop, 5.0))
    }
  }
  
  list(
    chain_full = chain,
    samples = chain[(n_burn + 1):n_total],
    burn_in = chain[1:n_burn],
    acc_rate = n_accept / n_total,
    sigma_final = sigma_prop
  )
}

# -------------------------------------------------------------
# Wrapper MCMC per parametrizzazione TAU (via z = atanh)
# -------------------------------------------------------------

mcmc_tau_z <- function(u, v, n_iter = 50000, n_burn = 5000, 
                       sigma_init = 0.6, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Inizializzazione da Kendall empirico
  tau_init <- cor(u, v, method = "kendall")
  tau_init <- max(min(tau_init, 0.95), 0.05)
  z_init <- z_from_tau(tau_init)
  
  # Log-target in spazio z
  log_target_z <- function(z) {
    tau <- tau_from_z(z)
    if (tau <= -1 || tau >= 1) return(-Inf)
    
    theta <- tau_to_theta_clayton(tau)
    if (theta <= 0) return(-Inf)
    
    ll <- loglik_clayton(theta, u, v)
    lp <- log_prior_uniform_tau(tau)
    jacobian <- log(1 - tau^2)
    
    ll + lp + jacobian
  }
  
  t_start <- Sys.time()
  fit <- mcmc_rw_mh(log_target_z, z_init, n_iter, n_burn, sigma_init)
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  
  tau_samples <- tau_from_z(fit$samples)
  
  tau_mcmc <- mcmc(tau_samples)
  ess <- effectiveSize(tau_mcmc)
  
  list(
    samples_tau = tau_samples,
    chain_full_tau = tau_from_z(fit$chain_full),
    burn_in_tau = tau_from_z(fit$burn_in),
    ess = as.numeric(ess),
    acc_rate = fit$acc_rate,
    elapsed = elapsed,
    ess_per_sec = as.numeric(ess) / elapsed
  )
}

# -------------------------------------------------------------
# Wrapper MCMC per parametrizzazione THETA (via w = log)
# -------------------------------------------------------------

mcmc_theta_w <- function(u, v, n_iter = 50000, n_burn = 5000,
                         sigma_init = 0.4, use_induced_prior = TRUE,
                         seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  tau_init <- cor(u, v, method = "kendall")
  tau_init <- max(min(tau_init, 0.95), 0.05)
  theta_init <- tau_to_theta_clayton(tau_init)
  w_init <- log(theta_init)
  
  log_target_w <- function(w) {
    theta <- exp(w)
    if (theta <= 0) return(-Inf)
    
    ll <- loglik_clayton(theta, u, v)
    
    if (use_induced_prior) {
      lp <- log_prior_induced_theta(theta)
    } else {
      lp <- log_prior_jeffreys_theta(theta)
    }
    
    jacobian <- w
    
    ll + lp + jacobian
  }
  
  t_start <- Sys.time()
  fit <- mcmc_rw_mh(log_target_w, w_init, n_iter, n_burn, sigma_init)
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  
  theta_samples <- exp(fit$samples)
  tau_samples <- theta_to_tau_clayton(theta_samples)
  
  tau_mcmc <- mcmc(tau_samples)
  ess <- effectiveSize(tau_mcmc)
  
  list(
    samples_tau = tau_samples,
    samples_theta = theta_samples,
    ess = as.numeric(ess),
    acc_rate = fit$acc_rate,
    elapsed = elapsed,
    ess_per_sec = as.numeric(ess) / elapsed
  )
}


###############################################################################
# SEZIONE 2.6.1: VERIFICA EQUIVALENZA PARAMETRIZZAZIONI
###############################################################################

# -------------------------------------------------------------
# Parametri di simulazione
# -------------------------------------------------------------
n <- 200
tau_true <- 0.5
theta_true <- tau_to_theta_clayton(tau_true)

# -------------------------------------------------------------
# Generazione dati dalla copula di Clayton
# -------------------------------------------------------------
clay_cop <- claytonCopula(theta_true, dim = 2)
data_sim <- rCopula(n, clay_cop)
u <- data_sim[, 1]
v <- data_sim[, 2]

cat(sprintf("Dataset: n = %d, tau_true = %.2f, theta_true = %.2f\n\n", 
            n, tau_true, theta_true))

# -------------------------------------------------------------
# Calcolo posteriori: parametrizzazione diretta su TAU
# -------------------------------------------------------------
tau_grid <- seq(0.01, 0.95, length.out = 500)

log_post_tau_direct <- sapply(tau_grid, function(tau) {
  theta <- tau_to_theta_clayton(tau)
  if (theta <= 0) return(-Inf)
  ll <- loglik_clayton(theta, u, v)
  lp <- log_prior_uniform_tau(tau)
  ll + lp
})

post_tau_direct <- exp(log_post_tau_direct - max(log_post_tau_direct))
post_tau_direct <- post_tau_direct / (sum(post_tau_direct) * diff(tau_grid)[1])

mean_tau_direct <- sum(tau_grid * post_tau_direct) * diff(tau_grid)[1]
sd_tau_direct <- sqrt(sum((tau_grid - mean_tau_direct)^2 * post_tau_direct) * 
                        diff(tau_grid)[1])
ic_tau_direct <- quantile_from_grid(tau_grid, post_tau_direct, 
                                    c(0.025, 0.5, 0.975))

# -------------------------------------------------------------
# Calcolo posteriori: parametrizzazione su THETA (prior indotta)
# -------------------------------------------------------------
theta_grid <- seq(0.05, 10, length.out = 500)

log_post_theta_induced <- sapply(theta_grid, function(theta) {
  ll <- loglik_clayton(theta, u, v)
  lp <- log_prior_induced_theta(theta)
  ll + lp
})

post_theta_induced <- exp(log_post_theta_induced - max(log_post_theta_induced))
post_theta_induced <- post_theta_induced / (sum(post_theta_induced) * 
                                              diff(theta_grid)[1])

tau_from_theta <- theta_to_tau_clayton(theta_grid)

mean_tau_from_theta <- sum(tau_from_theta * post_theta_induced) * 
  diff(theta_grid)[1]
sd_tau_from_theta <- sqrt(sum((tau_from_theta - mean_tau_from_theta)^2 * 
                                post_theta_induced) * diff(theta_grid)[1])
ic_tau_from_theta <- quantile_from_grid(tau_from_theta, post_theta_induced,
                                        c(0.025, 0.5, 0.975))

# -------------------------------------------------------------
# Statistiche riassuntive
# -------------------------------------------------------------
cat("\n=== Confronto statistiche posteriori ===\n\n")

cat(sprintf("%-25s %15s %15s\n", "Statistica", "Param. τ", "Param. θ"))
cat(strrep("-", 60), "\n")
cat(sprintf("%-25s %15.4f %15.4f\n", "Media", 
            mean_tau_direct, mean_tau_from_theta))
cat(sprintf("%-25s %15.4f %15.4f\n", "Mediana", 
            ic_tau_direct[2], ic_tau_from_theta[2]))
cat(sprintf("%-25s %15.4f %15.4f\n", "Dev. standard", 
            sd_tau_direct, sd_tau_from_theta))
cat(sprintf("%-25s [%.3f, %.3f] [%.3f, %.3f]\n", "IC 95%", 
            ic_tau_direct[1], ic_tau_direct[3],
            ic_tau_from_theta[1], ic_tau_from_theta[3]))
cat(strrep("-", 60), "\n")
cat(sprintf("%-25s %15.6f\n", "Differenza media", 
            abs(mean_tau_direct - mean_tau_from_theta)))

# Test Kolmogorov-Smirnov
sample_tau_direct <- sample(tau_grid, 10000, replace = TRUE, 
                            prob = post_tau_direct)
sample_tau_from_theta <- sample(tau_from_theta, 10000, replace = TRUE,
                                prob = post_theta_induced)
ks_test <- ks.test(sample_tau_direct, sample_tau_from_theta)
cat(sprintf("%-25s %15.4f\n", "KS test p-value", ks_test$p.value))
cat("(p > 0.05 indica equivalenza)\n\n")

# -------------------------------------------------------------
# Visualizzazione: overlay posteriori
# -------------------------------------------------------------
df_overlay <- data.frame(
  tau = c(tau_grid, tau_from_theta),
  density = c(post_tau_direct, post_theta_induced),
  Parametrizzazione = rep(c("Diretta su τ", "Indotta da θ"), 
                          c(length(tau_grid), length(theta_grid)))
)

p_overlay <- ggplot(df_overlay, aes(x = tau, y = density, 
                                    color = Parametrizzazione,
                                    linetype = Parametrizzazione)) +
  geom_line(linewidth = 1.3) +
  geom_vline(xintercept = tau_true, linetype = "dashed", 
             color = "black", linewidth = 0.8, alpha = 0.6) +
  annotate("text", x = tau_true + 0.08, y = max(df_overlay$density) * 0.95,
           label = paste0("tau[0] == ", tau_true), parse = TRUE, size = 5) +
  scale_color_manual(values = c("darkgreen", "blue")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(
    title = "Equivalenza posteriori con prior coerenti",
    x = expression(tau),
    y = "Densità a posteriori"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

print(p_overlay)


###############################################################################
# SEZIONE 2.6.2: CONFRONTO EFFICIENZA MCMC
###############################################################################

# -------------------------------------------------------------
# Parametri MCMC
# -------------------------------------------------------------
n_iter_comp <- 50000
n_burn_comp <- 5000

# -------------------------------------------------------------
# Esecuzione MCMC: parametrizzazione tau
# -------------------------------------------------------------
cat("\nEsecuzione MCMC per confronto efficienza...\n")
cat("(Richiede alcuni minuti)\n\n")

fit_tau <- mcmc_tau_z(u, v, n_iter_comp, n_burn_comp, 
                      sigma_init = 0.6, seed = 100)

# -------------------------------------------------------------
# Esecuzione MCMC: parametrizzazione theta (prior indotta)
# -------------------------------------------------------------
fit_theta_ind <- mcmc_theta_w(u, v, n_iter_comp, n_burn_comp,
                              sigma_init = 0.4, use_induced_prior = TRUE,
                              seed = 101)

# -------------------------------------------------------------
# Calcolo autocorrelazione
# -------------------------------------------------------------
acf_tau <- acf(fit_tau$samples_tau, plot = FALSE, lag.max = 50)
acf_theta <- acf(fit_theta_ind$samples_tau, plot = FALSE, lag.max = 50)

# -------------------------------------------------------------
# Statistiche di efficienza
# -------------------------------------------------------------
cat("\n=== Metriche di efficienza MCMC ===\n\n")

cat(sprintf("%-20s %15s %15s\n", "Metrica", "Param. τ", "Param. θ"))
cat(strrep("-", 55), "\n")
cat(sprintf("%-20s %15.0f %15.0f\n", "ESS", 
            fit_tau$ess, fit_theta_ind$ess))
cat(sprintf("%-20s %15.2f %15.2f\n", "ESS/sec", 
            fit_tau$ess_per_sec, fit_theta_ind$ess_per_sec))
cat(sprintf("%-20s %15.3f %15.3f\n", "Acceptance rate", 
            fit_tau$acc_rate, fit_theta_ind$acc_rate))
cat(sprintf("%-20s %15.4f %15.4f\n", "ACF(lag=50)", 
            acf_tau$acf[51], acf_theta$acf[51]))
cat(sprintf("%-20s %15.1f %15.1f\n", "Tempo (sec)", 
            fit_tau$elapsed, fit_theta_ind$elapsed))
cat(strrep("-", 55), "\n\n")

# -------------------------------------------------------------
# Visualizzazione: traceplot
# -------------------------------------------------------------
df_trace <- data.frame(
  iteration = rep(1:n_iter_comp, 2),
  tau = c(fit_tau$samples_tau, fit_theta_ind$samples_tau),
  Parametrizzazione = rep(c("τ (via z = atanh)", "θ (via w = log)"),
                          each = n_iter_comp)
)

p_trace <- ggplot(df_trace, aes(x = iteration, y = tau)) +
  geom_line(alpha = 0.6, linewidth = 0.3) +
  geom_hline(yintercept = tau_true, color = "red", 
             linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ Parametrizzazione, ncol = 1) +
  labs(
    title = "Traceplot: confronto mixing",
    x = "Iterazione (post burn-in)",
    y = expression(tau)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA)
  )

print(p_trace)

# -------------------------------------------------------------
# Visualizzazione: ACF
# -------------------------------------------------------------
df_acf <- data.frame(
  lag = rep(acf_tau$lag, 2),
  acf = c(acf_tau$acf, acf_theta$acf),
  Parametrizzazione = rep(c("τ (via z)", "θ (via w)"),
                          each = length(acf_tau$lag))
)

p_acf <- ggplot(df_acf, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_hline(yintercept = c(-1.96, 1.96) / sqrt(n_iter_comp),
             color = "blue", linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 1.2) +
  facet_wrap(~ Parametrizzazione, ncol = 2) +
  labs(
    title = "Autocorrelazione",
    subtitle = "Linee blu = IC 95% per rumore bianco",
    x = "Lag",
    y = "ACF"
  ) +
  ylim(-0.1, 1) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA)
  )

print(p_acf)


###############################################################################
# SEZIONE 2.6.3: ANALISI SENSIBILITÀ PRIOR
###############################################################################

# -------------------------------------------------------------
# Calcolo log-likelihood (comune a tutti i modelli)
# -------------------------------------------------------------
tau_grid_sens <- seq(0.01, 0.95, length.out = 600)
dt_sens <- diff(tau_grid_sens)[1]

loglik_tau_sens <- sapply(tau_grid_sens, function(tau) {
  theta <- tau_to_theta_clayton(tau)
  if (theta <= 0) return(-Inf)
  loglik_clayton(theta, u, v)
})

# -------------------------------------------------------------
# Prior 1: Uniforme su tau
# -------------------------------------------------------------
logprior_unif <- sapply(tau_grid_sens, log_prior_uniform_tau)
post_unif <- exp((loglik_tau_sens + logprior_unif) - 
                   max(loglik_tau_sens + logprior_unif))
post_unif <- normalize_on_grid(post_unif, dt_sens)

mean_unif <- sum(tau_grid_sens * post_unif) * dt_sens
sd_unif <- sqrt(sum((tau_grid_sens - mean_unif)^2 * post_unif) * dt_sens)
ic_unif <- quantile_from_grid(tau_grid_sens, post_unif, c(0.025, 0.975))

# -------------------------------------------------------------
# Prior 2: Beta(2,2) riscalata
# -------------------------------------------------------------
logprior_beta <- sapply(tau_grid_sens, function(t) 
  log_prior_beta_tau(t, 2, 2))
post_beta <- exp((loglik_tau_sens + logprior_beta) - 
                   max(loglik_tau_sens + logprior_beta))
post_beta <- normalize_on_grid(post_beta, dt_sens)

mean_beta <- sum(tau_grid_sens * post_beta) * dt_sens
sd_beta <- sqrt(sum((tau_grid_sens - mean_beta)^2 * post_beta) * dt_sens)
ic_beta <- quantile_from_grid(tau_grid_sens, post_beta, c(0.025, 0.975))

# -------------------------------------------------------------
# Prior 3: Jeffreys su theta (indotta su tau)
# -------------------------------------------------------------
logprior_jeff_tau <- sapply(tau_grid_sens, function(t) 
  log_prior_jeffreys_induced_tau(t, eps = 0.01))
post_jeff <- exp((loglik_tau_sens + logprior_jeff_tau) - 
                   max(loglik_tau_sens + logprior_jeff_tau))
post_jeff <- normalize_on_grid(post_jeff, dt_sens)

mean_jeff <- sum(tau_grid_sens * post_jeff) * dt_sens
sd_jeff <- sqrt(sum((tau_grid_sens - mean_jeff)^2 * post_jeff) * dt_sens)
ic_jeff <- quantile_from_grid(tau_grid_sens, post_jeff, c(0.025, 0.975))

# -------------------------------------------------------------
# Statistiche riassuntive
# -------------------------------------------------------------
cat("\n=== Sensibilità alle prior ===\n\n")

cat(sprintf("%-20s %12s %12s %12s\n", "Prior", "Mean", "SD", "IC 95%"))
cat(strrep("-", 60), "\n")
cat(sprintf("%-20s %12.3f %12.3f [%.3f, %.3f]\n", 
            "Unif(τ)", mean_unif, sd_unif, ic_unif[1], ic_unif[2]))
cat(sprintf("%-20s %12.3f %12.3f [%.3f, %.3f]\n", 
            "Beta(2,2)", mean_beta, sd_beta, ic_beta[1], ic_beta[2]))
cat(sprintf("%-20s %12.3f %12.3f [%.3f, %.3f]\n", 
            "Jeffreys(θ)", mean_jeff, sd_jeff, ic_jeff[1], ic_jeff[2]))
cat(strrep("-", 60), "\n")

diff_max <- max(abs(c(mean_unif - mean_jeff, mean_unif - mean_beta,
                      mean_jeff - mean_beta)))
cat(sprintf("Differenza massima tra medie: %.4f\n", diff_max))
cat("(Impatto prior limitato con n=200)\n\n")

# -------------------------------------------------------------
# Visualizzazione: prior separate
# -------------------------------------------------------------
prior_unif <- rep(0.5, length(tau_grid_sens))
prior_beta <- exp(logprior_beta)
prior_beta <- normalize_on_grid(prior_beta, dt_sens)
prior_jeff <- exp(logprior_jeff_tau)
prior_jeff <- normalize_on_grid(prior_jeff, dt_sens)

df_prior <- rbind(
  data.frame(tau = tau_grid_sens, density = prior_unif, model = "Unif(τ)"),
  data.frame(tau = tau_grid_sens, density = prior_beta, model = "Beta(2,2)"),
  data.frame(tau = tau_grid_sens, density = prior_jeff, model = "Jeffreys(θ)")
)

p_prior <- ggplot(df_prior, aes(x = tau, y = density, color = model)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = tau_true, linetype = "dotted",
             color = "black", linewidth = 0.7) +
  labs(
    title = "Confronto tra prior (scala τ)",
    x = expression(tau),
    y = "Densità"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

print(p_prior)

# -------------------------------------------------------------
# Visualizzazione: posteriori separate
# -------------------------------------------------------------
df_post <- rbind(
  data.frame(tau = tau_grid_sens, density = post_unif, model = "Unif(τ)"),
  data.frame(tau = tau_grid_sens, density = post_beta, model = "Beta(2,2)"),
  data.frame(tau = tau_grid_sens, density = post_jeff, model = "Jeffreys(θ)")
)

p_post <- ggplot(df_post, aes(x = tau, y = density, color = model)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = tau_true, linetype = "dotted",
             color = "black", linewidth = 0.7) +
  labs(
    title = "Confronto tra posteriori (scala τ)",
    subtitle = paste0("Clayton: tau0 = ", tau_true, ", n = ", n),
    x = expression(tau),
    y = "Densità a posteriori"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

print(p_post)


###############################################################################
# SEZIONE 2.7: STUDIO MONTE CARLO
###############################################################################

# -------------------------------------------------------------
# Worker per singola replica Monte Carlo
# -------------------------------------------------------------

mc_worker <- function(tau_true, n, n_iter_mc = 10000, n_burn_mc = 2000,
                      seed_offset = 0) {
  set.seed(seed_offset)
  theta_true <- tau_to_theta_clayton(tau_true)
# Simula dati
  data <- rCopula(n, claytonCopula(theta_true, dim = 2))
  u <- data[, 1]
  v <- data[, 2]
# Esegui MCMC
  fit <- mcmc_tau_z(u, v, n_iter_mc, n_burn_mc, sigma_init = 0.6)
# Statistiche posteriori
  tau_mean <- mean(fit$samples_tau)
  tau_hpd <- HPDinterval(mcmc(fit$samples_tau), prob = 0.95)
  coverage <- (tau_true >= tau_hpd[1]) && (tau_true <= tau_hpd[2])
  c(
    mean = tau_mean,
    coverage = as.numeric(coverage),
    ess = fit$ess,
    acc_rate = fit$acc_rate
  )
}
#-------------------------------------------------------------
#  Configurazione studio Monte Carlo
#-------------------------------------------------------------
n_grid_mc <- c(50, 100, 200)
R_mc <- 200
tau_values <- c(0.3, 0.5, 0.7)
#-------------------------------------------------------------
#  Parallelizzazione (opzionale)
#-------------------------------------------------------------
  use_parallel <- TRUE
if (use_parallel) {
  n_cores <- max(1, detectCores() - 1)
  cat(sprintf("\nStudio Monte Carlo usando %d core\n\n", n_cores))
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c(
    "tau_to_theta_clayton", "theta_to_tau_clayton",
    "loglik_clayton", "log_prior_uniform_tau",
    "tau_from_z", "z_from_tau",
    "mcmc_rw_mh", "mcmc_tau_z", "mc_worker"
  ))
  clusterEvalQ(cl, {
    library(copula)
    library(coda)
  })
}
#-------------------------------------------------------------
#  Esecuzione studio Monte Carlo
#-------------------------------------------------------------
  all_results <- list()
for (tau_true_current in tau_values) {
  cat(sprintf("\nTau = %.1f\n", tau_true_current))
  results_mc <- list()
  for (n_current in n_grid_mc) {
    cat(sprintf("  n = %d...\n", n_current))
    
    if (use_parallel && exists("cl")) {
      clusterExport(cl, c("n_current", "tau_true_current"), 
                    envir = environment())
      
      res_list <- parLapply(cl, 1:R_mc, function(r) {
        seed_offset <- r * 100000 + round(tau_true_current * 1000) * 100 + 
          n_current
        mc_worker(tau_true_current, n_current,
                  n_iter_mc = 10000, n_burn_mc = 2000,
                  seed_offset = seed_offset)
      })
    } else {
      res_list <- lapply(1:R_mc, function(r) {
        seed_offset <- r * 100000 + round(tau_true_current * 1000) * 100 + 
          n_current
        mc_worker(tau_true_current, n_current,
                  n_iter_mc = 10000, n_burn_mc = 2000,
                  seed_offset = seed_offset)
      })
    }
    
    res_matrix <- do.call(rbind, res_list)
    
    results_mc[[as.character(n_current)]] <- data.frame(
      tau_true = tau_true_current,
      n = n_current,
      Bias = mean(res_matrix[, "mean"] - tau_true_current),
      RMSE = sqrt(mean((res_matrix[, "mean"] - tau_true_current)^2)),
      Coverage = mean(res_matrix[, "coverage"]),
      ESS = mean(res_matrix[, "ess"]),
      AccRate = mean(res_matrix[, "acc_rate"])
    )
  }
  all_results[[as.character(tau_true_current)]] <- do.call(rbind, results_mc)
}
#-------------------------------------------------------------
#  Tabella riassuntiva
#-------------------------------------------------------------
  table_mc_all <- do.call(rbind, all_results)
rownames(table_mc_all) <- NULL
cat("\n=== Risultati Studio Monte Carlo ===\n\n")
print(kable(table_mc_all, digits = 3))
#-------------------------------------------------------------
 # Visualizzazione: RMSE vs n
#-------------------------------------------------------------
  baseline <- table_mc_all %>%
  group_by(tau_true) %>%
  arrange(n) %>%
  slice(1) %>%
  transmute(tau_true, n0 = n, rmse0 = RMSE)
plot_df <- table_mc_all %>%
  left_join(baseline, by = "tau_true") %>%
  mutate(RMSE_theory = rmse0 * sqrt(n0 / n))
p_rmse <- ggplot(plot_df, aes(x = n, y = RMSE,
                              color = factor(tau_true),
                              group = factor(tau_true))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_line(aes(y = RMSE_theory), linetype = "dashed", alpha = 0.6) +
  scale_x_log10(breaks = n_grid_mc) +
  scale_y_log10() +
  labs(
    title = "Convergenza RMSE",
    subtitle = expression("Linee tratteggiate: rate teorico" ~
                            RMSE %prop% n^{-0.5}),
    x = "Dimensione campione (n)",
    y = "RMSE",
    color = expression(tau[0])
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
print(p_rmse)
#-------------------------------------------------------------
  #Visualizzazione: Coverage vs n
#-------------------------------------------------------------
  p_coverage <- ggplot(table_mc_all, aes(x = n, y = Coverage,
                                         color = factor(tau_true),
                                         group = factor(tau_true))) +
  geom_hline(yintercept = 0.95, linetype = "dashed",
             linewidth = 0.8, alpha = 0.6) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = n_grid_mc) +
  ylim(0.88, 1.0) +
  labs(
    title = "Copertura IC 95%",
    subtitle = paste0("Basata su ", R_mc, " replicazioni Monte Carlo"),
    x = "Dimensione campione (n)",
    y = "Copertura empirica",
    color = expression(tau[0])
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
print(p_coverage)
###############################################################################
# CONFRONTO FREQUENTISTA vs BAYESIANO
###############################################################################
#-------------------------------------------------------------
  #Worker per confronto metodi
#-------------------------------------------------------------
  worker_comparison <- function(r, tau_true, n, n_iter, n_burn, sigma, seed) {
    set.seed(seed + r)
    theta_true <- tau_to_theta_clayton(tau_true)
    Simula dati
    data <- rCopula(n, claytonCopula(theta_true, dim = 2))
    u <- data[, 1]
    v <- data[, 2]
   # --- Kendall ---
      tau_k <- cor(u, v, method = "kendall")
   # --- ML ---
      t_ml <- system.time({
        fit_ml <- fitCopula(claytonCopula(), cbind(u, v), method = "ml")
        theta_hat <- coef(fit_ml)
        tau_ml <- theta_to_tau_clayton(theta_hat)
        se_theta <- sqrt(as.numeric(vcov(fit_ml)))
        ci_theta <- theta_hat + c(-1, 1) * qnorm(0.975) * se_theta
        ci_tau <- theta_to_tau_clayton(ci_theta)
      })
    # --- Bayes ---
      t_b <- system.time({
        tau_init <- max(min(tau_k, 0.95), -0.95)
        z <- numeric(n_iter + n_burn)
        z[1] <- z_from_tau(tau_init)
        log_post_z_local <- function(z) {
          tau <- tau_from_z(z)
          if (tau <= -1 || tau >= 1) return(-Inf)
          log_lik_tau(tau, u, v) - log(2) + log(1 - tau^2)
        }
        
        for (t in 2:length(z)) {
          prop <- rnorm(1, z[t-1], sigma)
          loga <- log_post_z_local(prop) - log_post_z_local(z[t-1])
          z[t] <- if (log(runif(1)) < loga) prop else z[t-1]
        }
        
        tau_samp <- tau_from_z(z[(n_burn + 1):length(z)])
        tau_b <- mean(tau_samp)
        hpd_b <- HPDinterval(mcmc(tau_samp), 0.95)
      })
    c(
      tau_k = tau_k,
      tau_ml = tau_ml,
      ml_low = ci_tau[1],
      ml_up = ci_tau[2],
      time_ml = t_ml["elapsed"],
      tau_b = tau_b,
      b_low = hpd_b[1],
      b_up = hpd_b[2],
      time_b = t_b["elapsed"]
    )
  }
#-------------------------------------------------------------
#  Configurazione confronto
#-------------------------------------------------------------
  TAU_TRUE_COMP <- 0.5
N_COMP <- 200
R_COMP <- 500
N_ITER_COMP <- 10000
N_BURN_COMP <- 2000
SIGMA_COMP <- 0.6
SEED_COMP <- 123
N_CORES_COMP <- max(1, detectCores() - 1)
#-------------------------------------------------------------
#  Esecuzione parallela
#-------------------------------------------------------------
  cat("\nConfronto Frequentista vs Bayesiano\n")
cat(sprintf("Usando %d core\n\n", N_CORES_COMP))
cl_comp <- makeCluster(N_CORES_COMP)
on.exit(stopCluster(cl_comp), add = TRUE)
clusterExport(
  cl_comp,
  varlist = c(
    "worker_comparison", "log_lik_tau", "log_clayton_density",
    "tau_from_z", "z_from_tau", "tau_to_theta_clayton",
    "theta_to_tau_clayton",
    "TAU_TRUE_COMP", "N_COMP", "N_ITER_COMP",
    "N_BURN_COMP", "SIGMA_COMP", "SEED_COMP"
  ),
  envir = environment()
)
clusterEvalQ(cl_comp, {
  library(copula)
  library(coda)
})
res_list_comp <- parLapply(
  cl_comp, 1:R_COMP, worker_comparison,
  tau_true = TAU_TRUE_COMP,
  n = N_COMP,
  n_iter = N_ITER_COMP,
  n_burn = N_BURN_COMP,
  sigma = SIGMA_COMP,
  seed = SEED_COMP
)
res_comp <- do.call(rbind, res_list_comp)
res_comp <- as.data.frame(res_comp)
#-------------------------------------------------------------
#  Statistiche riassuntive
#-------------------------------------------------------------
  summ <- function(est, low, up, true_val) {
    c(
      Bias = mean(est - true_val),
      RMSE = sqrt(mean((est - true_val)^2)),
      Coverage = mean(true_val >= low & true_val <= up)
    )
  }
tab_comp <- rbind(
  "Kendall (MM)" = c(
    Bias = mean(res_comp[, "tau_k"] - TAU_TRUE_COMP),
    RMSE = sqrt(mean((res_comp[, "tau_k"] - TAU_TRUE_COMP)^2)),
    Coverage = NA,
    Time = 0.01
  ),
  "ML" = c(
    summ(res_comp[, "tau_ml"], res_comp[, "ml_low"],
         res_comp[, "ml_up"], TAU_TRUE_COMP),
    Time = mean(res_comp[, "time_ml"])
  ),
  "Bayes (Unif)" = c(
    summ(res_comp[, "tau_b"], res_comp[, "b_low"],
         res_comp[, "b_up"], TAU_TRUE_COMP),
    Time = mean(res_comp[, "time_b"])
  )
)
cat("\n=== Confronto Frequentista vs Bayesiano ===\n\n")
print(round(tab_comp, 3))