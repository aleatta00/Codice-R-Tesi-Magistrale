###############################################################################
# ESEMPI CAPITOLO 2: LOG-VEROSIMIGLIANZA E DISTRIBUZIONI A POSTERIORI
# 
# Il presente codice genera le figure utilizzate negli Esempi del Capitolo 2
# per illustrare le differenze tra le parametrizzazioni in theta e tau
# della copula di Clayton.
#
# Esempio 3 (Sez. 2.2.2): Profilo di log-verosimiglianza
# Esempio 4 (Sez. 2.2.3): Confronto tra posteriori
###############################################################################

# -------------------------------------------------------------
# Librerie
# -------------------------------------------------------------
library(copula)
library(ggplot2)
library(gridExtra)

###############################################################################
# ESEMPIO 3: ASIMMETRIA DELLA LOG-VEROSIMIGLIANZA
# Sezione 2.2.2 - Figura 2.1
###############################################################################

# -------------------------------------------------------------
# Parametri di simulazione
# -------------------------------------------------------------
set.seed(123)

n_ex3 <- 500
tau_true_ex3 <- 0.5
theta_true_ex3 <- 2  # Corrispondente a tau = 0.5 per Clayton

# -------------------------------------------------------------
# Generazione dati dalla copula di Clayton
# -------------------------------------------------------------
cop_ex3 <- claytonCopula(theta_true_ex3)
data_ex3 <- rCopula(n_ex3, cop_ex3)

# -------------------------------------------------------------
# Calcolo della log-verosimiglianza su griglia di valori
# -------------------------------------------------------------
theta_grid_ex3 <- seq(0.05, 6, length.out = 300)

loglik_ex3 <- sapply(theta_grid_ex3, function(theta) {
  sum(dCopula(data_ex3, claytonCopula(theta), log = TRUE))
})

# Identificazione del massimo
idx_max <- which.max(loglik_ex3)
theta_mle <- theta_grid_ex3[idx_max]
loglik_max <- loglik_ex3[idx_max]

# -------------------------------------------------------------
# Visualizzazione
# -------------------------------------------------------------
df_ex3 <- data.frame(theta = theta_grid_ex3, loglik = loglik_ex3)

p_ex3 <- ggplot(df_ex3, aes(theta, loglik)) +
  geom_line(linewidth = 0.9, color = "black") +
  geom_vline(xintercept = theta_true_ex3, 
             linetype = "dashed", color = "red", linewidth = 1) +
  geom_point(data = data.frame(theta = theta_mle, loglik = loglik_max),
             aes(theta, loglik), color = "red", size = 3) +
  annotate("text", x = theta_true_ex3 + 0.6, y = loglik_max - 20,
           label = paste0("theta[true] == ", theta_true_ex3),
           parse = TRUE, color = "red", size = 4.5) +
  labs(
    title = "Log-verosimiglianza per la copula di Clayton",
    x = expression(theta),
    y = "Log-verosimiglianza"
  ) +
  theme_minimal(base_size = 16)

print(p_ex3)


###############################################################################
# ESEMPIO 4: CONFRONTO TRA POSTERIORI IN θ E τ
# Sezione 2.2.3 - Figura 2.2
###############################################################################

# -------------------------------------------------------------
# Parametri di simulazione
# -------------------------------------------------------------
set.seed(123)

n_ex4 <- 200
tau_true_ex4 <- 0.5
theta_true_ex4 <- 2

# -------------------------------------------------------------
# Generazione dati dalla copula di Clayton
# -------------------------------------------------------------
cop_ex4 <- claytonCopula(theta_true_ex4, dim = 2)
data_ex4 <- rCopula(n_ex4, cop_ex4)

u_ex4 <- data_ex4[, 1]
v_ex4 <- data_ex4[, 2]

# -------------------------------------------------------------
# Funzioni ausiliarie per il calcolo della log-verosimiglianza
# -------------------------------------------------------------

# Log-verosimiglianza in scala theta
loglik_theta <- function(theta, u, v) {
  if (theta <= 0) return(-Inf)
  cop <- claytonCopula(theta, dim = 2)
  sum(log(dCopula(cbind(u, v), cop)))
}

# Log-verosimiglianza in scala tau
loglik_tau <- function(tau, u, v) {
  if (tau <= -1 || tau >= 1) return(-Inf)
  theta <- 2 * tau / (1 - tau)  # Trasformazione tau -> theta per Clayton
  if (theta <= 0) return(-Inf)
  cop <- claytonCopula(theta, dim = 2)
  sum(log(dCopula(cbind(u, v), cop)))
}

# -------------------------------------------------------------
# Calcolo della distribuzione a posteriori in scala theta
# -------------------------------------------------------------
theta_grid <- seq(0.1, 6, length.out = 200)

loglik_theta_vals <- sapply(theta_grid, function(th) {
  loglik_theta(th, u_ex4, v_ex4)
})

# Normalizzazione (assumendo prior uniforme)
posterior_theta <- exp(loglik_theta_vals - max(loglik_theta_vals))
posterior_theta <- posterior_theta / (sum(posterior_theta) * diff(theta_grid)[1])

# -------------------------------------------------------------
# Calcolo della distribuzione a posteriori in scala tau
# -------------------------------------------------------------
tau_grid <- seq(0.05, 0.95, length.out = 200)

loglik_tau_vals <- sapply(tau_grid, function(t) {
  loglik_tau(t, u_ex4, v_ex4)
})

# Normalizzazione
posterior_tau <- exp(loglik_tau_vals - max(loglik_tau_vals))
posterior_tau <- posterior_tau / (sum(posterior_tau) * diff(tau_grid)[1])

# -------------------------------------------------------------
# Calcolo delle statistiche riassuntive
# -------------------------------------------------------------

# Funzione per calcolare media, varianza e skewness
compute_stats <- function(grid, density) {
  dt <- diff(grid)[1]
  mean_val <- sum(grid * density) * dt
  var_val <- sum((grid - mean_val)^2 * density) * dt
  skew_val <- sum((grid - mean_val)^3 * density) * dt / (var_val^(3/2))
  
  list(mean = mean_val, sd = sqrt(var_val), skewness = skew_val)
}

stats_theta <- compute_stats(theta_grid, posterior_theta)
stats_tau <- compute_stats(tau_grid, posterior_tau)

# -------------------------------------------------------------
# Visualizzazione: Posteriori in theta
# -------------------------------------------------------------
df_theta <- data.frame(theta = theta_grid, density = posterior_theta)

p_theta <- ggplot(df_theta, aes(theta, density)) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_vline(xintercept = theta_true_ex4, 
             linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = theta_true_ex4 + 0.7, 
           y = max(df_theta$density) * 0.9,
           label = paste0("theta[true] == ", theta_true_ex4),
           parse = TRUE, color = "red", size = 4) +
  labs(
    title = "Distribuzione a posteriori di θ",
    x = expression(theta),
    y = "Densità a posteriori"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# -------------------------------------------------------------
# Visualizzazione: Posteriori in tau
# -------------------------------------------------------------
df_tau <- data.frame(tau = tau_grid, density = posterior_tau)

p_tau <- ggplot(df_tau, aes(tau, density)) +
  geom_line(color = "darkgreen", linewidth = 1.2) +
  geom_vline(xintercept = tau_true_ex4, 
             linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = tau_true_ex4 - 0.12, 
           y = max(df_tau$density) * 0.9,
           label = paste0("tau[true] == ", tau_true_ex4),
           parse = TRUE, color = "red", size = 4) +
  labs(
    title = "Distribuzione a posteriori di τ",
    x = expression(tau),
    y = "Densità a posteriori"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# -------------------------------------------------------------
# Combinazione dei grafici
# -------------------------------------------------------------
combined_plot <- grid.arrange(p_theta, p_tau, ncol = 2)

# -------------------------------------------------------------
# Output delle statistiche riassuntive
# -------------------------------------------------------------
cat("\n=== STATISTICHE RIASSUNTIVE ===\n\n")

cat("Parametrizzazione in θ:\n")
cat(sprintf("  Media a posteriori: %.3f\n", stats_theta$mean))
cat(sprintf("  Deviazione standard: %.3f\n", stats_theta$sd))
cat(sprintf("  Valore vero: %.3f\n\n", theta_true_ex4))

cat("Parametrizzazione in τ:\n")
cat(sprintf("  Media a posteriori: %.3f\n", stats_tau$mean))
cat(sprintf("  Deviazione standard: %.3f\n", stats_tau$sd))
cat(sprintf("  Valore vero: %.3f\n\n", tau_true_ex4))

cat("Indici di asimmetria:\n")
cat(sprintf("  Skewness θ: %.3f\n", stats_theta$skewness))
cat(sprintf("  Skewness τ: %.3f\n", stats_tau$skewness))