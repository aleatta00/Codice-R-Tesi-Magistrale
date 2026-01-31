###############################################################
# Appendice – Codice R per le visualizzazioni del Capitolo 1
# Copule: funzioni di ripartizione C(u,v) e densità c(u,v)
# Griglia uniforme su (0,1)^2 con trimming ai bordi per stabilità.
###############################################################

library(copula)
library(ggplot2)
library(plotly)

# -------------------------------------------------------------
# Griglia su (0,1)^2 (evita 0 e 1 per stabilità delle densità)
# -------------------------------------------------------------
n_points <- 100
eps <- 1e-4

u <- seq(eps, 1 - eps, length.out = n_points)
v <- seq(eps, 1 - eps, length.out = n_points)

grid <- expand.grid(u = u, v = v)
uv <- as.matrix(grid)

as_surface_matrix <- function(values) {
  matrix(values, nrow = n_points, ncol = n_points)  # ordine coerente con expand.grid
}

# -------------------------------------------------------------
# Utility: pulizia valori non finiti (utile per densità)
# -------------------------------------------------------------
clean_values <- function(x) {
  x[!is.finite(x)] <- NA_real_
  x
}

# -------------------------------------------------------------
# Plot 2D (contour filled); opzionale trasformazione log per densità
# -------------------------------------------------------------
plot2D <- function(values, title, log_scale = FALSE) {
  df <- grid
  values <- clean_values(values)
  if (log_scale) values <- log1p(values)
  df$z <- values
  
  ggplot(df, aes(u, v, z = z)) +
    geom_contour_filled(na.rm = TRUE) +
    labs(title = title, x = "u", y = "v") +
    theme_minimal()
}

# -------------------------------------------------------------
# Plot 3D
# -------------------------------------------------------------
plot3D <- function(values, title, log_scale = FALSE) {
  values <- clean_values(values)
  if (log_scale) values <- log1p(values)
  
  plot_ly(
    x = u, y = v, z = as_surface_matrix(values),
    type = "surface"
  ) %>%
    layout(
      title = title,
      scene = list(
        xaxis = list(title = "u"),
        yaxis = list(title = "v"),
        zaxis = list(title = ifelse(log_scale, "log(1+z)", "z"))
      )
    )
}

# -------------------------------------------------------------
# 1) Copula indipendenza e bounds di Fréchet–Hoeffding (solo C)
# -------------------------------------------------------------
C_ind <- grid$u * grid$v
C_M   <- pmin(grid$u, grid$v)
C_W   <- pmax(grid$u + grid$v - 1, 0)

plot2D(C_ind, "Copula di indipendenza")
plot2D(C_M,   "Copula M (comonotona)")
plot2D(C_W,   "Copula W (contromonotona)")

plot3D(C_ind, "Copula di indipendenza – 3D")
plot3D(C_M,   "Copula M – 3D")
plot3D(C_W,   "Copula W – 3D")

# -------------------------------------------------------------
# 2) Copule ellittiche: Gaussiana e t-Student
# -------------------------------------------------------------
rho <- 0.8
nu  <- 3

cop_gauss <- normalCopula(param = rho)
cop_t     <- tCopula(param = rho, df = nu)

C_gauss <- pCopula(uv, cop_gauss)
D_gauss <- dCopula(uv, cop_gauss)

C_t <- pCopula(uv, cop_t)
D_t <- dCopula(uv, cop_t)

plot2D(C_gauss, "Copula Gaussiana")
plot2D(D_gauss, "Densità della copula Gaussiana", log_scale = TRUE)

plot2D(C_t, "Copula t-Student")
plot2D(D_t, "Densità della copula t-Student", log_scale = TRUE)

plot3D(C_gauss, "Copula Gaussiana – 3D")
plot3D(D_gauss, "Densità della copula Gaussiana – 3D", log_scale = TRUE)

plot3D(C_t, "Copula t-Student – 3D")
plot3D(D_t, "Densità della copula t-Student – 3D", log_scale = TRUE)

# -------------------------------------------------------------
# 3) Copule archimedee (parametrizzate via tau di Kendall)
# -------------------------------------------------------------
tau <- 0.5

theta_clayton <- 2 * tau / (1 - tau)
theta_gumbel  <- 1 / (1 - tau)
theta_frank   <- iTau(frankCopula(), tau)

arch <- list(
  Clayton = claytonCopula(theta_clayton),
  Gumbel  = gumbelCopula(theta_gumbel),
  Frank   = frankCopula(theta_frank)
)

for (name in names(arch)) {
  cop <- arch[[name]]
  Cval <- pCopula(uv, cop)
  Dval <- dCopula(uv, cop)
  
  plot2D(Cval, paste("Copula", name))
  plot2D(Dval, paste("Densità della copula", name), log_scale = TRUE)
  
  plot3D(Cval, paste("Copula", name, "– 3D"))
  plot3D(Dval, paste("Densità della copula", name, "– 3D"), log_scale = TRUE)
}
