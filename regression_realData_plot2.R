#################################
### This code execute orthogonal regression method for interval value curve for 
## reconstructed center and range method and center method 
## Result: This generate curves based on climate dataset
#################################

set.seed(100)

legendre_basis <- function(x, k) {
  x_scaled <- (x - xmin_global) / (xmax_global - xmin_global)
  x_scaled <- pmin(pmax(x_scaled, 0), 1)
  z <- 2*x_scaled - 1
  n <- length(x)
  Phi <- matrix(0, n, k+1)
  Phi[,1] <- 1
  if (k >= 1) Phi[,2] <- sqrt(3) * z
  if (k >= 2) Phi[,3] <- sqrt(5) * (3*z^2 - 1)/2
  if (k >= 3) {
    P_prev <- z; P_curr <- (3*z^2 - 1)/2
    for (j in 2:k) {
      P_next <- ((2*j-1)*z*P_curr - (j-1)*P_prev)/j
      Phi[,j+1] <- sqrt(2*j+1) * P_next
      P_prev <- P_curr; P_curr <- P_next
    }
  }
  return(Phi)
}

# Special basis for already-scaled inputs (no global rescaling)
legendre_basis_scaled <- function(x_scaled, k) {
  x_scaled <- pmin(pmax(x_scaled, 0), 1)
  z <- 2*x_scaled - 1
  n <- length(x_scaled)
  Phi <- matrix(0, n, k+1)
  Phi[,1] <- 1
  if (k >= 1) Phi[,2] <- sqrt(3) * z
  if (k >= 2) Phi[,3] <- sqrt(5) * (3*z^2 - 1)/2
  if (k >= 3) {
    P_prev <- z; P_curr <- (3*z^2 - 1)/2
    for (j in 2:k) {
      P_next <- ((2*j-1)*z*P_curr - (j-1)*P_prev)/j
      Phi[,j+1] <- sqrt(2*j+1) * P_next
      P_prev <- P_curr; P_curr <- P_next
    }
  }
  return(Phi)
}

select_k_lcv <- function(X, Y, k_max = 10) {
  n <- length(X)
  lcv_vals <- numeric(k_max)
  for (k in 1:k_max) {
    errors <- numeric(n)
    for (j in 1:n) {
      X_train <- X[-j]; Y_train <- Y[-j]
      Phi_train <- legendre_basis_scaled(X_train, k)
      theta_hat <- tryCatch(
        qr.solve(Phi_train, Y_train),
        error = function(e) {
          lambda <- 1e-6
          solve(t(Phi_train) %*% Phi_train + lambda*diag(ncol(Phi_train))) %*% t(Phi_train) %*% Y_train
        }
      )
      if(any(!is.finite(theta_hat))) next
      Phi_test <- legendre_basis_scaled(X[j], k)
      errors[j] <- (Y[j] - Phi_test %*% theta_hat)^2
    }
    lcv_vals[k] <- mean(errors)
  }
  which.min(ifelse(is.na(lcv_vals), Inf, lcv_vals))
}

orthogonal_regression_scaled <- function(X_sc, Y, k_max = 10) {
  k_cv  <- select_k_lcv(X_sc, Y, k_max)
  k_opt <- min(k_cv, floor(length(X_sc)^(1/3)))
  Phi   <- legendre_basis_scaled(X_sc, k_opt)
  theta_hat <- tryCatch(
    qr.solve(Phi, Y),
    error = function(e) {
      lambda <- 1e-6
      solve(crossprod(Phi) + lambda*diag(ncol(Phi)), crossprod(Phi, Y))
    }
  )
  list(theta = theta_hat, k = k_opt)
}

# CRM
crm_interval_regression <- function(XL_sc, XU_sc, YL, YU, k_max = 10, M = 100) {
  n   <- length(YL)
  Cy  <- (YL + YU) / 2
  Ry  <- (YU - YL) / 2
  Xc  <- (XL_sc + XU_sc) / 2
  Xr  <- (XU_sc - XL_sc) / 2
  Xr_sc2 <- (Xr - min(Xr))/(max(Xr)-min(Xr))
  ymin <- min(Cy); ymax <- max(Cy); yrng <- ymax - ymin
  if(yrng == 0) yrng <- 1
  yminr <- min(Ry); ymaxr <- max(Ry); yrngr <- ymaxr - yminr
  if(yrngr == 0) yrngr <- 1
  
  Cy_sc <- (Cy - ymin) / yrng
  Ry_sc <- pmax((Ry - yminr) / yrngr, 1e-6)
  
  fit_c0   <- orthogonal_regression_scaled(Xc, Cy_sc, k_max)
  Cy_hat_s <- as.vector(scale(legendre_basis_scaled(Xc, fit_c0$k) %*% fit_c0$theta))
  fit_r0 <- orthogonal_regression_scaled(Xr_sc2, Ry_sc, k_max)
  Ry_hat_s <- as.vector(scale(abs(legendre_basis_scaled(Xr, fit_r0$k) %*% fit_r0$theta)))
  
  w_c <- abs(Cy_hat_s) + 1e-6; w_c <- w_c / sum(w_c)
  w_r <- abs(Ry_hat_s) + 1e-6; w_r <- w_r / sum(w_r)
  
  # Grid in scaled [0,1] space
  x_grid_sc <- seq(min(Xc_sc), max(Xc_sc), length.out=300)
  YL_mat     <- matrix(0,  M, n)
  YU_mat     <- matrix(0,  M, n)
  grid_mat_c   <- matrix(NA, M, 300)
  grid_mat_r <- matrix(NA, M, 300)
  h_hat_l <- matrix(NA, M, 300)
  h_hat_u <- matrix(NA, M, 300)
  
  for(m in 1:M) {
    ic <- sample(n, replace=TRUE, prob=w_c)
    ir <- sample(n, replace=TRUE, prob=w_r)
    Ry_star <- pmax(Ry_sc[ir], 1e-6)
    Y_star  <- Cy_sc[ic] + runif(n, -Ry_star, Ry_star)
    
    fc <- orthogonal_regression_scaled(Xc[ic], Y_star, k_max)
    fr <- orthogonal_regression_scaled(Xr_sc2[ir], Y_star, k_max)
    
    Yhat_c <- legendre_basis_scaled(Xc, fc$k) %*% fc$theta
    Yhat_r <- legendre_basis_scaled(Xr, fr$k) %*% fr$theta
    YLh <- Yhat_c - abs(Yhat_r)/2
    YUh <- Yhat_c + abs(Yhat_r)/2
    if(any(!is.finite(YLh)) || any(!is.finite(YUh))) next
    
    YL_mat[m,] <- pmin(YLh, YUh)
    YU_mat[m,] <- pmax(YLh, YUh)
    
    # Predict on grid in same scaled space, back-transform to Y
    Phi_gc        <- legendre_basis_scaled(x_grid_sc, fc$k)
    xr_grid_sc2 <- seq(min(Xr_sc), max(Xr_sc), length.out=300)
    Phi_gr <- legendre_basis_scaled(xr_grid_sc2, fr$k)
    grid_mat_c[m,] <- ymin + yrng * as.vector(Phi_gc %*% fc$theta)
    grid_mat_r[m,] <- yminr + yrngr * as.vector(Phi_gr %*% fr$theta)
    h_hat_l[m,] <- grid_mat_c[m,]-grid_mat_r[m,]
    h_hat_u[m,] <- grid_mat_c[m,]+grid_mat_r[m,]
  }
  
  valid  <- apply(is.finite(YL_mat) & YL_mat != 0, 1, all)
  mean_L <- ymin + yrng * colMeans(YL_mat[valid,,drop=FALSE])
  mean_U <- ymin + yrng * colMeans(YU_mat[valid,,drop=FALSE])
  x_grid_sc <- seq(0,1,length.out=300)
  x_grid_real <- xmin_global + x_grid_sc*(xmax_global - xmin_global)
  y_grid_avg_l <- colMeans(h_hat_l[valid,,drop=FALSE], na.rm=TRUE)
  y_grid_avg_u <- colMeans(h_hat_u[valid,,drop=FALSE], na.rm=TRUE)
  
  list(
    YL_hat = as.vector(pmin(mean_L, mean_U)),
    YU_hat = as.vector(pmax(mean_L, mean_U)),
    x_grid = x_grid_real, 
    y_grid_l = y_grid_avg_l,    # length 300, on original Y scale
    y_grid_u = y_grid_avg_u,
    y_grid_avg = (y_grid_avg_l+y_grid_avg_u)/2
  )
}



#############################################
###############  Real life dataset
#############################################


dat <- read.csv("daily_intervals_q1_q3.csv")

YL <- dat$pressure_q1
YU <- dat$pressure_q3

######################## for wind speed
XL <- dat$wind_q1
XU <- dat$wind_q3

######################## for temperature speed
# XL <- dat$temp_q1
# XU <- dat$temp_q3

Xc11 <- (XL + XU)/2
# scale X to [0,1]
xmin_global <- min(XL)
xmax_global <- max(XU)
### for wind speed only 
XL_sc <- (XL - xmin_global)/(xmax_global - xmin_global)
XU_sc <- (XU - xmin_global)/(xmax_global - xmin_global)

Xc_sc <- (XL_sc + XU_sc)/2
Xr_sc <- (XU_sc - XL_sc)/2

Yc <- (YL + YU)/2

# CM
fit_cm_obj <- orthogonal_regression_scaled(Xc_sc, Yc)
# CRM
fit_crm <- crm_interval_regression(XL_sc, XU_sc, YL, YU)

x_grid <- seq(xmin_global, xmax_global, length.out = 300)
x_grid1 <- seq(min(Xc11), max(Xc11), length.out = 300)

# CM curve
Phi_cm_grid <- legendre_basis(x_grid1, fit_cm_obj$k)
y_cm_grid   <- as.vector(Phi_cm_grid %*% fit_cm_obj$theta)


par(mar = c(4.5, 4.5, 3, 2), bg = "white")

plot(NA,
     xlim = range(fit_crm$x_grid),
     ylim = range(c(YL, YU), na.rm=TRUE),
     xlab = "Wind Speed",
     ylab = "Surface Pressure",
     main = "",
     las = 1)

# rectangles
for(i in 1:length(YL)) {
  rect(XL[i], YL[i], XU[i], YU[i],
       col    = rgb(0.7,0.7,0.7,0.05),
       border = rgb(0.5,0.5,0.5,0.2))
}

# points
points((XL+XU)/2, (YL+YU)/2,
       pch = 16, cex = 0.5, col = "grey40")

# CRM
lines(fit_crm$x_grid, fit_crm$y_grid_avg,
      col = "blue", lwd = 2)
#cm
lines(x_grid, y_cm_grid,
      col = "red", lwd = 2, lty = 2)

legend("topright",
       legend = c("CRM", "CM"),
       col    = c("blue", "red"),
       lty    = c(1, 2),
       lwd    = 2,
       bty    = "n")

######### Real data diagnostic for outlier detection

library(ggplot2)

r_w <- dat$wind_q3 - dat$wind_q1
r_p <- dat$pressure_q3 - dat$pressure_q1

hist(r_p, main = "", xlab = "Range of Surface Pressure")
hist(r_w, main = "", xlab = "Range of Wind Speed")
hist(dat$temp_q3-dat$temp_q1, main = "", xlab = "Range of Temperature")
df_plot <- data.frame(
  WindRange = r_w,
  PressureRange = r_p
)

ggplot(df_plot, aes(x = WindRange, y = PressureRange)) +
  geom_point(size = 1.5, alpha = 0.5, color = "blue") +
  geom_smooth(method = "loess", span = 0.6,
              color = "black", linewidth = 1.2, se = FALSE) +
  labs(
    title = "",
    x = "Wind Range",
    y = "Pressure Range"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )
