#################################
### This code execute orthogonal regression method for interval value data using 
## reconstructed center and range method.
## Result: We examine the MSE for center method and center and range method.
#################################

legendre_basis <- function(x, k) {
  
  z <- 2*x - 1
  n <- length(x)
  
  Phi <- matrix(0, n, k+1)
  
  # Orthonormal Legendre basis
  Phi[,1] <- 1
  
  if (k >= 1) Phi[,2] <- sqrt(3) * z
  if (k >= 2) Phi[,3] <- sqrt(5) * (3*z^2 - 1)/2
  
  if (k >= 3) {
    P_prev <- z
    P_curr <- (3*z^2 - 1)/2
    
    for (j in 2:k) {
      P_next <- ((2*j-1)*z*P_curr - (j-1)*P_prev)/j
      
      # normalize
      Phi[,j+1] <- sqrt(2*j+1) * P_next
      
      P_prev <- P_curr
      P_curr <- P_next
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
      
      # Leave j out
      X_train <- X[-j]
      Y_train <- Y[-j]
      
      Phi_train <- legendre_basis(X_train, k)
      theta_hat <- tryCatch(
        qr.solve(Phi_train, Y_train),
        error = function(e) {
          # Ridge fallback
          lambda <- 1e-6
          solve(t(Phi_train) %*% Phi_train + lambda * diag(ncol(Phi_train))) %*% t(Phi_train) %*% Y_train
        }
      )
      if(any(!is.finite(theta_hat))) next
      # Predict left-out point
      Phi_test <- legendre_basis(X[j], k)
      Y_hat_j  <- Phi_test %*% theta_hat
      
      errors[j] <- (Y[j] - Y_hat_j)^2
    }
    
    lcv_vals[k] <- mean(errors)
  }
  
  which.min(ifelse(is.na(lcv_vals), Inf, lcv_vals))
}

scale01 <- function(x){
  rng <- max(x) - min(x)
  if(rng == 0) return(rep(0, length(x)))
  (x - min(x)) / rng
}


orthogonal_regression <- function(X, Y, k_max = 10) {
  
  n <- length(X)
  k_cv <- select_k_lcv(X, Y, k_max)
  #k_opt <- min(k_cv, 3)
  k_opt = 2 ## for simplicity
  Phi <- legendre_basis(X, k_opt)
  
  theta_hat <- tryCatch(
    qr.solve(Phi, Y),
    error = function(e) {
      lambda <- 1e-6
      solve(crossprod(Phi) + lambda * diag(ncol(Phi)),
            crossprod(Phi, Y))
    }
  )
  
  list(
    theta = theta_hat,
    k = k_opt
  )
}


estimate_error_from_predictions <- function(YL, YU, YL_hat, YU_hat) {
  
  e_L <- YL - YL_hat
  e_U <- YU - YU_hat
  
  e_center <- (e_L + e_U)/2
  e_range  <- (e_U - e_L)/2
  
  e_bar <- mean(e_center)
  sigma_hat <- sd(e_center)
  
  list(
    e_L = e_L,
    e_U = e_U,
    e_center = e_center,
    e_range = e_range,
    e_bar = e_bar,
    sigma_hat = sigma_hat
  )
}

center_method_np <- function(XL, XU, YL, YU, k_max = 10){
  Xc_raw <- (XL + XU)/2
  Yc_raw <- (YL + YU)/2
  xmin <- min(c(XL, XU))
  xmax <- max(c(XL, XU))
  xrange <- xmax - xmin
  if(xrange == 0) xrange <- 1
  
  Xc <- (Xc_raw - xmin) / xrange
  ymin <- min(Yc_raw)
  ymax <- max(Yc_raw)
  yrng <- ymax - ymin
  if(yrng == 0) yrng <- 1
  
  Yc <- (Yc_raw - ymin) / yrng
  fit <- orthogonal_regression(Xc, Yc, k_max)
  XL_s <- pmin(pmax((XL - xmin) / xrange, 0), 1)
  XU_s <- pmin(pmax((XU - xmin) / xrange, 0), 1)
  Phi_XL <- legendre_basis(XL_s, fit$k)
  Phi_XU <- legendre_basis(XU_s, fit$k)
  
  YL_scaled <- Phi_XL %*% fit$theta
  YU_scaled <- Phi_XU %*% fit$theta
  YL_raw <- ymin + yrng * YL_scaled
  YU_raw <- ymin + yrng * YU_scaled
  YL_hat <- pmin(YL_raw, YU_raw)
  YU_hat <- pmax(YL_raw, YU_raw)
  YL_hat <- pmin(pmax(YL_hat, min(YL)), max(YL))
  YU_hat <- pmin(pmax(YU_hat, min(YU)), max(YU))
  list(
    YL_hat = as.vector(YL_hat),
    YU_hat = as.vector(YU_hat),
    theta = fit$theta,
    k = fit$k
  )
}


crm_interval_regression <- function(XL, XU, YL, YU,
                                    k_max = 10,
                                    M = 100){
  
  n <- length(YL)
  
  Cy <- (YL + YU)/2
  Ry <- (YU - YL)
  
  xmin <- min(c(XL, XU))
  xmax <- max(c(XL, XU))
  
  xrange <- xmax - xmin
  if(xrange == 0) xrange <- 1
  
  XL_scaled <- pmin(pmax((XL - xmin)/xrange, 0), 1)
  XU_scaled <- pmin(pmax((XU - xmin)/xrange, 0), 1)
  
  Xc <- (XL_scaled + XU_scaled)/2
  Xr <- (XU_scaled - XL_scaled)/2
  
  ymin <- min(Cy)
  ymax <- max(Cy)
  yminr <- min(Ry)
  ymaxr <- max(Ry)
  yrng <- ymax - ymin
  yrngr <- ymaxr - yminr
  if(yrng == 0) yrng <- 1
  if(yrngr == 0) yrngr <- 1
  
  Cy_scaled <- (Cy - ymin)/yrng
  
  fit_center <- orthogonal_regression(Xc, Cy_scaled, k_max)
  Phi_Xc <- legendre_basis(Xc, fit_center$k)
  Cy_hat <- Phi_Xc %*% fit_center$theta
  Ry_scaled <- (Ry - yminr) / yrngr
  Ry_scaled <- pmax(Ry_scaled, 1e-6)
  fit_range <- orthogonal_regression(Xr, Ry_scaled, k_max)
  Phi_Xr <- legendre_basis(Xr, fit_range$k)
  Ry_hat <- abs(Phi_Xr %*% fit_range$theta)
  
  Cy_hat_s <- as.vector(scale(Cy_hat))
  Ry_hat_s <- as.vector(scale(Ry_hat))
  
  w_c <- abs(Cy_hat_s) + 1e-6
  w_c <- w_c / sum(w_c)
  w_r <- abs(Ry_hat_s) + 1e-6
  w_r <- w_r / sum(w_r)
  
  YL_hat_mat <- matrix(0, M, n)
  YU_hat_mat <- matrix(0, M, n)
  
  for(m in 1:M){
    
    idx_c <- sample(1:n, n, replace = TRUE, prob = w_c)
    idx_r <- sample(1:n, n, replace = TRUE, prob = w_r)
    
    Cy_star <- Cy_scaled[idx_c]
    #Ry_star <- Ry[idx]
    Xc_s <- Xc[idx_c]
    Xr_s <- Xr[idx_r]
    
    Ry_star <- Ry_scaled[idx_r]
    Ry_star <- pmax(Ry_star, 1e-6)
    
    Y_star <- Cy_star + runif(n, -Ry_star, Ry_star)
    fit_c <- orthogonal_regression(Xc_s, Y_star, k_max)
    
    
    fit_r <- orthogonal_regression(Xr_s, Y_star, k_max)
    
    
    Phi_c <- legendre_basis(Xc, fit_c$k)
    Phi_r <- legendre_basis(Xr, fit_r$k)
    
    Yhat_c <- Phi_c %*% fit_c$theta
    Yhat_r <- Phi_r %*% fit_r$theta
    YL_hat <- Yhat_c - Yhat_r/2
    YU_hat <- Yhat_c + Yhat_r/2
    if(any(!is.finite(YL_hat)) || any(!is.finite(YU_hat))) next
    YL_hat_mat[m,] <- pmin(YL_hat, YU_hat)
    YU_hat_mat[m,] <- pmax(YL_hat, YU_hat)
  }
  
  valid_rows <- apply(is.finite(YL_hat_mat), 1, all)
  mean_L <- colMeans(YL_hat_mat[valid_rows, , drop=FALSE])
  mean_U <- colMeans(YU_hat_mat[valid_rows, , drop=FALSE])
  
  mean_L <- ymin + yrng * mean_L
  mean_U <- ymin + yrng * mean_U
  
  YL_hat <- pmin(mean_L, mean_U)
  YU_hat <- pmax(mean_L, mean_U)
  
  list(
    YL_hat = as.vector(YL_hat),
    YU_hat = as.vector(YU_hat)
  )
}


predict_orthogonal <- function(x_new,
                               fit_object,
                               disturbance_object,
                               alpha = 0.05) {
  
  # Direct Legendre evaluation
  Phi_new <- legendre_basis(x_new, fit_object$k)
  h_hat   <- Phi_new %*% fit_object$theta
  
  # forecast mean (paper eq 11)
  u <- h_hat + disturbance_object$e_bar
  
  # Normal CI
  z_alpha <- qnorm(1 - alpha/2)
  
  lower <- u - z_alpha * disturbance_object$sigma_hat
  upper <- u + z_alpha * disturbance_object$sigma_hat
  
  list(
    mean_prediction = as.vector(u),
    lower = as.vector(lower),
    upper = as.vector(upper),
    sigma_hat = disturbance_object$sigma_hat
  )
}

#############################################
# Monte Carlo MSE_L and MSE_U comparison
#############################################



true_function <- function(x) {
  6+4*sin(0.25*pi*x)
}


compute_MSE <- function(y_true, y_pred) {
  mean((y_true - y_pred)^2)
}

sample_sizes <- c(100, 200, 300, 400, 500)

B <- 1000  # Monte Carlo repetitions
set.seed(100)
results <- data.frame()

for (n in sample_sizes) {
  
  MSE_L_CM  <- numeric(B)
  MSE_U_CM  <- numeric(B)
  
  MSE_L_CRM <- numeric(B)
  MSE_U_CRM <- numeric(B)
  k_max <- min(10, floor(n^(1/3)))
  for (b in 1:B) {
    
    ####################################
    # Generate interval X
    ####################################
    X_raw <- rnorm(n, 5, 2^2)
    X <- scale01(X_raw)
    r1x <- runif(n, 0.01, 3)
    r2x <- runif(n, 0.01, 3)
    
    XL <- X - r1x
    XU <- X + r2x
    
    XL <- pmax(0, XL)
    XU <- pmin(1, XU)
    ####################################
    # Generate true Y
    ####################################
    
    Y_true <- true_function(X)
    
    eps <- rnorm(n, 0, 0.3)
    Y <- Y_true+eps
    
    r1y <- runif(n, 0.01, 3)
    r2y <- runif(n, 0.01, 3)
    
    YL <- Y - r1y
    YU <- Y + r2y
    
    ####################################
    # Fit CM
    fit_cm <- center_method_np(XL, XU, YL, YU, k_max = k_max)
    
    ####################################
    # Fit CRM
    ####################################
    
    fit_crm <- crm_interval_regression(XL, XU, YL, YU, k_max = k_max)
    
    ####################################
    # Compute MSE_L and MSE_U
    ####################################
    # CM error
    err_cm <- estimate_error_from_predictions(
      YL, YU,
      fit_cm$YL_hat,
      fit_cm$YU_hat
    )
    
    # CRM error
    err_crm <- estimate_error_from_predictions(
      YL, YU,
      fit_crm$YL_hat,
      fit_crm$YU_hat
    )
    
    MSE_L_CM[b]  <- compute_MSE(YL, fit_cm$YL_hat)
    MSE_U_CM[b]  <- compute_MSE(YU, fit_cm$YU_hat)
    MSE_L_CRM[b] <- compute_MSE(YL, fit_crm$YL_hat)
    MSE_U_CRM[b] <- compute_MSE(YU, fit_crm$YU_hat)
  }
  
  results <- rbind(results,
                   data.frame(
                     n = n,
                     MSE_L_CM  = mean(MSE_L_CM),
                     MSE_U_CM  = mean(MSE_U_CM),
                     MSE_L_CRM = mean(MSE_L_CRM),
                     MSE_U_CRM = mean(MSE_U_CRM)
                   ))
}

results

library(tidyverse)

write_csv(results, "result.csv")

