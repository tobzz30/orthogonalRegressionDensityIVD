#################################
### This code execute orthogonal density method for interval value data using 
## reconstructed center and range method.
## Result: We reveal the result using the MISE
#################################

orthogonal_density <- function(X,
                               method = c("truncation","oracle"),
                               J_max = 200,
                               grid_size = 512,
                               cT = sqrt(2)) {
  
  method <- match.arg(method)
  n <- length(X)
  
  # FOR CONSISTENCY
  J_max <- min(J_max, floor(n^(1/2)))
  
  phi <- function(j, x) {
    if (j == 0) return(rep(1, length(x)))
    sqrt(2) * cos(pi * j * x)
  }
  
  theta_hat <- sapply(0:J_max, function(j) mean(phi(j,X)))
  
  phi_vals <- sapply(0:J_max, function(j) phi(j,X))
  d_j <- rep(1, length(theta_hat))
  
  # DATA-DRIVEN CUT-OFF (Hart)
  cutoff_objective <- function(J) {
    sum(2*d_j[2:(J+1)]/n - theta_hat[2:(J+1)]^2)
  }
  
  J_hat <- which.min(sapply(1:J_max, cutoff_objective))
  J_hat <- min(max(J_hat, floor(n^(1/3)), floor(sqrt(n))))
  #J_hat <- max(J_hat, floor(n^(1/3)))
  
  theta_final <- theta_hat
  
  if(method == "truncation"){
    if(J_hat < J_max){
      theta_final[(J_hat+2):(J_max+1)] <- 0
    } ### since r =2
  }
  
  # if(method == "hardthresholding"){
  #   
  #   for(j in 1:J_hat){
  #     
  #     threshold <- cT * sqrt(d_j[j+1] * log(n) / n)
  #     
  #     if(abs(theta_hat[j+1]) <= threshold){
  #       theta_final[j+1] <- 0
  #     }
  #   }
  #   
  #   theta_final[(J_hat+2):(J_max+1)] <- 0
  # }
  
  if(method == "oracle"){
    for(j in 1:J_hat){
      theta_sq_hat <- theta_hat[j+1]^2
      theta_sq_unbiased <- max(0, theta_sq_hat - d_j[j+1]/n)
      
      denom <- theta_sq_unbiased + d_j[j+1]/n
      
      if(denom > 0){
        w_hat <- theta_sq_unbiased / denom
      } else {
        w_hat <- 0
      }
      
      theta_final[j+1] <- w_hat * theta_hat[j+1]
    }
    
    if(J_hat < J_max){
      theta_final[(J_hat+2):(J_max+1)] <- 0
    }
  }
  
  grid <- seq(0,1,length.out=grid_size)
  
  f_hat <- rowSums(sapply(0:J_max,
                          function(j) theta_final[j+1]*phi(j,grid)))
  
  # L2 PROJECTION
  dx <- grid[2]-grid[1] 
  
  
  projection <- function(c) sum(pmax(0,f_hat-c))*dx - 1
  
  c_hat <- uniroot(projection,
                   lower=min(f_hat)-1,
                   upper=max(f_hat))$root
  
  f_proj <- pmax(0,f_hat-c_hat)
  f_proj <- f_proj / sum(f_proj*dx)
  
  list(grid=grid,
       density=f_proj,
       J_hat=J_hat)
}

# Interval-valued orthogonal density estimator
orthogonal_interval_density_avg <- function(YL, YU,
                                            method = c("truncation",
                                                       "oracle"),
                                            J_max = 200,
                                            block_growth = 0.5,
                                            grid_size = 512,
                                            M = 100) {
  
  method <- match.arg(method)
  stopifnot(length(YL) == length(YU))
  
  n <- length(YL)
  
  # Center and Range
  
  C <- (YL + YU) / 2
  R <- (YU - YL) / 2
  stopifnot(all(R > 0))
  
  # Estimate Center Density (pilot for weight)
  
  C_min <- min(C)
  C_max <- max(C)
  C_rng <- C_max - C_min
  if(C_rng == 0) C_rng <- 1
  C_scaled <- (C - C_min) / C_rng
  
  fit_C <- orthogonal_density(C_scaled,method = method,
                              J_max = J_max,
                              grid_size = grid_size)
  
  # Estimate Range Density (log-scale)
  Z <- log(R)
  Z_min <- min(Z)
  Z_max <- max(Z)
  Z_rng <- Z_max - Z_min
  if(Z_rng == 0) Z_rng <- 1
  Z_scaled <- (Z - Z_min) / Z_rng
  
  fit_Z <- orthogonal_density(Z_scaled,
                              method = method,
                              J_max = J_max,
                              grid_size = grid_size)
  
  # Monte Carlo Averaging
  
  density_mat <- matrix(0, M, grid_size)
  Y_global_min <- min(YL)
  Y_global_max <- max(YU)
  range_Y      <- Y_global_max - Y_global_min
  dx_C <- fit_C$grid[2] - fit_C$grid[1]
  dx_Z <- fit_Z$grid[2] - fit_Z$grid[1]
  probC <- pmax(fit_C$density, 0)
  probC <- probC / sum(probC * dx_C)
  
  probZ <- pmax(fit_Z$density, 0)
  probZ <- probZ / sum(probZ * dx_Z)
  
  for (m in 1:M) {
    # Draw pseudo-sample of size n
    C_star_scaled <- sample(fit_C$grid, n, replace=TRUE, prob=probC)
    Z_star_scaled <- sample(fit_Z$grid, n, replace=TRUE, prob=probZ)
    
    # Back-transform
    C_star <- C_min + C_star_scaled * (C_max - C_min)
    R_star <- exp(Z_min + Z_star_scaled * (Z_max - Z_min))
    
    Y_star <- C_star + runif(n, -R_star, R_star)
    Y_star <- pmin(pmax(Y_star, Y_global_min), Y_global_max) 
    Y_scaled <- (Y_star - Y_global_min) / range_Y
    
    fit_tmp <- orthogonal_density(
      Y_scaled,
      method = method,
      J_max = J_max,
      grid_size = grid_size
    )
    
    # Back-transform density using GLOBAL Jacobian
    density_tmp <- fit_tmp$density / range_Y
    # Back-transform grid 
    y_grid_tmp  <- Y_global_min + fit_tmp$grid * range_Y
    
    density_interp <- approx(
      y_grid_tmp,
      density_tmp,
      xout = seq(Y_global_min, Y_global_max, length.out = grid_size),
      yleft = 0,
      yright = 0
    )$y
    
    density_mat[m, ] <- density_interp
  }
  
  # Average the M density estimates
  y_grid <- seq(Y_global_min, Y_global_max, length.out = grid_size)
  
  f_y <- colMeans(density_mat)
  
  dx <- y_grid[2] - y_grid[1]
  f_y <- f_y / sum(f_y * dx)
  
  y_grid <- seq(Y_global_min, Y_global_max, length.out = grid_size)
  
  # Center density backtransformation
  center_grid <- C_min + fit_C$grid * (C_max - C_min)
  center_density <- fit_C$density / (C_max - C_min)
  
  list(
    grid_interval = y_grid,
    density_interval = f_y,
    grid_center = center_grid,
    density_center = center_density,
    J_hat_center = fit_C$J_hat,
    J_hat_range  = fit_Z$J_hat,
    method = method
  )
}

MSE_integrated <- function(f_hat, f_true, dx) {
  y <- (f_hat - f_true)^2
  sum((y[-1] + y[-length(y)]) / 2) * dx
} ## Trapezoidal rule

plot_raw_intervals <- function(YL, YU, type_label = "") {
  
  n <- length(YL)
  
  plot(1:n, YL,
       type = "n",
       ylim = range(c(YL, YU)),
       xlab = "Observation index",
       ylab = "Value",
       main = paste("Raw Interval Data:", type_label))
  
  segments(1:n, YL, 1:n, YU,
           col = "blue", lwd = 2)
  
  points(1:n, (YL+YU)/2,
         col = "red", pch = 16)
}

corner_density <- function(x, type) {
  
  if (type == "uniform") {
    return(rep(1, length(x)))
  }
  
  if (type == "normal") {
    return(dnorm(x, 0.5, 0.15))
  }
  
  if (type == "bimodal") {
    return(0.5*dnorm(x,0.4,0.12) +
             0.5*dnorm(x,0.7,0.08))
  }
  
  if (type == "strata") {
    return(0.5*dnorm(x,0.2,0.06) +
             0.5*dnorm(x,0.7,0.08))
  }
  
  if (type == "delta") {
    return(dnorm(x,0.5,0.02))
  }
  
  if (type == "angle") {
    f <- numeric(length(x))
    f[x <= 0.5] <- dnorm(x[x <= 0.5],1,0.7)
    f[x > 0.5]  <- dnorm(x[x > 0.5],0,0.7)
    return(f)
  }
  
  if (type == "monotone") {
    norm_const <- 1/3 + 0.8 + 0.64   # exact integral
    return((x + 0.8)^2 / norm_const)
  }
  
  if (type == "steps") {
    f <- numeric(length(x))
    f[x < 1/3] <- 0.6
    f[x >= 1/3 & x < 3/4] <- 0.9
    f[x >= 3/4] <- 204/120
    return(f)
  }
}

sample_corner <- function(n, type) {
  
  grid <- seq(0,1,length.out=512)
  dens <- corner_density(grid,type)
  dens <- dens / sum(dens)
  
  sample(grid, n, replace=TRUE, prob=dens)
}

corner_types <- c("uniform","normal","bimodal",
                  "strata","delta","angle",
                  "monotone","steps")

sample_sizes <- c(100, 200, 300, 400, 500)
methods <- c("truncation","oracle")
B <- 1000

set.seed(100)

results_table <- data.frame()

for (type in corner_types) {
  
  for (method in methods) {
    
    for (n in sample_sizes) {
      
      cat("Running:", type, "| Method:", method, "| n =", n, "\n")
      
      MSE_center   <- numeric(B)
      MSE_interval <- numeric(B)
      
      for (b in 1:B) {
        
        Y_true <- sample_corner(n, type)
        r1 <- runif(n, 0.02, 0.12)
        r2 <- runif(n, 0.02, 0.12)
        YL <- pmax(0, Y_true - r1)
        YU <- pmin(1, Y_true + r2)
        
        fit <- orthogonal_interval_density_avg(
          YL, YU,
          method = method,
          M = 100
        )
        eval_grid <- seq(0, 1, length.out = 512)
        dx <- eval_grid[2] - eval_grid[1]
        
        f_true_vals <- corner_density(eval_grid, type)
        f_true_vals <- f_true_vals / (sum(f_true_vals) * dx)
        
        f_interval_fixed <- approx(
          fit$grid_interval,
          fit$density_interval,
          xout = eval_grid,
          yleft = 0, yright = 0   # density = 0 outside estimator support
        )$y
        
        f_center_fixed <- approx(
          fit$grid_center,
          fit$density_center,
          xout = eval_grid,
          yleft = 0, yright = 0
        )$y
        if(sum(f_interval_fixed)*dx > 0)
          f_interval_fixed <- f_interval_fixed / (sum(f_interval_fixed)*dx)
        
        if(sum(f_center_fixed)*dx > 0)
          f_center_fixed <- f_center_fixed / (sum(f_center_fixed)*dx)
        MSE_interval[b] <- MSE_integrated(f_interval_fixed, f_true_vals, dx)
        MSE_center[b]   <- MSE_integrated(f_center_fixed,   f_true_vals, dx)
      }
      
      # STORE RESULTS
      results_table <- rbind(
        results_table,
        data.frame(
          SampleSize = n,
          Density = type,
          Method = method,
          MISE_center = mean(MSE_center),
          MISE_interval = mean(MSE_interval)
        )
      )
    }
  }
}

print(results_table)


method_selected <- "truncation"   # change to "truncation" or "oracle"
results_m <- results_table[results_table$Method == method_selected, ]
densities <- unique(results_m$Density)

par(mfrow = c(2,4),
    mar = c(3,4,3,1),
    oma = c(3,0,4,0))   

for (d in densities) {
  
  sub <- results_m[results_m$Density == d, ]
  
  # sort for proper lines
  sub <- sub[order(sub$SampleSize), ]
  
  y_range <- range(c(sub$MISE_center,
                     sub$MISE_interval), na.rm = TRUE)
  
  plot(sub$SampleSize,
       sub$MISE_center,
       type = "b",
       pch = 16,
       lwd = 2,
       col = "red",
       ylim = y_range,
       xlab = "",
       ylab = "MISE",
       main = paste(d))
  
  lines(sub$SampleSize,
        sub$MISE_interval,
        type = "b",
        pch = 17,
        lwd = 2,
        col = "blue")
  
}

par(xpd = NA)
fig <- par("fig")
par(fig = c(0,1,0,1), new = TRUE)
plot.new()
mtext("n", 
      side = 1, 
      outer = TRUE, 
      line = 0.7,   # instead of 2.5
      cex = 1.0)
legend("top",
       legend = c("CM", "CRM"),
       col = c("red", "blue"),
       lwd = 3,
       pch = c(16,17),
       horiz = TRUE,
       bty = "n",
       cex = 1.2,
       inset = c(0, -0.20))
