#################################
### This code execute orthogonal regression method for interval value curve for 
## reconstructed center and range method and center method 
## We compare this with KDE as benchmark since we don't know the true curve for real data
## Result: This generate regression function curves based on climate dataset and simulation settings
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
  #d_j <- apply(phi_vals, 2, var) * n
  d_j <- rep(1, length(theta_hat))
  # -----------------------------
  # DATA-DRIVEN CUT-OFF (Hart)
  # -----------------------------
  cutoff_objective <- function(J) {
    sum(2*d_j[2:(J+1)]/n - theta_hat[2:(J+1)]^2)
  }
  
  J_hat <- which.min(sapply(1:J_max, cutoff_objective))
  J_hat <- max(J_hat, floor(n^(1/3)))
  
  theta_final <- theta_hat
  
  # =============================
  # TRUNCATION
  # =============================
  if(method == "truncation"){
    if(J_hat < J_max){
      theta_final[(J_hat+2):(J_max+1)] <- 0
    } ### since r =2
  }
  
  # =============================
  # ORACLE (CORRECT FORM)
  # =============================
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
  
  # =============================
  # L2 PROJECTION
  # =============================
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

############################################################
# Interval-valued orthogonal density estimator
# Averaged Monte Carlo version
############################################################
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
  
  ########################################################
  # Center and Range
  ########################################################
  
  C <- (YL + YU) / 2
  R <- (YU - YL) / 2
  stopifnot(all(R > 0))
  
  ########################################################
  # Estimate Center Density
  ########################################################
  
  C_min <- min(C)
  C_max <- max(C)
  C_rng <- C_max - C_min
  if(C_rng == 0) C_rng <- 1
  C_scaled <- (C - C_min) / C_rng
  
  fit_C <- orthogonal_density(C_scaled,method = method,
                              J_max = J_max,
                              grid_size = grid_size)
  
  ########################################################
  # Estimate Range Density (log-scale)
  ########################################################
  
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
  
  ########################################################
  # Monte Carlo Averaging
  ########################################################
  
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
    
    # Back-transform grid using GLOBAL support
    y_grid_tmp  <- Y_global_min + fit_tmp$grid * range_Y
    
    density_interp <- approx(
      y_grid_tmp,
      density_tmp,
      xout = seq(Y_global_min, Y_global_max, length.out = grid_size),
      rule = 2
    )$y
    
    density_mat[m, ] <- density_interp
  }
  
  ########################################################
  # Average the M density estimates
  ########################################################
  
  y_grid <- seq(Y_global_min, Y_global_max, length.out = grid_size)
  
  f_y <- colMeans(density_mat)
  
  dx <- y_grid[2] - y_grid[1]
  f_y <- f_y / sum(f_y * dx)
  
  # Use grid from last iteration (all identical length)
  #y_grid <- seq(min(YL), max(YU), length.out = grid_size)
  y_grid <- seq(Y_global_min, Y_global_max, length.out = grid_size)
  #y_grid <- seq(Y_global_min, Y_global_max, length.out = grid_size)
  ########################################################
  # Center density (baseline for comparison)
  ########################################################
  
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


interval_kde <- function(YL, YU, grid, h = NULL, var_floor = 1e-3) {
  
  n <- length(YL)
  
  m <- (YL + YU) / 2
  v <- (YU - YL)^2 / 12
  
  # Silverman bandwidth 
  if (is.null(h)) {
    h_silverman <- 1.06 * min(sd(m), IQR(m)/1.34) * n^(-1/5)
    h <- 0.5 * h_silverman
  }
  
  f_kde <- rep(0, length(grid))
  
  for (i in 1:n) {
    sigma <- sqrt(pmax(v[i], var_floor) + h^2)
    f_kde <- f_kde + dnorm(grid, mean = m[i], sd = sigma)
  }
  
  f_kde <- f_kde / n
  
  dx <- grid[2] - grid[1]
  f_kde <- f_kde / sum(f_kde * dx)
  
  return(list(f = f_kde, h = h))
}

plot_kde_vs_estimators <- function(YL, YU, method, type_label, n_label) {
  
  fit <- orthogonal_interval_density_avg(YL, YU, method = method, M = 100)
  
  grid <- fit$grid_interval
  dx   <- grid[2] - grid[1]
  
  #  KDE (truth proxy) 
  f_true <- corner_density(grid, type_label)
  dx <- grid[2] - grid[1]
  f_true <- f_true / (sum(f_true) * dx)
  # ── CRM ──
  f_crm <- fit$density_interval
  
  # ── CM (interpolated) ──
  f_cm <- approx(fit$grid_center,
                 fit$density_center,
                 xout = grid,
                 rule = 2)$y
  
  f_cm <- pmax(f_cm, 0)
  f_cm <- f_cm / sum(f_cm * dx)
  
  ##################################################
  # KDE vs CRM
  ##################################################
  
  plot(grid, f_true,
       type = "l", lwd = 3, lty = 2, col = "black",
       ylim = c(0, max(f_true, f_crm) * 1.1),
       main = paste(type_label, "| n =", n_label, "|", method, "\nTRUE vs CRM"),
       xlab = "x", ylab = "Density")
  
  lines(grid, f_crm, col = "blue", lwd = 2)
  
  legend("topright",
         legend = c("True density", "CRM"),
         col = c("black", "blue"),
         lty = c(2,1),
         lwd = 2,
         bty = "n")
  
  ##################################################
  # KDE vs CM
  ##################################################
  
  plot(grid, f_true,
       type = "l", lwd = 3, lty = 2, col = "black",
       ylim = c(0, max(f_true, f_cm) * 1.1),
       main = paste(type_label, "| n =", n_label, "|", method, "\nTRUE vs CM"),
       xlab = "x", ylab = "Density")
  
  lines(grid, f_cm, col = "red", lwd = 2)
  
  legend("topright",
         legend = c("True density", "CM"),
         col = c("black", "red"),
         lty = c(2,1),
         lwd = 2,
         bty = "n")
}

set.seed(100)

# Create output directory, make sure to use your directory if run locally.
# by default folder is created under the current directory
#out_dir <- "/Users/oluwatobiakinbode/Documents/PHD Classes/Research/Dissertation Docs/Spring 2026/density_plot_orthorgonal"
out_dir <- file.path(path.expand("./"), "density_plot_orthorgonal")
dir.create(out_dir, showWarnings = FALSE)

set.seed(100)
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

method_selected <- "oracle"

for (type in corner_types) {
  
  file_name <- file.path(out_dir,
                         paste0("density_", type, "_", method_selected, ".png"))
  
  png(file_name,
      width = 2000, height = 2500, res = 150)
  
  par(mfrow = c(length(sample_sizes), 2),
      mar = c(4,4,3,1),
      oma = c(0,0,3,0))
  
  first_plot_done <- FALSE  
  
  for (n in sample_sizes) {
    
    Y_true <- sample_corner(n, type)
    
    r1 <- runif(n, 0.02, 0.12)
    r2 <- runif(n, 0.02, 0.12)
    YL <- pmax(0, Y_true - r1)
    YU <- pmin(1, Y_true + r2)
    
    plot_kde_vs_estimators(YL, YU,
                           method = method_selected,
                           type_label = type,
                           n_label = n)
    
    if (!first_plot_done) {
      mtext(paste("Density:", toupper(type),
                  "| Method:", toupper(method_selected)),
            outer = TRUE, cex = 1.5, font = 2)
      
      first_plot_done <- TRUE
    }
  }
  
  dev.off()
}



######################################
#### Real Life Dataset
######################################

dat <- read.csv("daily_intervals_q1_q3.csv")

YL <- dat$temp_q1 
YU <- dat$temp_q3


########## Put them together 
dat <- read.csv("daily_intervals_q1_q3.csv")

YL <- dat$pressure_q1  ####### for temperature and wind speed, you need to call them out here
YU <- dat$pressure_q3

plot_realdata_kde_vs_estimators <- function(YL, YU, method_label) {
  
  fit <- orthogonal_interval_density_avg(YL, YU,
                                         method = method_label,
                                         M = 100)
  
  grid <- fit$grid_interval
  dx   <- grid[2] - grid[1]
  
  # KDE (reference)
  kde_obj <- interval_kde(YL, YU, grid)
  f_kde <- kde_obj$f
  
  # CRM
  f_crm <- fit$density_interval
  
  # CM
  f_cm <- approx(fit$grid_center,
                 fit$density_center,
                 xout = grid,
                 rule = 2)$y
  
  f_cm <- pmax(f_cm, 0)
  f_cm <- f_cm / sum(f_cm * dx)
  
  y_max <- max(f_kde, f_crm, f_cm)
  
  ##################################################
  # COMBINED PLOT
  ##################################################
  bar_centers <- grid
  bar_width   <- (grid[2] - grid[1]) * 0.9
  
  idx <- seq(1, length(bar_centers), length.out = 40)
  
  plot(bar_centers[idx], f_kde[idx],
       type = "h",
       lwd = 5,
       lend = "butt",
       col = "black",
       ylim = c(0, y_max*1.1),
       main = paste(method_label),
       xlab = "Surface Pressure",
       ylab = "Density")
  
  # KDE (dashed)
  lines(grid, f_kde, col = "black", lwd = 2, lty = 2)
  
  # CRM
  lines(grid, f_crm, col = "blue", lwd = 2)
  
  # CM
  lines(grid, f_cm, col = "red", lwd = 2)
  
  legend("topright",
         legend = c("KDE", "CRM", "CM"),
         col = c("black", "blue", "red"),
         lty = c(2, 1, 1),
         lwd = 2,
         bty = "n")
}


pdf(file.path(out_dir, "realdata_sp_combined.pdf"),
    width = 10, height = 4)   # wide layout for 2 plots

par(mfrow = c(1, 2),         
    mar = c(4,4,3,1))

# Call function for each method
plot_realdata_kde_vs_estimators(YL, YU, "truncation")
plot_realdata_kde_vs_estimators(YL, YU, "oracle")

dev.off()
