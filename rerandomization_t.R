#---- Packages ----
if (!require("pacman")) 
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load("tidyverse", "MASS", "Matrix", "matrixcalc")

#Compute Mahalanobis Distance according to Morgan dissertation (pg 33)
mahal_dist <- function(n, pw, vector, samp_cov){
  return(n*pw*(1 - pw)*t(vector)%*%solve(samp_cov)%*%vector)
}

rerandomization_ttest <- function(samp_size, num_cov, pw, a){
  #Generate Correlation Matrices
  corr_A <- diag(num_cov + 1)
  
  corr_B <- diag(num_cov + 1)
  corr_B[corr_B == 0] <- 0.20
  
  corr_C <- diag(num_cov + 1)
  corr_C[corr_C == 0] <- 0.40
  
  corr_D <- diag(num_cov + 1)
  corr_D[corr_D == 0] <- 0.60
  
  corr_E <- diag(num_cov + 1)
  corr_E[corr_E == 0] <- 0.80
  
  corr_F <- Matrix(runif((num_cov + 1)*(num_cov + 1), min = 0, max = 1), 
                   (num_cov + 1)) %>% 
    forceSymmetric() %>% as.matrix()
  diag(corr_F) <- 1
  corr_F = nearPD(corr_F, keepDiag = TRUE)$mat %>% as.matrix
  
  corr_list <- list(corr_A, corr_B, corr_C, corr_D, corr_E, corr_F)
  
  #Compute variance-covariance matrices
  cov_mat_list <- vector("list", length(corr_list))
  for(i in 1:length(corr_list)){
    corr_mat = corr_list[[i]]
    S <- diag(nrow(corr_mat))               #Diagonal matrix of covariate SDs; assumes standardized covariates
    cov_mat_list[[i]] <- S%*%corr_mat%*%S   #Compute Variance/Covariance matrix 
  }
  
  results <- vector("list", length(cov_mat_list))
  
  for(i in 1:length(cov_mat_list)){
    cov_mat = cov_mat_list[[i]]
    
    #Generate the data
    data <- mvrnorm(n = samp_size, mu = rep(0, num_cov + 1), Sigma = cov_mat)
    samp_cov <- cov(data[, 1:num_cov])
    
    #Find an acceptable rerandomization
    Tx_ind <- sample(seq_len(nrow(data)), size = pw*nrow(data))
    Tx <- data[Tx_ind, ]
    Tc <- data[-Tx_ind, ]
    differences <- (colMeans(Tx)[1:num_cov] - colMeans(Tc)[1:num_cov])
    M = mahal_dist(samp_size, pw, differences, samp_cov)
    
    #Store Mean Differences
    mean_diffs = mean(Tx[, (num_cov + 1)]) - mean(Tc[, (num_cov + 1)])
    
    while(M > a){
      Tx_ind <- sample(seq_len(nrow(data)), size = pw*nrow(data)) 
      Tx <- data[Tx_ind, ]
      Tc <- data[-Tx_ind, ]
      differences <- (colMeans(Tx)[1:num_cov] - colMeans(Tc)[1:num_cov])
      M = mahal_dist(samp_size, pw, differences, samp_cov)
    }
    p_val = t.test(Tx[, (num_cov + 1)], Tc[, (num_cov + 1)])$p.value
    mean_diffs = c(mean_diffs, 
                   mean(Tx[, (num_cov + 1)]) - mean(Tc[, (num_cov + 1)]))
    results[[i]] <- list("p_val" = p_val, "mean_diffs" = mean_diffs)
  }
  return(results)
}
