library(snowfall)
library(MASS)
library(mvtnorm)
library(dplyr)
library(tidyr)

# 시뮬레이션 설정
n_replicates <- 200    
n_cpus <- parallel::detectCores() - 1 

# 이상치 기준 설정 
ALPHA_MIN <- 0.0001
ALPHA_MAX <- 2
MU_MIN <- -2
MU_MAX <- 8
SIGMA_MIN <- 0.0000001
SIGMA_MAX <- 5  # Naive에서 2000+ 나오는 거 방지

run_simulation_worker <- function(iter) {
  ecount <- 0  
  eflag  <- 1   
  
  while(eflag > 0) {
    
    # seed 설정 (재시도시 ecount만큼 증가)
    set.seed(iter + ecount * 10000)
    
    # sample size
    n <- 1000     # total_sample_size
    n1 <- 200     # P_sample_size
    pi1 <- n1/n
    
    # parameter
    p <- 3                    # dimension of X
    mu.yp <- 0.1              # P(Y=1)
    theta_true <- 0.9         # Q(Y=1) 
    mu1_vec_true <- rep(3, p) # mean of X give Y=1
    sigma_sq_true <- 1        # variance of X give Y
    
    # =======================================================
    
    # P 데이터 생성
    yp <- rbinom(n1, 1, mu.yp)
    xp <- matrix(NA, nrow=n1, ncol=p)
    for(i in 1:n1) {
      if(yp[i] == 1) {
        xp[i, ] <- mvrnorm(1, mu1_vec_true, sigma_sq_true * diag(p)) 
      } else { 
        xp[i, ] <- mvrnorm(1, rep(0, p), sigma_sq_true * diag(p)) 
      }
    }
    P <- data.frame(yp = factor(yp), xp)
    
    # Q 데이터 생성
    yq_true <- rbinom(n - n1, 1, theta_true)
    xq <- matrix(NA, nrow=n-n1, ncol=p)
    for(i in 1:(n-n1)) {
      if(yq_true[i] == 1) {
        xq[i, ] <- mvrnorm(1, mu1_vec_true, sigma_sq_true * diag(p))
      } else {
        xq[i, ] <- mvrnorm(1, rep(0, p), sigma_sq_true * diag(p))
      }
    }
    
    df_Q <- data.frame(yq = factor(yq_true*0), xq)
    colnames(df_Q) <- colnames(P) 
    Q <- df_Q
    
    # Initial Classifier
    fit.lda <- classifier_P(P)
    
    # =======================================================
    
    # 초기값
    start_vals <- c(
      log(theta_true / (1 - theta_true)),
      mu1_vec_true,
      log(sigma_sq_true)
    )
    
    results <- list()
    
    # Efficient 모델 =========================================
    tryCatch({
      estimating_function_eff <- function(params_transformed) {
        p1 <- params_transformed[1] 
        p2 <- params_transformed[2:(p+1)] 
        p3 <- params_transformed[p+2] 
        
        current_alpha <- 1 / (1 + exp(-p1)) 
        current_mu1 <- p2
        current_sigma_sq <- exp(p3)
        
        current_thetas <- list(
          alpha    = current_alpha,
          mu1      = current_mu1,
          mu0      = rep(0, p), 
          sigma_sq = current_sigma_sq
        )
        
        phi_values <- eff.f(model = fit.lda, P = P, Q = Q, pi1 = pi1, 
                            thetas = current_thetas)
        return(colMeans(phi_values))
      }
      
      sol <- optim(start_vals,
                   fn = function(param) sum(estimating_function_eff(param)^2),
                   method = "BFGS",
                   control=list(reltol=1e-5))
      
      final_p <- sol$par
      results$Efficient <- c(
        1 / (1 + exp(-final_p[1])),
        final_p[2:(p+1)],
        exp(final_p[p+2])
      )
    }, error = function(e) {
      results$Efficient <<- rep(NA, 5) 
    })
    
    # Naive 모델 =============================================
    tryCatch({
      estimating_function_naive <- function(params_transformed) {
        p1 <- params_transformed[1] 
        p2 <- params_transformed[2:(p+1)] 
        p3 <- params_transformed[p+2] 
        
        current_alpha <- 1 / (1 + exp(-p1)) 
        current_mu1 <- p2
        current_sigma_sq <- exp(p3)
        
        current_thetas <- list(
          alpha    = current_alpha,
          mu1      = current_mu1,
          mu0      = rep(0, p), 
          sigma_sq = current_sigma_sq
        )
        
        phi_values <- naive.f(model = fit.lda, P = P, Q = Q, pi1 = pi1, 
                              thetas = current_thetas)
        return(colMeans(phi_values)) 
      }
      
      sol <- optim(start_vals,
                   fn = function(param) sum(estimating_function_naive(param)^2),
                   method = "BFGS",
                   control=list(reltol=1e-5))
      
      final_p <- sol$par
      results$Naive <- c(
        1 / (1 + exp(-final_p[1])),
        final_p[2:(p+1)],
        exp(final_p[p+2])
      )
    }, error = function(e) {
      results$Naive <<- rep(NA, 5) 
    })
    
    # BBSE 모델 ==============================================
    tryCatch({
      estimating_function_bbse <- function(params_transformed) {
        p1 <- params_transformed[1] 
        p2 <- params_transformed[2:(p+1)] 
        p3 <- params_transformed[p+2] 
        
        current_alpha <- 1 / (1 + exp(-p1)) 
        current_mu1 <- p2
        current_sigma_sq <- exp(p3)
        
        current_thetas <- list(
          alpha    = current_alpha,
          mu1      = current_mu1,
          mu0      = rep(0, p), 
          sigma_sq = current_sigma_sq
        )
        
        phi_values <- bbse.f(model = fit.lda, P = P, Q = Q, pi1 = pi1, 
                             thetas = current_thetas)
        return(colMeans(phi_values)) 
      }
      
      sol <- optim(start_vals,
                   fn = function(param) sum(estimating_function_bbse(param)^2),
                   method = "BFGS",
                   control=list(reltol=1e-5))
      
      final_p <- sol$par
      results$BBSE <- c(
        1 / (1 + exp(-final_p[1])),
        final_p[2:(p+1)],
        exp(final_p[p+2])
      )
    }, error = function(e) {
      results$BBSE <<- rep(NA, 5) 
    })
    
    # RLLS 모델 ==============================================
    tryCatch({
      estimating_function_rlls <- function(params_transformed) {
        p1 <- params_transformed[1] 
        p2 <- params_transformed[2:(p+1)] 
        p3 <- params_transformed[p+2] 
        
        current_alpha <- 1 / (1 + exp(-p1)) 
        current_mu1 <- p2
        current_sigma_sq <- exp(p3)
        
        current_thetas <- list(
          alpha    = current_alpha,
          mu1      = current_mu1,
          mu0      = rep(0, p), 
          sigma_sq = current_sigma_sq
        )
        
        phi_values <- rlls.f(model = fit.lda, P = P, Q = Q, pi1 = pi1, 
                             thetas = current_thetas)
        return(colMeans(phi_values)) 
      }
      
      sol <- optim(start_vals,
                   fn = function(param) sum(estimating_function_rlls(param)^2),
                   method = "BFGS",
                   control=list(reltol=1e-5))
      
      final_p <- sol$par
      results$RLLS <- c(
        1 / (1 + exp(-final_p[1])),
        final_p[2:(p+1)],
        exp(final_p[p+2])
      )
    }, error = function(e) {
      results$RLLS <<- rep(NA, 5) 
    })
    
    # 결과 정리 ===============================================
    
    result_vec <- c(results$Efficient, results$Naive, results$BBSE, results$RLLS)  
    is_outlier <- FALSE #while 루프 > 매번 초기화 필요  
    
    if(sum(is.na(result_vec)) > 0) { #NA 안전장치
      is_outlier <- TRUE
    } else {
      # 각 모델별로 체크
      for(i in 1:4) {
        model_result <- result_vec[(1 + (i-1)*5):(5 + (i-1)*5)]
        # i=1 > 1:5   (Efficient)
        # i=2 > 6:10  (Naive)
        # i=3 > 11:15 (BBSE)
        # i=4 > 16:20 (RLLS)
        
        if(model_result[1] < ALPHA_MIN || model_result[1] > ALPHA_MAX) {
          is_outlier <- TRUE
          break
        }
        
        # Mu 체크 (2, 3, 4번째 값)
        if(any(model_result[2:4] < MU_MIN) || any(model_result[2:4] > MU_MAX)) {
          is_outlier <- TRUE
          break
        }
        
        # Sigma_sq 체크 (5번째 값)
        if(model_result[5] < SIGMA_MIN || model_result[5] > SIGMA_MAX) {
          is_outlier <- TRUE
          break
        }
      }
    }
    
    if(is_outlier) {
      ecount <- ecount + 1
    } else {
      eflag <- 0
    }
    
  }  # while문 종료
  
  return(c(result_vec, ecount))
}  # 함수 종료

# 병렬 실행 ===============================================

sfInit(parallel = TRUE, cpus = n_cpus)
sfLibrary(MASS)
sfLibrary(mvtnorm)
sfExport("classifier_P", "U.f", "theta.f", "pyx.f", "rho.f", "w.f", "E.f", 
         "a.f", "E_star.f", "eff.f", "naive.f", "CM.f", "Q.f", "rho.bbse.f", 
         "bbse.f", "pred_dist.f", "CM_joint.f", "norm.f", "rho.rlls.f", "rlls.f",
         "ALPHA_MIN", "ALPHA_MAX", "MU_MIN", "MU_MAX", "SIGMA_MIN", "SIGMA_MAX")

results_list <- sfLapply(1:n_replicates, run_simulation_worker)
sfStop()

# 데이터 프레임으로 변환
results_df <- do.call(rbind, results_list) %>% as.data.frame()
param_names <- c("Alpha", "Mu1_1", "Mu1_2", "Mu1_3", "Sigma_Sq")
model_names <- c("Efficient", "Naive", "BBSE", "RLLS")

col_names <- c()
for(model in model_names) {
  col_names <- c(col_names, paste0(model, "_", param_names))
}
col_names <- c(col_names, "ecount")  
colnames(results_df) <- col_names

results_df_analysis <- results_df[, 1:20]

# 결과 요약 ===============================================

# True values
true_values <- data.frame(
  Parameter = param_names,
  True_Value = c(0.9, 3, 3, 3, 1.0)
)

# 각 모델별로 요약 통계 계산
summary_list <- list()

for(model in model_names) {
  model_cols <- paste0(model, "_", param_names)
  
  model_data <- results_df_analysis[, model_cols]
  colnames(model_data) <- param_names
  
  model_long <- model_data %>%
    mutate(Iter = row_number()) %>%
    pivot_longer(cols = -Iter, names_to = "Parameter", values_to = "Estimate")
  
  model_summary <- model_long %>%
    left_join(true_values, by = "Parameter") %>%
    group_by(Parameter) %>%
    summarise(
      Model = model,
      True_Value = mean(True_Value),
      Est_Mean = mean(Estimate, na.rm = TRUE),
      Bias = mean(Estimate, na.rm = TRUE) - mean(True_Value),
      SE = sd(Estimate, na.rm = TRUE),
      MSE = mean((Estimate - True_Value)^2, na.rm = TRUE),
      RMSE = sqrt(mean((Estimate - True_Value)^2, na.rm = TRUE)),
      .groups = 'drop'
    ) %>%
    mutate(across(where(is.numeric), ~ round(., 4)))
  
  summary_list[[model]] <- model_summary
}

# overview
all_summary <- bind_rows(summary_list)
print(all_summary)

#ecount 
ecount_summary <- sum(results_df$ecount)
cat("\n=== ecount===\n")
ecount_summary


# Classification Accuracy Evaluation 

eval_performance <- function(xq_data, yq_actual, est_alpha, est_mu1, est_mu0, est_Sigma, p){
  log_f1_hat <- dmvnorm(xq_data, mean = est_mu1, sigma = est_Sigma, log = TRUE)
  log_f0_hat <- dmvnorm(xq_data, mean = est_mu0, sigma = est_Sigma, log = TRUE)
  
  log_prob1_hat <- log(est_alpha) + log_f1_hat
  log_prob0_hat <- log(1 - est_alpha) + log_f0_hat
  
  y_pred <- ifelse(log_prob1_hat > log_prob0_hat, 1, 0)
  
  TP <- as.numeric(sum(yq_actual == 1 & y_pred == 1))
  TN <- as.numeric(sum(yq_actual == 0 & y_pred == 0))
  FP <- as.numeric(sum(yq_actual == 0 & y_pred == 1))
  FN <- as.numeric(sum(yq_actual == 1 & y_pred == 0))
  
  Pos <- TP + FN
  Neg <- TN + FP
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(denom == 0, 0, (TP * TN - FP * FN) / denom) 
  
  return(c(Pos = Pos, Neg = Neg, TP = TP, TN = TN, FP = FP, FN = FN, 
           Accuracy = accuracy, MCC = mcc))
}

#Q_data set 재생성
set.seed(999)
n_test <- 10000
p_dim  <- 3
theta_true_q <- 0.9
mu1_true_vec <- rep(3, p_dim)
mu0_true_vec <- rep(0, p_dim)
sigma_sq_true_val <- 1

yq_test <- rbinom(n_test, 1, theta_true_q)
xq_test <- matrix(NA, nrow=n_test, ncol=p_dim)
for(i in 1:n_test) {
  if(yq_test[i] == 1) {
    xq_test[i, ] <- mvrnorm(1, mu1_true_vec, sigma_sq_true_val * diag(p_dim))
  } else {
    xq_test[i, ] <- mvrnorm(1, mu0_true_vec, sigma_sq_true_val * diag(p_dim))
  }
}

# base line 
# Oracle
oracle_thetas <- list(
  alpha = theta_true_q, mu1 = mu1_true_vec, 
  mu0 = mu0_true_vec, Sigma = sigma_sq_true_val * diag(p_dim)
)

# Source
source_thetas <- list(
  alpha = 0.1, mu1 = mu1_true_vec, 
  mu0 = mu0_true_vec, Sigma = sigma_sq_true_val * diag(p_dim)
)

# evaluation
performance_list <- list()

# Oracle 및 Source는 고정값이므로 한 번만 평가
res_oracle <- eval_performance(xq_test, yq_test, oracle_thetas$alpha, oracle_thetas$mu1, oracle_thetas$mu0, oracle_thetas$Sigma, p_dim)
performance_list[["Oracle"]] <- data.frame(Model = "Oracle", t(round(res_oracle, 4)))

res_source <- eval_performance(xq_test, yq_test, source_thetas$alpha, source_thetas$mu1, source_thetas$mu0, source_thetas$Sigma, p_dim)
performance_list[["Source"]] <- data.frame(Model = "Source", t(round(res_source, 4)))

# 시뮬레이션으로 추정된 4가지 모델 성능 평가
for (model in model_names) {
  col_alpha <- paste0(model, "_Alpha")
  col_mu    <- paste0(model, "_Mu1_", 1:p_dim)
  col_sigma <- paste0(model, "_Sigma_Sq")
  
  model_metrics <- list()
  
  for (i in 1:nrow(results_df_analysis)) {
    est_alpha <- results_df_analysis[i, col_alpha]
    est_mu1   <- as.numeric(results_df_analysis[i, col_mu])
    est_sigma <- results_df_analysis[i, col_sigma]
    
    if (is.na(est_alpha)) next 
    res <- eval_performance(
      xq_data   = xq_test,
      yq_actual = yq_test,
      est_alpha = est_alpha,
      est_mu1   = est_mu1,
      est_mu0   = rep(0, p_dim), 
      est_Sigma = est_sigma * diag(p_dim),
      p         = p_dim
    )
    model_metrics[[i]] <- res
  }
  
  # 200회 측정된 성능 지표들의 평균 계산
  avg_metrics <- colMeans(do.call(rbind, model_metrics), na.rm = TRUE)
  performance_list[[model]] <- data.frame(Model = model, t(round(avg_metrics, 4)))
}

performance_df <- do.call(rbind, performance_list)
rownames(performance_df) <- NULL

cat("\n=== [3] 분류 예측 성능 (시뮬레이션 200회 평균) ===\n")
performance_df
