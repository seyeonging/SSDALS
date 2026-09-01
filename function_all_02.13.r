library(MASS)
library(mvtnorm)

# working classifier ---------------------------------------------------------

classifier_P <- function(P) {
  fit.lda <- lda(yp ~ ., data = P) 
  return(fit.lda)
}

# score function -------------------------------------------------------------
 
U.f <- function(y, x, thetas) {
  x <- as.matrix(x)
  alpha <- as.numeric(thetas$alpha)
  mu1   <- as.numeric(thetas$mu1) 
  mu0   <- as.numeric(thetas$mu0) 
  sig   <- as.numeric(thetas$sigma_sq)
  
  Sigma_mat <- sig * diag(length(mu1))
  
  f1 <- dmvnorm(x, mean = mu1, sigma = Sigma_mat)
  f0 <- dmvnorm(x, mean = mu0, sigma = Sigma_mat)
  
  denom <- (alpha*f1 + (1-alpha)*f0)
  q1 <- alpha*f1 / denom
  
  u1 <- (y - q1) / (alpha * (1 - alpha)) 
  x_minus_mu <- t(t(x) - mu1)  
  u2 <- (y - q1) * x_minus_mu / sig
  u3 <- (y - q1) * (sum(mu1^2) - 2*(x %*% mu1)) / (2*sig^2) 
  
  result <- cbind(u1, u2, u3)
  return(as.data.frame(result))
}

# ============================================================================
# eff model
# ============================================================================

theta.f <- function(model, P) {
  p <- ncol(P) - 1
  alpha_hat <- model$prior[2]
  mu1_hat <- model$means["1", ] 
  mu0_hat <- rep(0, p)
  
  x_data <- P[, -1] 
  y_data <- P$yp
  
  x_group0 <- x_data[y_data == "0", ] 
  x_group1 <- x_data[y_data == "1", ] 
  
  SS_0 <- 0 
  SS_1 <- 0
  
  for(j in 1:p) {
    xj_0 <- x_group0[, j] 
    muj_0 <- mu0_hat[j]   
    ss_part_0 <- sum((xj_0 - muj_0)^2)
    SS_0 <- SS_0 + ss_part_0
    
    xj_1 <- x_group1[, j]
    muj_1 <- mu1_hat[j]
    ss_part_1 <- sum((xj_1 - muj_1)^2)
    SS_1 <- SS_1 + ss_part_1
  }
  
  SS_total <- SS_0 + SS_1
  sigma_sq_hat <- SS_total / (nrow(P)*p - 2*p) 
  
  thetas_list <- list(
    alpha    = alpha_hat,
    mu1      = mu1_hat,   
    mu0      = mu0_hat,   
    sigma_sq = sigma_sq_hat
  )
  
  return(thetas_list)
}

pyx.f <- function(model, Q) {
  new_data <- Q[, -1, drop=FALSE]
  pred_p <- predict(model, newdata=new_data)$posterior
  p1 <- pred_p[,"1"] 
  return(data.frame("0"= 1-p1, "1"=p1))
}

# working importance weight  -------------------------------------------------

rho.f <- function(model, P, Q) {
  p_hat <- mean(P$yp == "1")
  q_pred_prob <- pyx.f(model, Q)[,2]
  q_hat <- mean(q_pred_prob)
  return(c(`0` = (1 - q_hat)/(1 - p_hat), `1` =  q_hat/p_hat))
}

# conditional expectation term -----------------------------------------------

w.f <- function(model, Q, P, pi1) {
  X_all <- data.frame(rbind(P[, -1], Q[, -1]))
  A <- data.frame(y_dummy = NA, X_all)  
  
  pyx_prob <- pyx.f(model, A)
  rho_0 <- rho.f(model, P, Q)
  ratio <- pi1 / (1 - pi1)
  
  w0 <- (rho_0[1]^2 + (ratio * rho_0[1]^1)) * pyx_prob[, 1]
  w1 <- (rho_0[2]^2 + (ratio * rho_0[2]^1)) * pyx_prob[, 2]
  w <- 1/ (w0+w1)
  return(w)
}

E.f <- function(model, P, Q, thetas) {
  X_all <- data.frame(rbind(P[, -1], Q[, -1]))
  A <- data.frame(y_dummy = NA, X_all)  
  
  rho_vec <- rho.f(model, P, Q)
  rho_0 <- rho_vec[1]
  rho_1 <- rho_vec[2]
  
  pyx_probs <- pyx.f(model, A)
  p_y0.x <- pyx_probs[,1]
  p_y1.x <- pyx_probs[,2]
  
  U_y0 <- U.f(y = 0, x = X_all, thetas) 
  U_y1 <- U.f(y = 1, x = X_all, thetas)
  
  Term0 <- (rho_0^2 * U_y0) * p_y0.x
  Term1 <- (rho_1^2 * U_y1) * p_y1.x
  
  E_values <- Term0 + Term1 
  return(E_values)
}

a.f <- function(model, P, Q, pi1, thetas) {
  yp_numeric <- as.numeric(as.character(P$yp))
  
  n0 <- sum(yp_numeric==0)
  n1 <- sum(yp_numeric==1)
  np <- sum(n0+n1)
  
  X_all <- data.frame(rbind(P[, -1], Q[, -1]))
  A <- data.frame(y_dummy = NA, X_all)
  w_all <- w.f(model, Q, P, pi1)
  E_all <- E.f(model, P, Q, thetas)
  pxy_all <- pyx.f(model, A)
  rho_vec <- rho.f(model, P, Q)
  
  u_P <- U.f(yp_numeric, P[, -1], thetas)
  w_P <- w_all[1:np]
  pxy_all_P <- pxy_all[1:np, ]
  E_P <- E_all[1:np,]
  
  M_00 <- mean((w_P[yp_numeric==0]) * (rho_vec[1]) * (pxy_all_P[yp_numeric==0,1]))
  M_01 <- mean((w_P[yp_numeric==0]) * (rho_vec[2]) * (pxy_all_P[yp_numeric==0,2]))
  M_10 <- mean((w_P[yp_numeric==1]) * (rho_vec[1]) * (pxy_all_P[yp_numeric==1,1]))
  M_11 <- mean((w_P[yp_numeric==1]) * (rho_vec[2]) * (pxy_all_P[yp_numeric==1,2]))
  M <- matrix(c(M_00, M_01, M_10, M_11), nrow=2, byrow=TRUE)
  
  R_inside <- u_P - (w_P * E_P)
  R_row1 <- colMeans(R_inside[yp_numeric==0, , drop=FALSE], na.rm=TRUE)
  R_row2 <- colMeans(R_inside[yp_numeric==1, , drop=FALSE], na.rm=TRUE)
  
  R <- rbind(R_row1, R_row2)
  A_mat <- solve(M, R)
  
  return(A_mat)
}

E_star.f <- function(model, P, Q, thetas, pi1) {
  X_all <- data.frame(rbind(P[, -1], Q[, -1]))
  A <- data.frame(y_dummy = NA, X_all)
  
  rho_vec <- rho.f(model, P, Q)
  rho_0 <- rho_vec[1]
  rho_1 <- rho_vec[2]
  
  pyx_probs <- pyx.f(model, A)
  p_y0.x <- pyx_probs[,1]
  p_y1.x <- pyx_probs[,2]
  
  a <- a.f(model, P, Q, pi1, thetas)
  a_0 <- a[1, ]
  a_1 <- a[2, ]
  
  U_y0 <- U.f(y = 0, x = X_all, thetas) 
  U_y1 <- U.f(y = 1, x = X_all, thetas)
  
  Term0 <- (rho_0^2 * U_y0) + rho_0 * matrix(rep(a_0, nrow(U_y0)), nrow=nrow(U_y0), byrow=TRUE)
  Term0 <- Term0 * p_y0.x
  
  Term1 <- (rho_1^2 * U_y1) + rho_1 * matrix(rep(a_1, nrow(U_y1)), nrow=nrow(U_y1), byrow=TRUE)
  Term1 <- Term1 * p_y1.x
  
  E_values <- Term0 + Term1 
  return(E_values)
}

# working estimating function : influence function ---------------------------

eff.f <- function(model, P, Q, pi1, thetas) {
  n1 <- nrow(P)
  n <- n1 + nrow(Q)
  
  w_star_values <- w.f(model, Q, P, pi1)  
  E_star_values <- E_star.f(model, P, Q, thetas, pi1) 
  rho_vec <- rho.f(model, P, Q) 
  
  w_P <- as.vector(w_star_values[1:n1])
  E_P <- as.matrix(E_star_values[1:n1, ])
  rho_P <- ifelse(P$yp == "0", rho_vec[1], rho_vec[2]) 
  
  yp <- as.numeric(as.character(P$yp))
  U_P <- as.matrix(U.f(yp, P[, -1], thetas))
  
  term_P_in <- U_P - (w_P * E_P)
  phi_P <- (1/pi1) * (rho_P * term_P_in)
  
  w_Q <- as.vector(w_star_values[(n1 + 1):n])
  E_Q <- as.matrix(E_star_values[(n1 + 1):n, ])
  phi_Q <- (1/(1-pi1)) * (w_Q * E_Q)
  
  phi_all <- rbind(phi_P, phi_Q)
  
  return(phi_all)
}

# ============================================================================
# naive model
# ============================================================================

# working estimating function : naive function -------------------------------

naive.f <- function(model, P, Q, pi1, thetas) {
  rho_vec <- rho.f(model, P, Q) 
  rho_P <- ifelse(P$yp == "0", rho_vec[1], rho_vec[2]) 
  
  yp <- as.numeric(as.character(P$yp))
  U_P <- as.matrix(U.f(yp, P[, -1], thetas))
  
  naive_values <- (1/pi1)*(rho_P * U_P)
  return(naive_values)
}

# ============================================================================
# bbse model
# ============================================================================

# confunsion matrix function -------------------------------------------------

CM.f <- function(model, P) {
  pred_class <- predict(model, P)$class
  class_matrix <- table(
    factor(pred_class, levels = c("0","1")),
    factor(P$yp,       levels = c("0","1")))
  C_hat <- prop.table(class_matrix, margin = 2)
  return(C_hat)
}

# 생략가능
Q.f <- function(model, Q) {
  new_data <- Q[, -1, drop = FALSE]
  post <- predict(model, newdata = new_data)$posterior
  
  qhat1 <- mean(post[, "1"])
  qhat0 <- 1 - qhat1
  return(c(qhat0, qhat1))
}

rho.bbse.f <- function(model, P, Q) {
  C_hat <- as.matrix(CM.f(model, P))
  
  #pred_Q_prob <- Q.f(model, Q)
  pred_Q_prob <- colMeans(pyx.f(model, Q)) 

  theta_hat <- solve(C_hat) %*% pred_Q_prob
  q_hat <- theta_hat[2]
  p_hat <- mean(P$yp == "1")
  return(c(`0` = (1 - q_hat)/(1 - p_hat), `1` =  q_hat/p_hat))
}

# working estimating function : bbse function --------------------------------

bbse.f <- function(model, P, Q, pi1, thetas) {
  rho_bbse_vec <- rho.bbse.f(model, P, Q) 
  rho_bbse <- ifelse(P$yp == "0", rho_bbse_vec[1], rho_bbse_vec[2]) 
  
  yp <- as.numeric(as.character(P$yp))
  U_P <- as.matrix(U.f(yp, P[, -1], thetas))
  
  bbse_values <- (1/pi1)*(rho_bbse * U_P)
  return(bbse_values)
}

# ============================================================================
# rlls model
# ============================================================================

#생략가능
pred_dist.f <- function(model, P) {
  new_data <- P[, -1, drop = FALSE]
  post <- predict(model, newdata = new_data)$posterior
  
  hat1 <- mean(post[, "1"])
  hat0 <- 1 - hat1
  return(c(hat0, hat1))
}

# confusion_matrix_joint version ---------------------------------------------

CM_joint.f <- function(model, P) {
  new_data <- P[, -1, drop = FALSE]
  pred_class <- predict(model, newdata = new_data)$class
  M <- table(
    factor(pred_class, levels = c("0","1")),
    factor(P$yp,       levels = c("0","1"))
  )
  as.matrix(M) / nrow(P)
}

norm.f <- function(x) sqrt(sum(x^2))

# RLLS rho estimator: w = 1 + theta_hat --------------------------------------

rho.rlls.f <- function(model, P, Q, delta = 0.05, alpha_rlls = 0.01) {
  np <- nrow(P)
  K <- 2
  
  CM_joint <- CM_joint.f(model, P)
  
  mu_p_hat <- colMeans(pyx.f(model, P))
  mu_q_hat <- colMeans(pyx.f(model, Q))
  
  #mu_p_hat <- pred_dist.f(model, P) 
  #mu_q_hat <- pred_dist.f(model, Q)
  
  b <- mu_q_hat - mu_p_hat
  
  rho <- 3 * (2*log(2*K/delta)/(3*np) + sqrt(2*log(2*K/delta)/np))
  lambda <- alpha_rlls * rho
  
  opt <- optim(
    par = rep(0, K),
    fn = function(theta) norm.f(CM_joint %*% theta - b) + lambda * norm.f(theta),
    method = "L-BFGS-B",
    lower = rep(-1, K)
  )
  
  theta_hat <- opt$par
  rho_rlls_vec <- 1 + theta_hat
  names(rho_rlls_vec) <- c("0","1")
  return(rho_rlls_vec)
}

# working estimating function : rlls function --------------------------------

rlls.f <- function(model, P, Q, pi1, thetas) {
  rho_rlls_vec <- rho.rlls.f(model, P, Q, delta = 0.05, alpha_rlls = 0.01)
  rho_rlls <- ifelse(P$yp == "0", rho_rlls_vec[1], rho_rlls_vec[2])
  
  yp <- as.numeric(as.character(P$yp))
  U_P <- as.matrix(U.f(yp, P[, -1], thetas))
  rlls_values <- (1/pi1) * (rho_rlls * U_P)
  return(rlls_values)
  }