## tau_delta

tau_from_delta <- function(delta) {
  m <- ncol(delta)

  if (m == 1) return(cbind(-1e10, delta[, 1], 1e10))

  E <- exp(delta[, 2:m, drop = FALSE])
  U <- matrix(0, m - 1, m - 1)
  U[upper.tri(U, diag = TRUE)] <- 1
  tau_fin <- cbind(delta[, 1, drop = FALSE], delta[, 1] + E %*% U)

  cbind(-1e10, tau_fin, 1e10)
}

delta_from_tau <- function(tau) {
  tau_fin <- tau[, 2:(ncol(tau) - 1), drop = FALSE]
  m <- ncol(tau_fin)

  if (m == 1) return(tau_fin)

  gaps <- tau_fin[, 2:m, drop = FALSE] - tau_fin[, 1:(m - 1), drop = FALSE]
  gaps <- pmax(gaps, 1e-6)
  cbind(tau_fin[, 1, drop = FALSE], log(gaps))
}

#####################
#-----Gradients-----#
#####################

grad_tau <- function(Y, J, C, PFf, clip_value) {
  G <- matrix(0, nrow = J, ncol = C - 1)

  for (c in 2:C) {
    fc <- PFf$f_all[, seq(c - 1, J * (C - 1), by = C - 1)]
    den_left  <- PFf$P_all[, seq(c - 1, J * C, by = C)]
    den_right <- PFf$P_all[, seq(c,     J * C, by = C)]

    G[, c - 1] <- colSums(
      (Y == (c - 1)) * (fc / den_left) -
        (Y == c)       * (fc / den_right)
    )
  }

  G <- -G
  G <- pmax(pmin(G, clip_value), -clip_value)

  return(G)
}

grad_delta <- function(grad_t, tau) {
  m <- ncol(grad_t)

  M <- matrix(0, m, m)
  M[lower.tri(M, diag = TRUE)] <- 1
  grad_delta <- grad_t %*% M

  gaps <- tau[, 3:(m + 1), drop = FALSE] - tau[, 2:m, drop = FALSE]
  grad_delta[, 2:m] <- gaps * grad_delta[, 2:m, drop = FALSE]

  return(grad_delta)
}

fisher_tau <- function(J,C,PFf){
  f_all_squared <- PFf$f_all^2
  P_left <- PFf$P_all[,-seq(C,J*C,by=C)]
  P_right <- PFf$P_all[,-seq(1,J*C,by=C)]
  diag_vec <- colSums(f_all_squared/P_left+f_all_squared/P_right)

  f_left  <- PFf$f_all[,-seq(C-1,J*(C-1),by=C-1)]
  f_right <- PFf$f_all[,-seq(1,  J*(C-1),by=C-1)]
  P_mid <- PFf$P_all[,-c(seq(1,J*C,by=C),
                         seq(C,J*C,by=C))]
  off_vec <- -colSums((f_left*f_right)/P_mid)

  return(list(diag_vec=diag_vec,off_vec=off_vec))
}

solve_blocks <- function(J,C,fisher_tau_list,grad_t,tau){
  diag_vec <- fisher_tau_list$diag_vec
  off_vec <- fisher_tau_list$off_vec
  I_block <- matrix(0, nrow=C-1,ncol=C-1)
  x_mat <- matrix(NA, nrow=nrow(grad_t), ncol=ncol(grad_t))

  for (i in 1:J){
    grad_t_vec <- grad_t[i,]
    tau_vec <- tau[i,]
    diag(I_block) <- diag_vec[(1+(i-1)*(C-1)):(i*(C-1))]+ 1e-6
    I_block[row(I_block)+1==col(I_block)] <- off_vec[(1+(i-1)*(C-2)):(i*C-2*i)]
    I_block[col(I_block)+1==row(I_block)] <- I_block[row(I_block)+1==col(I_block)]

    y <- solve(I_block,grad_t_vec)
    Ly <- c(y[1], diff(y))
    gaps <- tau_vec[3:C]-tau_vec[2:(C-1)]
    x_mat[i,] <- c(Ly[1], Ly[2:(C-1)]/gaps)
  }

  return(x_mat)
}

# X (N x J linear predictor) replaces the old PFf argument: the gradient
# and Fisher blocks now come from the fused tau_stats_grm() kernel, which
# streams over rows instead of materializing P_all / f_all and the
# full-width temporaries grad_tau()/fisher_tau() used to build. Results are
# bitwise identical to the old grad_tau() + fisher_tau() pair.
update_delta <- function(Params, Y, J, C, alpha, X, obj_prev, obj_fn,
                         max_halve = 10) {
  delta_old <- Params$delta
  st            <- tau_stats_grm(Y, X, Params$tau, clip_value = 0.1 * nrow(Y))
  grad_t        <- st$G
  fisher_tau_list <- list(diag_vec = st$diag_vec, off_vec = st$off_vec)
  x_mat         <- solve_blocks(J, C, fisher_tau_list, grad_t, Params$tau)
  Params$delta  <- Params$delta - alpha * x_mat
  Params$tau    <- tau_from_delta(Params$delta)
  obj           <- obj_fn(Params)
  step_halve_num <- 0
  while (!is.finite(obj) || isTRUE(obj > obj_prev)) {
    Params$delta <- (delta_old + Params$delta) / 2
    Params$tau   <- tau_from_delta(Params$delta)
    obj          <- obj_fn(Params)
    step_halve_num <- step_halve_num + 1
    if (step_halve_num == max_halve) break
  }
  list(Params = Params, obj = obj)
}

# Convergence tracker with mutable state (environment-based).
# Usage:
#   tracker <- convergence_tracker(max_iter, init_value, label = "Theta+A",
#                                    verbose = TRUE, print_init = TRUE)
#   for (i in 1:max_iter) {
#     # ... compute new_value ...
#     if (check_converged(tracker, i, new_value, tol)) break
#   }
#   obj_vals <- tracker$obj_vals  # retrieve full history if needed
convergence_tracker <- function(max_iter, init_value, label = NULL,
                                  verbose = TRUE, print_init = FALSE) {
  env <- new.env(parent = emptyenv())
  env$obj_vals <- rep(NA, max_iter + 1)
  env$obj_vals[1] <- init_value
  env$label <- label
  env$verbose <- verbose
  if (verbose && print_init) {
    print(sprintf("The initial value is %.2f", init_value))
  }
  env
}

check_converged <- function(tracker, i, new_value, tol) {
  tracker$obj_vals[i + 1] <- new_value
  if (tracker$verbose) {
    lbl <- if (is.null(tracker$label)) "" else sprintf(" (%s)", tracker$label)
    print(sprintf("After %.d update%s value is %.2f", i, lbl, new_value))
  }
  abs(new_value - tracker$obj_vals[i]) / tracker$obj_vals[i] < tol
}


Linv <- function(u, m) {
  stopifnot(length(u) %% m == 0)
  J <- length(u) / m

  U <- matrix(u, nrow = m, ncol = J)
  B <- rbind(U[1, , drop=FALSE], U[-1, , drop=FALSE] - U[-m, , drop=FALSE])
  as.vector(B)
}

#############################
#-----Numerical helpers-----#
#############################

# solve() with a pseudo-inverse fallback for exactly singular systems --
# e.g. crossprod(A) when a very sparse fit zeroes out an entire factor
# column of A ("Lapack routine dgesv: system is exactly singular").
# Identical to solve() whenever the system is nonsingular; in the singular
# case the minimum-norm solution leaves dead factors at zero.
safe_solve <- function(M, b = NULL) {
  if (is.null(b)) {
    tryCatch(solve(M), error = function(e) MASS::ginv(M))
  } else {
    tryCatch(solve(M, b), error = function(e) MASS::ginv(M) %*% b)
  }
}

###########################
#-----Standardization-----#
###########################

standardize <- function(Theta,A,tau=NULL,C=NULL,
                        model=c("GRM","GPCM")){
  model <- match.arg(model, c("GRM","GPCM"))

  if (model=="GRM") {
    Theta_sd <- apply(Theta, 2, sd)
    Theta_sd[Theta_sd < 1e-12] <- 1   # dead factor: leave column untouched
    Theta <- t(t(Theta)/Theta_sd)
    A <- t(t(A)*Theta_sd)

    Theta_mean <- colMeans(Theta)
    sd_c <- apply(Theta, 2, sd)
    sd_c[sd_c < 1e-12] <- 1
    Theta <- t((t(Theta)-Theta_mean)/sd_c)
    tau[,2:C] <- tau[,2:C]-drop(A%*%Theta_mean)

    return(list(Theta=Theta,A=A,tau=tau))
  } else {
    Theta_sd <- apply(Theta, 2, sd)
    Theta_sd[Theta_sd < 1e-12] <- 1   # dead factor: leave column untouched
    Theta <- sweep(Theta, 2, Theta_sd, "/")
    A <- sweep(A, 2, Theta_sd, "*")

    Theta_mean <- colMeans(Theta)
    Theta <- sweep(Theta, 2, Theta_mean, "-")

    c <- 0:(C - 1)
    v <- as.vector(A %*% Theta_mean)
    Af_mu <- rep(v, each = C) * rep(c, times = length(v))
    # Paper convention X = (c-1)*eta + tau: centering shift on Theta adds
    # (c-1)*v_j to tau (sign opposite the old X = (c-1)*eta - tau code).
    tau <- tau + Af_mu

    return(list(Theta=Theta,A=A,tau=tau))
  }
}

##################################
#-----Shared GRM update steps-----#
##################################

# MM working response for the GRM: Z = X - 2H, where H is the majorizer
# gradient term. Collapses the PFf / compute_H_grm / Z three-liner that
# was repeated in every GRM fitter.
#
# Dispatches to the fused compute_Z_grm_cpp() kernel, which streams over
# rows and never allocates the N x J*C P_all array or the N x J H / 2*H
# temporaries. Bitwise identical to the former three-step version; C is
# taken from tau and kept in the signature for backward compatibility.
compute_Z_grm <- function(Y, X, tau, C = ncol(tau) - 1) {
  compute_Z_grm_cpp(Y, X, tau)
}

# Theta update shared by GRM_lasso_HiMIRT / GRM_CHiMIRT / GRM_card_HiMIRT.
# orthogonal = TRUE : Procrustes (SVD) step, then center+scale Theta and
#                     absorb the mean shift into tau.
# orthogonal = FALSE: least-squares step, then standardize().
# use_tilde: use A_tilde (ADMM copy) as the loading matrix in the Theta
#            step; the standardized loadings are written back to both
#            A and A_tilde (lasso behavior).
update_Theta_grm <- function(Params, Z, N, C, orthogonal, use_tilde = FALSE) {
  A_load <- if (use_tilde) Params$A_tilde else Params$A

  if (orthogonal) {
    svd_result <- svd(1 / sqrt(N) * Z %*% A_load)
    Params$Theta <- sqrt(N) * svd_result$u %*% t(svd_result$v)

    Theta_mean <- colMeans(Params$Theta)
    Params$Theta <- scale(Params$Theta)
    Params$tau[, 2:C] <- Params$tau[, 2:C] - drop(Params$A %*% Theta_mean)
  } else {
    Theta <- Z %*% A_load %*% safe_solve(crossprod(A_load))

    std <- standardize(Theta, A_load, tau = Params$tau, C = C, model = "GRM")
    Params$Theta <- std$Theta
    Params$A     <- std$A
    Params$tau   <- std$tau
    if (use_tilde) Params$A_tilde <- std$A
  }
  Params$delta <- delta_from_tau(Params$tau)

  Params
}

###################################
#-----Shared GPCM update steps-----#
###################################

# One MM step for tau. The full C-vector step under the 1/2 * E bound followed
# by the tau_1 = 0 anchoring is identically the Bohning free-block step, since
# the multinomial score sums to zero over categories.
update_tau_gpcm <- function(Y, eta, tau, N, J, C) {
  H_res <- compute_H_gpcm(Y, eta, tau, C)
  tau <- tau + 2/N * H_res$col_sums
  tau_mat <- matrix(tau, nrow = C, ncol = J)
  tau_mat <- sweep(tau_mat, 2, tau_mat[1, ], "-")
  as.vector(tau_mat)
}

# MM working response for the GPCM: Z = eta + Hw / w, where w = c' B c is the
# majorization constant (Bohning 1992, Lemma 2.2 restricted to the GPCM
# direction c = (0,...,C-1)). Returns Hw and eta as well: the card MM subloop
# needs both, and Hw alone gives the log-likelihood gradient wrt A as
# -crossprod(Hw, Theta).
compute_Z_gpcm <- function(Y, eta, tau, C, w) {
  H_res <- compute_H_gpcm(Y, eta, tau, C)
  list(Z = eta+H_res$Hw/w, Hw = H_res$Hw, eta = eta)
}

# Theta update shared by GPCM_lasso_HiMIRT / GPCM_CHiMIRT / GPCM_card_HiMIRT.
# Mirrors update_Theta_grm: Procrustes (orthogonal) vs least squares,
# followed by standardize(). use_tilde: the ADMM copy A_tilde is the
# loading matrix; the standardized loadings go back to A_tilde only.
update_Theta_gpcm <- function(Params, Z, N, C, orthogonal,
                              use_tilde = FALSE) {
  A_load <- if (use_tilde) Params$A_tilde else Params$A

  if (orthogonal) {
    svd_result <- svd(1 / sqrt(N) * Z %*% A_load)
    Theta <- sqrt(N) * svd_result$u %*% t(svd_result$v)
  } else {
    Theta <- Z %*% A_load %*% safe_solve(crossprod(A_load))
  }

  std <- standardize(Theta, A_load, tau = Params$tau, C = C, model = "GPCM")
  Params$Theta <- std$Theta
  if (use_tilde) Params$A_tilde <- std$A else Params$A <- std$A
  Params$tau <- std$tau

  Params
}

proj_K <- function(Params, C, K, item=F, person=F){
  M <- 5*sqrt(K)
  Params$proj=F

  if (person) {
    Theta <- Params$Theta
    vec_1 <- which(rowSums(Theta^2) + 1 > M^2)
    if (length(vec_1) != 0){
      s1 <- sqrt((M^2 - 1) / rowSums(Theta[vec_1, , drop=FALSE]^2))
      Theta[vec_1, ] <- s1 * Theta[vec_1, , drop=FALSE]
      Params$Theta <- Theta
      Params$proj=T
    }
  }

  if (item){
    tau_fin <- Params$tau[,2:C]
    Atau <- cbind(sqrt(C-1)*Params$A, tau_fin)
    vec_2 <- which(rowSums(Atau^2) > (C-1)*M^2)
    if (length(vec_2) != 0){
      s2 <- sqrt(C-1) * M / sqrt(rowSums(Atau[vec_2, , drop=FALSE]^2))
      Atau[vec_2, ] <- s2 * Atau[vec_2, , drop=FALSE]
      Params$A <- Atau[, 1:K] / sqrt(C-1)
      Params$tau[, 2:C] <- Atau[, (K+1):ncol(Atau)]
      Params$delta <- delta_from_tau(Params$tau)
      Params$proj=T
    }
  }

  return(Params)
}

################
#-----ADMM-----#
################

soft_threshold <- function(lambda, rho, x){
  x_abs <- abs(x)
  shrinkage_factor <- pmax(0, 1 - (lambda / (rho * x_abs)))
  return(shrinkage_factor*x)
}

getResiduals <- function(A,A_tilde,A_old,U){
  r=sum((A-A_tilde)^2)/sum(A^2)
  s=sum((A-A_old)^2)/sum(U^2)

  return(c(r,s))
}

subADMM <- function(Y,Z,A,A_tilde,Theta,C=NULL,
                    lambda,model=c("GRM","GPCM"),N_iter){
  model <- match.arg(model, c("GRM","GPCM"))
  N <- nrow(Y);J <- ncol(Y);K <- ncol(Theta)
  U=matrix(0,nrow=J,ncol=K)

  if (model=="GRM") {
    G=crossprod(Theta)/2
    rho=sum(diag(G))/K
    L=t(chol(G+rho*diag(1,nrow(G)),pivot=F))
    M=crossprod(Z,Theta)/2
  } else {
    w <- C*(C^2-1)/24
    G <- crossprod(Theta)*w
    rho=sum(diag(G))/K
    L=t(chol(G+rho*diag(1,nrow(G)),pivot=F))
    M <- crossprod(Z,Theta)*w
  }

  relative_res=c(1,1)

  for (i in 1:N_iter){
    A_tilde <- t(backsolve(t(L), forwardsolve(L, t(M+rho*(A+U)))))

    A_old=A
    A=soft_threshold(lambda,rho,A_tilde-U)

    U=U+A-A_tilde

    relative_res=getResiduals(A,A_tilde,A_old,U)

    if (all(relative_res<1e-4)){
      break
    }
  }

  return(list(A=A,A_tilde=A_tilde))
}

#######################
#-----Cardinality-----#
#######################

cardize <- function(B,kappaC){
  if (length(kappaC)==1){
    B_rank=order(abs(B))
    B[B_rank[1:kappaC]] <- 0
  }else{
    for (i in 1:K){
      B_rank <- order(abs(B[,i]))
      B[B_rank[1:kappaC[i]],i] <- 0
    }
  }

  return(B)
}

# A_update_constrained -------------------------------------------------------
#
# Constrained closed-form A update for the GPCM MM step:
#
#   minimize ||Z - Theta * A^T||_F^2   subject to  A * (1 - Lambda) = 0
#
# The majorization constant w scales the whole objective and so cancels out of
# the argmin; it does not appear here. Z is the working response eta + Hw/w.
#
# When Lambda is a scalar 1 (no mask) or a matrix of all 1s, returns the
# unconstrained global solution `crossprod(Z, Theta) %*% solve(crossprod(Theta))`.
# When Lambda has zeros, solves a per-item least squares over each row's
# active columns of Theta. This is the correct constrained MM minimizer for
# the non-orth case (where the closed-form * Lambda projection breaks the
# monotone-decrease guarantee).
#
# Args:
#   Z:      N x J working response.
#   Theta:  N x K predictor matrix.
#   Lambda: J x K binary mask (or scalar 1 for dense case).
#
# Returns: J x K matrix A with A[j, k] = 0 wherever Lambda[j, k] = 0.
A_update_constrained <- function(Z, Theta, Lambda) {
  J <- ncol(Z); K <- ncol(Theta)

  # Dense case: use the global closed-form solve (fastest path).
  if (length(Lambda) == 1 || all(Lambda == 1)) {
    return(crossprod(Z, Theta) %*% safe_solve(crossprod(Theta)))
  }

  # Sparse case: per-item LS over each item's active columns of Theta.
  A <- matrix(0, J, K)
  for (j in seq_len(J)) {
    act <- which(Lambda[j, ] == 1)
    if (length(act) == 0) next
    Th_act <- Theta[, act, drop = FALSE]
    A[j, act] <- as.vector(
      safe_solve(crossprod(Th_act), crossprod(Th_act, Z[, j]))
    )
  }
  A
}


### Initialization

paramsINIT_binary <- function(Y,cut_off,K,epsilon= .01){
  N=nrow(Y)
  Y <- 1 * (Y >= cut_off)
  p=1  # the proportion of observed responses
  Params_start <- list()

  Q <- 2*Y-1
  Q.SVD <- svd(Q)
  U.Q <- Q.SVD$u
  V.Q <- Q.SVD$v
  Sigma.Q <- diag(Q.SVD$d)
  J <- sum(Q.SVD$d>=2*sqrt(N*p))

  Q_tilde <- U.Q[,1:J,drop=F]%*%Sigma.Q[1:J,1:J,drop=F]%*%t(V.Q[,1:J,drop=F])

  # SVD 1 outputs are dead once Q_tilde exists; drop them so the transform
  # below can reuse the space instead of growing the heap.
  rm(Q.SVD, U.Q, V.Q, Sigma.Q, Q, Y)

  # Fused replacement for the nested ifelse() + apply(M,2,sum) + t(t(M)-d)
  # block: one N x J allocation instead of ~11. Bitwise identical.
  init_res <- init_M_tilde(Q_tilde, epsilon)
  Params_start$d <- init_res$d
  M_tilde        <- init_res$M_tilde
  rm(init_res, Q_tilde)

  M.SVD <- svd(M_tilde)
  U.M <- M.SVD$u
  V.M <- M.SVD$v
  Sigma.M <- diag(M.SVD$d)

  Params_start$Theta <- U.M[,1:K,drop=F]*sqrt(N)
  Params_start$A <- V.M[,1:K,drop=F]%*%Sigma.M[1:K,1:K,drop=F]/sqrt(N)

  return(Params_start)
}

paramsINIT <- function(Y,C,K,epsilon=.01,model=c("GRM","GPCM")){
  J <- ncol(Y)
  tau_fin <- matrix(NA, nrow=J,ncol=C-1)

  for (cut_off in 2:C){
    Params_c <- paramsINIT_binary(Y,cut_off,K,epsilon=.01)
    tau_fin[, cut_off-1] <- Params_c$d
  }

  cut_off <- C/2
  Params_mid <- paramsINIT_binary(Y,cut_off,K,epsilon=.01)
  Theta <- Params_mid$Theta
  A <- Params_mid$A

  if (model=="GRM") {
    tau <- cbind(-1e10,-tau_fin,1e10) # tau=-d
    delta <- delta_from_tau(tau)
    # Each cut point above is fitted by an INDEPENDENT binary model, so the
    # resulting thresholds are not guaranteed to be ordered -- about one item
    # in fifty comes back with tau_{c+1} < tau_c. A crossed pair gives the
    # category between them a near-zero (or negative) probability, so the
    # f^2 / P terms in fisher_tau() blow up to ~1e16 and solve_blocks() hits
    # "system is computationally singular". delta_from_tau() already floors the
    # gaps at 1e-6; regenerating tau from delta propagates that repair, so what
    # leaves this function is always ordered and consistent with delta.
    tau <- tau_from_delta(delta)

    return(list(Theta=Theta,A=A,delta=delta,tau=tau))
  } else {
    tau <- c(matrixStats::colCumsums(rbind(0, t(tau_fin))))

    return(list(Theta=Theta,A=A,tau=tau))
  }
}
