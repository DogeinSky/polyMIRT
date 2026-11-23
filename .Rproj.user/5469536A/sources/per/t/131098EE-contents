#################################
#-----Graded Response Model-----#
#################################

GRM_HiMMIRT <- function(Y, K, C, alpha = 1, rotation = c(TRUE, FALSE),
                         max_iter, tol = 1e-4) {
  N <- nrow(Y); J <- ncol(Y)

  Params <- paramsINIT(Y, C, K, model = "GRM")

  X <- with(Params, tcrossprod(Theta, A))
  tracker <- convergence_tracker(max_iter, logl_grm(Y, Params$tau, X, NULL),
                                   label = "Theta+A", verbose = TRUE, print_init = FALSE)

  for (i in 1:max_iter) {
    ud <- update_delta(Params, Y, J, C, alpha, X, tracker$obj_vals[i],
                       function(p) logl_grm(Y, p$tau, X, NULL))
    Params <- ud$Params

    Params <- proj_K(Params, C, K, item = TRUE)
    if (Params$proj) {
      X <- with(Params, tcrossprod(Theta, A))
    }

    Z <- compute_Z_grm(Y, X, Params$tau, C)

    Params$Theta <- with(Params, Z %*% A %*% safe_solve(crossprod(A)))
    Params$A     <- with(Params, t(Z) %*% Theta %*% safe_solve(crossprod(Theta)))

    Params <- proj_K(Params, C, K, person = TRUE, item = TRUE)

    X <- with(Params, tcrossprod(Theta, A))
    if (check_converged(tracker, i, logl_grm(Y, Params$tau, X, NULL), tol)) break
  }

  # Last recorded objective (obj_vals is NA-padded past convergence). Named
  # $obj rather than $value to keep it distinct from the sparse fitters'
  # $value, which includes a penalty term; this one is the unpenalized
  # negative log-likelihood. Held in a local because standardize() below
  # returns a FRESH list and would drop any field attached to Params here;
  # it is re-attached just before the return. standardize/SVD/rotation all
  # leave Theta %*% t(A) and the implied probabilities unchanged, so the
  # value stays valid for the returned parameters.
  obj_hist  <- tracker$obj_vals[!is.na(tracker$obj_vals)]
  obj_final <- obj_hist[length(obj_hist)]

  Params <- with(Params, standardize(Theta, A, tau = tau, C = C, model = "GRM"))
  Params$delta <- delta_from_tau(Params$tau)

  svd_result <- svd(Params$Theta)
  Params$Theta <- sqrt(N) * svd_result$u
  Params$A     <- Params$A %*% svd_result$v %*% diag(svd_result$d) / sqrt(N)

  if (rotation) {
    rotate_list <- GPArotation::GPFoblq(Params$A, method = "geomin")
    Params$A    <- rotate_list$loadings
    rotation_M  <- rotate_list$Th
    Params$Theta <- Params$Theta %*% rotation_M
  }

  Params$obj <- obj_final

  return(Params)
}

# Core MM loop shared by GRM_lasso_HiMIRT / GRM_CHiMIRT / GRM_card_HiMIRT.
# The three fitters differ only in their A update and objective penalty;
# both are supplied as functions:
#   A_step(Params, Z, N) -> Params with updated A (and A_tilde if used)
#   penalty(A)           -> scalar added to the tracked objective
# use_tilde: linearize/track around A_tilde (the ADMM copy) instead of A.
# Returns Params with $value (last tracked objective) and, if Z_output,
# the final working response $Z.
GRM_sparse_HiMIRT <- function(Y, K, C, Params, A_step,
                              penalty = function(A) 0,
                              orthogonal = c(TRUE, FALSE),
                              use_tilde = FALSE, Z_output = FALSE,
                              alpha = 1, tol = 1e-4, max_iter,
                              verbose = FALSE, label = NULL,
                              print_init = FALSE) {
  N <- nrow(Y); J <- ncol(Y)
  A_x <- function(p) if (use_tilde) p$A_tilde else p$A

  X <- tcrossprod(Params$Theta, A_x(Params))
  tracker <- convergence_tracker(max_iter,
                                   logl_grm(Y, Params$tau, X, NULL) +
                                     penalty(Params$A),
                                   label = label, verbose = verbose,
                                   print_init = print_init)

  for (i in 1:max_iter) {
    ud <- update_delta(Params, Y, J, C, alpha, X, tracker$obj_vals[i],
                       function(p) logl_grm(Y, p$tau, X, NULL) + penalty(p$A))
    Params <- ud$Params

    Params <- proj_K(Params, C, K, item = TRUE)
    if (Params$proj) {
      X <- with(Params, tcrossprod(Theta, A))
    }

    Z <- compute_Z_grm(Y, X, Params$tau, C)

    Params <- update_Theta_grm(Params, Z, N, C, orthogonal,
                               use_tilde = use_tilde)

    Params <- A_step(Params, Z, N)

    X <- tcrossprod(Params$Theta, A_x(Params))
    if (check_converged(tracker, i,
                        logl_grm(Y, Params$tau, X, NULL) +
                          penalty(Params$A), tol)) break
  }

  if (Z_output) {
    Params$Z <- compute_Z_grm(Y, X, Params$tau, C)
  }

  # Last recorded objective (obj_vals is NA-padded past convergence).
  obj_hist <- tracker$obj_vals[!is.na(tracker$obj_vals)]
  Params$value <- obj_hist[length(obj_hist)]

  return(Params)
}

GRM_lasso_HiMIRT <- function(Y, K, C, Params = NULL,
                              lambda = 0, orthogonal = c(TRUE, FALSE),
                              alpha = 1, max_iter) {
  if (is.null(Params)) {
    Params <- GRM_HiMMIRT(Y, K, C, alpha = 1, rotation = TRUE, max_iter = 200, tol = 1e-3)
  }
  Params$A_tilde <- Params$A

  A_step <- function(Params, Z, N) {
    ADMM_res <- with(Params,
                     subADMM(Y, Z, A, A_tilde, Theta,
                             lambda = lambda, model = "GRM", N_iter = 20))
    Params$A <- ADMM_res$A
    Params$A_tilde <- ADMM_res$A_tilde
    Params
  }

  GRM_sparse_HiMIRT(Y, K, C, Params, A_step,
                    penalty = function(A) lambda * sum(abs(A)),
                    orthogonal = orthogonal, use_tilde = TRUE,
                    alpha = alpha, tol = 1e-4, max_iter = max_iter)
}

GRM_CHiMIRT <- function(Y,Params,Lambda,K,
                        orthogonal=c(TRUE,FALSE),Z_output=TRUE,
                        alpha=1,tol=1e-4,max_iter) {
  C <- ncol(Params$delta)+1

  Params$A <- ifelse(Lambda==0,0,Params$A)

  A_step <- function(Params, Z, N) {
    if (orthogonal) {
      # Theta has approximately orthonormal columns (post-Procrustes), so
      # crossprod(Theta) is nearly diagonal; the global closed-form * Lambda
      # projection is the correct constrained minimizer.
      A <- with(Params, t(Z)%*%Theta%*%safe_solve(crossprod(Theta)))
      Params$A <- A*Lambda
    } else {
      # Non-orth: masking the unconstrained solution is NOT the constrained
      # minimizer (off-diagonals in crossprod(Theta) couple active and
      # inactive entries). Solve per-item LS over each item's active
      # columns of Theta instead (same as the GPCM path). The GRM MM target
      # is Z itself; the curvature constant cancels out of the argmin.
      Params$A <- A_update_constrained(Z, Params$Theta, Lambda)
    }
    Params
  }

  GRM_sparse_HiMIRT(Y, K, C, Params, A_step,
                    orthogonal = orthogonal, Z_output = Z_output,
                    alpha = alpha, tol = tol, max_iter = max_iter)
}

GRM_autolasso_HiMIRT <- function(Y,K,C,Params=NULL,
                                 lambda_range,try_out=25,alpha=1,
                                 orthogonal=c(TRUE,FALSE)){
  lambda_values=seq(lambda_range[2],lambda_range[1],length.out=try_out)

  # Unpenalized refit on the active set Lambda; repeats until the KKT
  # condition holds (no zero loading with |gradient| > lambda). Pass
  # lambda = NULL to skip the KKT check (single refit).
  refit_active_set <- function(result, Lambda, lambda = NULL) {
    repeat {
      result <- GRM_CHiMIRT(Y,result,Lambda,K,
                            orthogonal=orthogonal,alpha=alpha,
                            max_iter=200)
      grad_A <- with(result,1/2*(A%*%crossprod(Theta)-crossprod(Z,Theta)))

      if (is.null(lambda)) break
      violations <- abs(grad_A)>lambda & result$A==0
      if (!any(violations)) break
      Lambda[violations] <- 1
    }
    list(result=result,grad_A=grad_A)
  }

  cat("fitting the most sparse model!\n")
  result <- GRM_lasso_HiMIRT(Y,K,C,Params,
                             lambda=lambda_range[2],orthogonal=orthogonal,alpha=alpha,
                             max_iter=200)

  cat("removing shrinkage!\n")
  fit <- refit_active_set(result, ifelse(result$A==0,0,1))

  result_set <- vector("list",length(lambda_values))
  names(result_set) <- lambda_values
  fit$result$Z <- NULL
  result_set[[1]] <- fit$result

  cat("trying out all lambda values given\n")

  for (i in 2:length(lambda_values)) {
    lambda=lambda_values[i]
    threshold=2*lambda-lambda_values[i-1]
    discarded_mat=abs(fit$grad_A)>=threshold

    if (sum(discarded_mat)!=0) {
      Lambda=ifelse(discarded_mat|fit$result$A!=0,1,0)
      fit <- refit_active_set(fit$result, Lambda, lambda)
    }

    fit$result$Z <- NULL
    result_set[[i]] <- fit$result
  }


  result_set$lambda_values=lambda_values
  cat("\u2714 Whole procedure finished!!!\n")

  return(result_set)
}

GRM_card_HiMIRT <- function(Y,K,C,Params=NULL,random_start=F,
                            kappaC,orthogonal=c(TRUE,FALSE),
                            tol=1e-4,alpha=1,max_iter){
  J <- ncol(Y)

  if (is.null(Params)){
    Params <- GRM_HiMMIRT(Y,K,C,alpha=1,rotation=T,max_iter=200,tol=1e-3)
  }

  if (random_start) {
    a_max <- max(abs(Params$A))
    Params$A <- a_max*matrix(stats::runif(J*K,-1,1),nrow=J)
  }

  kappaC0 <- J*K-kappaC

  A_step <- function(Params, Z, N) {
    if (orthogonal) {
      B <- crossprod(Z,Params$Theta)/N

      Params$A <- cardize(B,kappaC0)
    } else {
      # MM subiterations on the cardinality-constrained LS problem.
      Theta_squared <- crossprod(Params$Theta)
      alpha_card <- max(eigen(Theta_squared)$values)

      for (sub_i in 1:20){
        MM_loss_before <- with(Params,sum((Z-Theta%*%t(A))^2))
        B <- with(Params, t(t(A)-(Theta_squared%*%t(A)-t(Theta)%*%Z)/alpha_card))

        Params$A <- cardize(B,kappaC0)
        MM_loss_after <- with(Params, sum((Z-Theta%*%t(A))^2))

        relative_change_sub <- ifelse(sub_i==1,
                                      10,
                                      (MM_loss_before-MM_loss_after)/MM_loss_before)

        if (relative_change_sub<1e-5){
          break
        }
      }
    }
    Params
  }

  GRM_sparse_HiMIRT(Y, K, C, Params, A_step,
                    orthogonal = orthogonal,
                    alpha = alpha, tol = tol, max_iter = max_iter,
                    verbose = TRUE, label = "Theta+A", print_init = TRUE)
}

################################
#-----Partial Credit Model-----#
################################

# Dense exploratory fitter (the confirmatory / masked role now lives in
# GPCM_CHiMIRT, mirroring the GRM layout).
GPCM_HiMMIRT <- function(Y,K,C,
                         Params=NULL,
                         rotation=c(TRUE,FALSE),
                         orthogonal=TRUE,
                         max_iter,
                         tol=1e-4) {
  N <- nrow(Y); J <- ncol(Y); w <- C*(C^2-1)/24

  if (is.null(Params)){
    Params <- paramsINIT(Y,C,K,model="GPCM")
  }

  eta <- with(Params, tcrossprod(Theta, A))
  tracker <- convergence_tracker(max_iter, logl_gpcm(Y, eta, Params$tau, C),
                                   label = "Theta+A", verbose = TRUE, print_init = TRUE)

  for (i in 1:max_iter){
    # one tau step
    Params$tau <- update_tau_gpcm(Y, eta, Params$tau, N, J, C)

    # one working response, shared by the A and Theta steps
    Z <- compute_Z_gpcm(Y, eta, Params$tau, C, w)$Z

    # one A step
    Params$A <- crossprod(Z, Params$Theta) %*% safe_solve(crossprod(Params$Theta))

    # one Theta step, against the same Z
    if (orthogonal) {
      svd_result <- svd(1/sqrt(N) * Z %*% Params$A)
      Params$Theta <- sqrt(N) * svd_result$u %*% t(svd_result$v)

      Params <- with(Params,
                      standardize(Theta, A, tau = tau,
                                  C = C, model = "GPCM"))
    } else {
      Params$Theta <- Z %*% Params$A %*% safe_solve(crossprod(Params$A))

      if (rotation==F) {
        Params <- with(Params,
                        standardize(Theta, A, tau = tau,
                                    C = C, model = "GPCM"))
      }
    }

    eta <- with(Params,tcrossprod(Theta, A))
    if (check_converged(tracker, i, logl_gpcm(Y, eta, Params$tau, C), tol)) break
  }

  # Last recorded objective (obj_vals is NA-padded past convergence). Named
  # $obj rather than $value to keep it distinct from the sparse fitters'
  # $value, which includes a penalty term; this one is the unpenalized
  # negative log-likelihood. The SVD/rotation step below leaves
  # Theta %*% t(A) and the implied probabilities unchanged, so $obj stays
  # valid for the returned parameters.
  obj_hist   <- tracker$obj_vals[!is.na(tracker$obj_vals)]
  Params$obj <- obj_hist[length(obj_hist)]

  if (rotation) {
    svd_result <- svd(Params$Theta)
    Params$Theta <- sqrt(N)*svd_result$u
    Params$A <- Params$A%*%svd_result$v%*%diag(svd_result$d)/sqrt(N)

    rotate_list <- GPArotation::GPFoblq(Params$A,method="geomin")
    Params$A <- rotate_list$loadings # use geomin rotation
    rotation_M <- rotate_list$Th
    Params$Theta <- Params$Theta %*% rotation_M
  }

  Params$beta <- -Linv(Params$tau,C)

  return(Params)
}

# Core MM loop shared by GPCM_lasso_HiMIRT / GPCM_CHiMIRT / GPCM_card_HiMIRT.
# Mirrors GRM_sparse_HiMIRT: the fitters differ only in their A update and
# objective penalty, supplied as functions:
#   A_step(Params, Z, N) -> Params with updated A (and A_tilde if used);
#                           Z = N x J working response from compute_Z_gpcm()
#   penalty(A)           -> scalar added to the tracked objective
#
# One tau step, one Z, one Theta step and one A step per outer iteration; the
# Theta and A steps share the same Z.
GPCM_sparse_HiMIRT <- function(Y, K, C, Params, A_step,
                               penalty = function(A) 0,
                               orthogonal = c(TRUE, FALSE),
                               use_tilde = FALSE, Z_output = FALSE,
                               tol = 1e-4, max_iter,
                               verbose = FALSE, label = NULL,
                               print_init = FALSE) {
  N <- nrow(Y); J <- ncol(Y); w <- C*(C^2-1)/24
  A_x <- function(p) if (use_tilde) p$A_tilde else p$A

  eta <- tcrossprod(Params$Theta, A_x(Params))
  tracker <- convergence_tracker(max_iter,
                                   logl_gpcm(Y, eta, Params$tau, C) +
                                     penalty(Params$A),
                                   label = label, verbose = verbose,
                                   print_init = print_init)

  for (i in 1:max_iter) {
    Params$tau <- update_tau_gpcm(Y, eta, Params$tau, N, J, C)

    Z <- compute_Z_gpcm(Y, eta, Params$tau, C, w)$Z

    Params <- update_Theta_gpcm(Params, Z, N, C, orthogonal,
                                use_tilde = use_tilde)

    Params <- A_step(Params, Z, N)

    eta <- tcrossprod(Params$Theta, A_x(Params))
    if (check_converged(tracker, i,
                        logl_gpcm(Y, eta, Params$tau, C) + penalty(Params$A),
                        tol)) break
  }

  if (Z_output) {
    # Store Hw, not Z. At the tangency point (eta = Theta A'), the surrogate
    # gradient w*(A'Theta'Theta - Z'Theta) collapses exactly to -crossprod(Hw,
    # Theta) -- the w terms cancel -- so Hw alone gives the log-likelihood
    # gradient wrt A. Same N x J footprint, no constant, and immune to any
    # A / A_tilde mismatch on the lasso path.
    Params$Hw <- compute_Z_gpcm(Y, eta, Params$tau, C, w)$Hw
  }

  # Last recorded objective (obj_vals is NA-padded past convergence).
  obj_hist <- tracker$obj_vals[!is.na(tracker$obj_vals)]
  Params$value <- obj_hist[length(obj_hist)]
  Params$beta <- -Linv(Params$tau, C)

  return(Params)
}

# Confirmatory fitter: fixed sparsity pattern Lambda (0 = fixed at zero).
# Mirrors GRM_CHiMIRT; extracted from the old GPCM_HiMMIRT Lambda path.
# One masked A step per outer iteration, against the Z built by the caller.
GPCM_CHiMIRT <- function(Y,Params,Lambda,K,
                         orthogonal=c(TRUE,FALSE),Z_output=TRUE,
                         tol=1e-4,max_iter) {
  J <- ncol(Y); C <- length(Params$tau)/J

  Params$A <- ifelse(Lambda==0,0,Params$A)

  A_step <- function(Params, Z, N) {
    if (orthogonal) {
      # Post-Procrustes Theta: closed-form * Lambda is the constrained
      # minimizer (crossprod(Theta) nearly diagonal).
      Params$A <- (crossprod(Z, Params$Theta) %*%
                     safe_solve(crossprod(Params$Theta))) * Lambda
    } else {
      # Non-orth: per-item LS over each item's active columns of Theta.
      Params$A <- A_update_constrained(Z, Params$Theta, Lambda)
    }
    Params
  }

  GPCM_sparse_HiMIRT(Y, K, C, Params, A_step,
                     orthogonal = orthogonal, Z_output = Z_output, verbose = F,
                     tol = tol, max_iter = max_iter)
}

GPCM_lasso_HiMIRT <- function(Y,K,C,Params=NULL,
                              lambda=NULL,orthogonal=c(TRUE,FALSE),
                              max_iter) {
  if (is.null(Params)){
    Params <- GPCM_HiMMIRT(Y,K,C,rotation=T,orthogonal=F,max_iter=200,tol=1e-3)
  }

  Params$A_tilde <- Params$A

  A_step <- function(Params, Z, N) {
    ADMM_res <- with(Params,
                     subADMM(Y, Z, A, A_tilde, Theta, C = C, lambda,
                        model = "GPCM", N_iter = 20))
    Params$A       <- ADMM_res$A
    Params$A_tilde <- ADMM_res$A_tilde
    Params
  }

  GPCM_sparse_HiMIRT(Y, K, C, Params, A_step,
                     penalty = function(A) lambda * sum(abs(A)),
                     orthogonal = orthogonal, use_tilde = TRUE,
                     tol = 1e-4, max_iter = max_iter,
                     verbose = F, label = "A", print_init = F)
}

GPCM_autolasso_HiMIRT <- function(Y,K,C,Params,
                                  lambda_range,try_out=25,
                                  orthogonal=c(TRUE,FALSE)) {
  w <- C*(C^2-1)/24
  lambda_values <- seq(lambda_range[2],lambda_range[1],length.out=try_out)

  # Unpenalized refit on the active set Lambda; repeats until the KKT
  # condition holds (no zero loading with |gradient| > lambda). Pass
  # lambda = NULL to skip the KKT check (single refit).
  refit_active_set <- function(result, Lambda, lambda = NULL) {
    repeat {
      result <- GPCM_CHiMIRT(Y,result,Lambda,K,
                             orthogonal=orthogonal,max_iter=200)
      # gradient of the LOSS wrt A (= -grad of the log-likelihood), which is
      # what the lasso KKT check below needs
      grad_A <- with(result, -crossprod(Hw, Theta))

      if (is.null(lambda)) break
      violations <- abs(grad_A)>lambda & result$A==0
      if (!any(violations)) break
      Lambda[violations] <- 1
    }
    list(result=result,grad_A=grad_A)
  }

  cat("fitting the most sparse model!\n")
  result <- GPCM_lasso_HiMIRT(Y,K,C,Params=Params,
                              lambda=lambda_range[2],
                              orthogonal=orthogonal,max_iter=200)

  cat("removing shrinkage!\n")
  fit <- refit_active_set(result, sign(result$A!=0))

  result_set <- vector("list",length(lambda_values))
  names(result_set) <- lambda_values
  fit$result$Hw <- NULL   # drop the N x J array before storing the fit
  result_set[[1]] <- fit$result

  cat("trying out all lambda values given\n")

  for (i in 2:length(lambda_values)) {
    lambda=lambda_values[i]
    threshold=2*lambda-lambda_values[i-1]
    discarded_mat=abs(fit$grad_A)>=threshold

    if (sum(discarded_mat)!=0) {
      Lambda=ifelse(discarded_mat|fit$result$A!=0,1,0)
      fit <- refit_active_set(fit$result, Lambda, lambda)
    }

    fit$result$Hw <- NULL   # drop the N x J array before storing the fit
    result_set[[i]] <- fit$result
  }

  result_set$lambda_values=lambda_values
  cat("\u2714 Whole procedure finished!!!\n")

  return(result_set)
}

GPCM_card_HiMIRT <- function(Y,K,C,Params=NULL,random_start=F,
                             kappaC=0,orthogonal=c(TRUE,FALSE),
                             tol=1e-4,
                             max_iter) {
  N <- nrow(Y); J <- ncol(Y);kappaC0=J*K-kappaC
  # NB: w is deliberately absent. The cardinality A-step is a constrained LS
  # argmin against Z, and the majorization constant scales the objective, so
  # it cancels from both the closed form and the proximal step.
  # GPCM_sparse_HiMIRT computes its own w for the Z it builds.

  if (is.null(Params)){
    Params <- GPCM_HiMMIRT(Y,K,C,rotation=T,orthogonal=T,max_iter=200,tol=1e-3)
  }

  if (random_start) {
    a_max <- max(abs(Params$A))
    Params$A <- a_max*matrix(stats::runif(J*K,-1,1),nrow=J)
  }

  A_step <- function(Params, Z, N) {
    if (orthogonal) {
      # crossprod(Theta) = N * I on the orthogonal path, so the unconstrained
      # LS solution is Z'Theta / N. w cancels out of the argmin.
      B <- crossprod(Z, Params$Theta) / N

      Params$A <- cardize(B, kappaC0)
    } else {
      # MM subiterations on the cardinality-constrained LS problem. Unlike the
      # orthogonal branch -- where crossprod(Theta) = N * I makes the problem
      # separable and a single cardize() of the closed form is exact -- here
      # the off-diagonals couple the active and inactive entries, so one
      # proximal step is NOT the constrained minimizer of the surrogate. The
      # loop has to stay.
      #
      # Step is A - grad/L with grad = w(A Theta'Theta - Z'Theta) and
      # L = w * lambda_max(Theta'Theta), so w cancels top and bottom.
      Theta_squared <- crossprod(Params$Theta)
      alpha <- max(eigen(Theta_squared)$values)

      # The MM surrogate ||Z - Theta A'||^2. The w/2 in front scales the whole
      # thing and cancels from both the argmin and the relative test below, so
      # it is omitted. Non-negative, which makes the ratio well defined.
      mm_var <- function(A_curr) {
        sum((Z - tcrossprod(Params$Theta, A_curr))^2)
      }

      for (sub_i in 1:20) {
        MM_loss_before <- mm_var(Params$A)

        B <- Params$A -
          (Params$A %*% Theta_squared - crossprod(Z, Params$Theta)) / alpha

        Params$A <- cardize(B, kappaC0)
        MM_loss_after <- mm_var(Params$A)

        # Relative decrease, matching the GRM path above: scale-free, so
        # unaffected by w, N, J or C.
        rel_change <- (MM_loss_before - MM_loss_after) /
          max(MM_loss_before, .Machine$double.eps)
        if (sub_i > 1 && rel_change < 1e-5) {
          break
        }
      }
    }
    Params
  }

  GPCM_sparse_HiMIRT(Y, K, C, Params, A_step,
                     orthogonal = orthogonal,
                     tol = tol, max_iter = max_iter,
                     verbose = TRUE, label = "Theta+A", print_init = TRUE)
}
