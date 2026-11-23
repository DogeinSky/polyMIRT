#generate_q <- function(K,zero_prob){
#  first <- rbinom(K-1, 1, 1-zero_prob)
#  last <- ifelse(sum(first)==0, 1, rbinom(1,1,1-zero_prob))
#  c(first, last)
#}

generate_q <- function(K, zero_prob) {
  q <- rbinom(K, 1, 1 - zero_prob)
  if (sum(q) == 0) q[sample.int(K, 1)] <- 1
  q
}

generate_d <- function(J,C){
  d_true <- unlist(replicate(J, c(0, runif(C-1,-2,2)), simplify = FALSE))
  return(d_true)
}

generate_beta <- function(J,C){
  beta_true <- unlist(replicate(J,
                                c(0, MASS::mvrnorm(n=1,mu=1:(C-1) - C/2,
                                             Sigma=diag(0.4,C-1))),
                                simplify=FALSE))

  return(beta_true)
}

# generate_tau ---------------------------------------------------------------
#
# Item-specific GRM thresholds, drawn at random instead of the fixed
#   matrix(c(-1e10,-2,-0.5,0.5,2,1e10), nrow = J, ncol = C+1, byrow = TRUE)
#
# Each of the C-1 interior thresholds is drawn independently per item from
# its own uniform interval. The default intervals (for C = 5) are
#
#   tau_1 ~ U(-2.5, -1.5)    tau_3 ~ U( 0.25, 1.25)
#   tau_2 ~ U(-1.25, -0.25)  tau_4 ~ U( 1.5,  2.5)
#
# These are disjoint and increasing, so every row comes out strictly
# ordered automatically -- no rejection or re-sorting needed. `sort_rows`
# is a safety net for custom `bounds` that do overlap.
#
# Args:
#   J        number of items (rows)
#   C        number of response categories; C-1 interior thresholds
#   bounds   (C-1) x 2 matrix of (lower, upper) per threshold. Defaults to
#            the intervals above, which require C = 5.
#   sort_rows  sort each row ascending after drawing (no-op for disjoint,
#            increasing bounds; keeps validity if bounds overlap)
#
# Returns: J x (C+1) matrix, -1e10 and 1e10 in the outer columns, matching
# the layout generate_data(model = "GRM") expects.
generate_tau <- function(J, C = 5, bounds = NULL, sort_rows = TRUE) {
  if (is.null(bounds)) {
    if (C != 5) {
      stop("default bounds are defined for C = 5 only; ",
           "supply a (C-1) x 2 `bounds` matrix for C = ", C)
    }
    bounds <- rbind(c(-2.50, -1.50),
                    c(-1.25, -0.25),
                    c( 0.25,  1.25),
                    c( 1.50,  2.50))
  }

  bounds <- as.matrix(bounds)
  if (ncol(bounds) != 2 || nrow(bounds) != C - 1) {
    stop("`bounds` must be a (C-1) x 2 matrix; got ",
         nrow(bounds), " x ", ncol(bounds), " for C = ", C)
  }
  if (any(bounds[, 1] > bounds[, 2])) stop("each row of `bounds` must be (lower, upper)")

  # one runif(J) per threshold, in threshold order
  tau_fin <- matrix(
    vapply(seq_len(C - 1),
           function(k) stats::runif(J, bounds[k, 1], bounds[k, 2]),
           numeric(J)),
    nrow = J, ncol = C - 1)

  if (sort_rows && C - 1 >= 2) {
    tau_fin <- t(apply(tau_fin, 1, sort))
  }

  cbind(-1e10, tau_fin, 1e10)
}

# min_pure: reject draws in which any factor has fewer than this many PURE
#   items (items loading on it and nothing else). Default 0 reproduces the
#   previous unconditional behaviour exactly, so existing scripts are unchanged.
#
#   Why gate. With zero_prob = 0.7 an item has Binomial(K, 0.3) nonzeros, so
#   P(exactly one) = K(0.3)(0.7^(K-1)); at J = 50, K = 5 that leaves only
#   ~3.6 pure items per factor in expectation, and over 20000 draws
#       P(min pure <= 2) = 39.7%,  P(min pure <= 1) = 12.1%,  P(min = 0) = 2.0%
#   -- i.e. 2% of draws are not even identified in the IC(5) sense, and a
#   weakly identified draw governs every replicate of a cell. min_pure = 3 is
#   the classical three-indicator rule for factor identification.
#
#   The gate binds ONLY at J = 50, K = 5 (acceptance 60% there, >=97% in every
#   other cell of the design), so it leaves the rest of the design untouched.
#   It does condition the loading distribution slightly: at J = 50, K = 5 the
#   mean number of nonzeros falls from 83.3 to 81.7, about 1.9% sparser.
#
#   The returned matrix carries attr "pure" (pure items per factor) and
#   attr "tries", so a simulation can log the realised structure.
generate_A <- function(J, K, zero_prob, min_pure = 0, max_try = 1000){
    draw <- function() t(vapply(1:J,
                       function(j) generate_q(K,zero_prob)*stats::runif(K, 0.5, 2.5),
                       numeric(K)))
    pure_count <- function(A)
        colSums(A[rowSums(A != 0) == 1, , drop = FALSE] != 0)

    if (min_pure <= 0) {
        A_true <- draw()
        return(structure(A_true, pure = pure_count(A_true), tries = 1L))
    }
    for (tr in seq_len(max_try)) {
        A_true <- draw()
        pure   <- pure_count(A_true)
        if (min(pure) >= min_pure)
            return(structure(A_true, pure = pure, tries = tr))
    }
    stop(sprintf(paste("generate_A: no draw with >= %d pure items per factor",
                       "in %d tries (J = %d, K = %d, zero_prob = %.2f).",
                       "Lower min_pure, raise J, or lower zero_prob."),
                 min_pure, max_try, J, K, zero_prob))
}

generate_data <- function(N,J,rho,A_true,model=c("GRM","GPCM"),theta_distr=NULL,
                          C=NULL,d_true=NULL,tau_true=NULL,beta_true=NULL){
  model <- match.arg(model, c("GRM","GPCM"))

  K <- ncol(A_true)

  # rho: scalar            -> constant correlation for all pairs (0 = orthogonal);
  #      length 2           -> interval: each of the choose(K,2) pairwise
  #                            correlations drawn from runif(rho[1], rho[2]);
  #      length choose(K,2) -> fixed pairwise correlations, used as-is
  #                            (fill order: lower triangle, column-major).
  # Note choose(K,2) is never 2, so the interval form is unambiguous.
  n_pairs <- K * (K - 1) / 2
  if (length(rho) == 1) {
    rho_vec <- rep(rho, n_pairs)
  } else if (length(rho) == n_pairs) {
    rho_vec <- rho
  } else if (length(rho) == 2) {
    rho_vec <- runif(n_pairs, rho[1], rho[2])
  } else {
    stop("rho must be a scalar, a length-2 interval c(lo, hi), ",
         "or a vector of choose(K,2) fixed correlations")
  }

  cov_mat <- diag(K)
  cov_mat[lower.tri(cov_mat)] <- rho_vec
  cov_mat[upper.tri(cov_mat)] <- t(cov_mat)[upper.tri(cov_mat)]

  if (min(eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values) < 1e-8) {
    stop("sampled correlation matrix is not positive definite; ",
         "narrow the rho interval")
  }

  if (is.null(theta_distr)||theta_distr=="normal"){
    Theta_true <- MASS::mvrnorm(n = N, mu = rep(0,K), Sigma = cov_mat)

    if (all(rho == 0)) {
      Theta_true=scale(svd(Theta_true)$u)
    } else {
      Theta_true=scale(Theta_true)
    }
  } else if (theta_distr=="gamma"){
    Theta_true <- scale(matrix(rgamma(N * K, shape = 2), N, K))
  } else if (theta_distr=="bimodal"){
    mu <- ifelse(rbinom(N * K, 1, 0.5) == 1, 2, -2)
    Theta_true <- scale(matrix(rnorm(N * K, mean = mu, sd = 1), N, K))
  }

  if (model=="GRM") {
    epsilon <- matrix(rlogis(N*J),nrow=N,ncol=J)
    Response_true <- Theta_true%*%t(A_true)+epsilon

    Y <- vapply(1:ncol(Response_true),
                function(i) findInterval(Response_true[,i],tau_true[i,]),
                numeric(nrow(Response_true)))

    return(list(Y=Y,Theta_true=Theta_true,A_true=A_true,tau_true=tau_true,
                cov_true=cov_mat))
  } else {
    B <- matrix(beta_true, nrow = C, ncol = J)
    tau_true <- as.vector(apply(B,2,cumsum))

    eta  <- tcrossprod(Theta_true, A_true)
    cvec <- t(0:(C - 1))

    # Row-blocked. The full N x J*C kronecker / exp arrays were the memory
    # peak of this function (~800 MB at N = 12500, J = 500, C = 5, i.e. more
    # than the fitter itself). Blocking keeps the vectorized exp() but holds
    # only block_size x J*C at a time.
    #
    # Rows are still visited in order 1..N and items in order 1..J, so the
    # rmultinom() draw sequence is untouched: for a given seed Y is bitwise
    # identical to the unblocked version, and the RNG is left in the same
    # state afterwards.
    block_size <- max(1L, min(N, as.integer(ceiling(2e6 / (J * C)))))
    Y <- matrix(0, nrow = N, ncol = J)

    for (start in seq(1L, N, by = block_size)) {
      idx <- start:min(start + block_size - 1L, N)
      exp_X <- exp(t(t(kronecker(eta[idx, , drop = FALSE], cvec)) - tau_true))

      for (b in seq_along(idx)) {
        Y[idx[b], ] <- vapply(1:J, function(j){
          prob <- exp_X[b,((j-1)*C+1):(j*C)]
          which.max(rmultinom(1,1,prob))
        },numeric(1))
      }
    }

    return(list(Y=Y,Theta_true=Theta_true,
                A_true=A_true,beta_true=beta_true,
                tau_true=tau_true,cov_true=cov_mat))
  }
}

ref.perp <- function(A_hat,A_true,Theta_hat){
  K=ncol(A_hat)
  matrix_data=t(vapply(1:K, function(i){
    colSums(ifelse(A_hat[,i]!=0 & A_true!=0, 1,0))},
    numeric(K)))
  assignment <- clue::solve_LSAP(t(matrix_data), maximum = TRUE)
  row_indices <- as.integer(assignment)
  A_best <- A_hat[,row_indices]
  Theta_best <- Theta_hat[,row_indices]

  cc <- diag(suppressWarnings(cor(A_true, A_best)))
  cc[!is.finite(cc)] <- 1
  corsign <- sign(cc); corsign[corsign == 0] <- 1


  A_best <- A_best %*% diag(corsign)
  Theta_best <- Theta_best %*% diag(corsign)
  return(list(A_best=A_best,Theta_best=Theta_best))
}

showResults <- function(result_list,Data){
  df <- data.frame(
    lambda=result_list$lambda_values,
    sens=NA,
    spec=NA,
    load_diff=NA
  )

  df[2:4] <- t(vapply(result_list[1:length(result_list$lambda_values)], function(result) {
    A_hat <- result$A
    Theta_hat <- result$Theta
    tau_hat <- result$tau

    A_true <- Data$A_true
    tau_true <- Data$tau_true
    Theta_true <- Data$Theta_true

    A_best=ref.perp(A_hat,A_true,Theta_hat)$A_best

    sens <- sum(ifelse(A_best!=0&A_true!=0, 1, 0))/sum(ifelse(A_true!=0,1,0))
    spec <- sum(ifelse(A_best==0&A_true==0, 1, 0))/sum(ifelse(A_true==0,1,0))

    load_diff <- sum((A_best-A_true)^2)/(J*K)

    #BIC = 2*result$logl_value+(sum(result$A!=0)+J+N*K)*log(N)
    return(c(sens,spec,load_diff))
  }, numeric(3)))

  pt_col <- rep("black", length(result_list$lambda_values))
  pt_col[which.min(df$BIC)] <- "red"

  #plot(df$lambda,df$prob_diff,col=pt_col,pch=16)
  #lines(df$lambda,df$prob_diff)

  plot(df$lambda,df$sens,col=pt_col,pch=16)
  lines(df$lambda,df$sens)

  plot(df$lambda,df$spec,col=pt_col,pch=16)
  lines(df$lambda,df$spec)

  plot(df$lambda,df$load_diff,col=pt_col,pch=16)
  lines(df$lambda,df$load_diff)

  #plot(df$lambda,df$BIC,col=pt_col,pch=16)
  #lines(df$lambda,df$BIC)
}

compute_P_gpcm <- function(eta, beta, C) {
  N <- nrow(eta); J <- ncol(eta)
  beta_mat <- matrix(beta, nrow = C, ncol = J)
  tau      <- as.vector(apply(beta_mat, 2, cumsum))
  tau_mat <- matrix(tau, nrow = C, ncol = J)     # tau_mat[c, j]

  # Pass 1: streaming log-sum-exp -> final running max M and normaliser S.
  # c = 1: X[, 1, ] = -tau_1 (broadcast across rows).
  M <- matrix(-tau_mat[1L, ], N, J, byrow = TRUE)
  S <- matrix(1, N, J)                           # exp(X_1 - M) = 1
  if (C >= 2L) for (c_idx in 2:C) {
    Xc    <- (c_idx - 1L) * eta -
      matrix(tau_mat[c_idx, ], N, J, byrow = TRUE)
    M_new <- pmax(M, Xc)
    S     <- S * exp(M - M_new) + exp(Xc - M_new)
    M     <- M_new
  }

  # Pass 2: write each softmax slice into the output P (N x J*C, viewed as
  # N x C x J via the O(1) dim<- so column-major lay-out matches the rest of
  # the codebase: P[i, (j-1)C + c] is at P[i, c, j]).
  P <- matrix(0, N, J * C)
  dim(P) <- c(N, C, J)
  Xc <- matrix(-tau_mat[1L, ], N, J, byrow = TRUE)
  P[, 1L, ] <- exp(Xc - M) / S
  if (C >= 2L) for (c_idx in 2:C) {
    Xc <- (c_idx - 1L) * eta -
      matrix(tau_mat[c_idx, ], N, J, byrow = TRUE)
    P[, c_idx, ] <- exp(Xc - M) / S
  }
  dim(P) <- c(N, J * C)
  P
}

# calSelect -----------------------------------------------------------------
#
# The two quantities needed to PICK a lambda off the regularization path, and
# nothing else. calResults() returns nine metrics, of which the path only ever
# uses card_diff and BIC -- the other seven are computed once per lambda and
# discarded. That waste is not free: calResults() builds two N x J x C
# probability arrays per call (11M doubles at J=300/sc=25/C=5, ~90 MB each),
# so 15 lambdas x 10 parallel workers is enough transient memory to get
# workers killed, which silently costs whole replicates.
#
# Neither quantity needs ref.perp():
#   * card_diff counts nonzeros of A_best, but ref.perp() only permutes the
#     columns of A_hat and multiplies them by +/-1, and neither operation can
#     change a nonzero count. So sum(A_best != 0) == sum(A_hat != 0) exactly.
#   * BIC already uses A_hat and result$value directly, not A_best.
#
# Returns the same two names calResults() does, so
#   which.min(df$card_diff) / which.min(df$BIC)
# are unchanged. Model-agnostic: nothing here depends on GRM vs GPCM.
calSelect <- function(result, Data, N, K) {
  nnz <- sum(result$A != 0)
  c(card_diff = abs(nnz - sum(Data$A_true != 0)),
    BIC       = 2*result$value + (nnz + N*K)*log(N))
}

calResults <- function(result,Data,N,J,K,C,model=c("GRM","GPCM")){
  A_hat <- result$A
  A_true <- Data$A_true
  Theta_true <- Data$Theta_true
  rho_true <- Data$cov_true[upper.tri(Data$cov_true)]

  best_res <- ref.perp(A_hat,A_true,result$Theta)
  A_best <- best_res$A_best
  Theta_best <- best_res$Theta_best

  cov_hat <- crossprod(Theta_best)/(N-1)
  rho_hat <- cov_hat[upper.tri(cov_hat)]

  sens <- sum(A_best!=0 & A_true!=0)/sum(A_true!=0)
  spec <- sum(A_best==0 & A_true==0)/sum(A_true==0)
  load_diff <- sum((A_best-A_true)^2)/(J*K)
  Theta_diff <- sum((Theta_best-Theta_true)^2)/(N*K)
  card_diff <- abs(sum(A_best!=0)-sum(A_true!=0))
  rho_diff <- sum((rho_true-rho_hat)^2)/(length(rho_true))

  if (model=="GRM") {
    tau_diff <- sum((result$tau[,2:C]-Data$tau_true[,2:C])^2)/(J*(C-1))

    X_true <- tcrossprod(Theta_true,A_true)
    PFf_true <- compute_PFf(X_true,Data$tau_true)
    P_true <- PFf_true$P_all
    X_hat <- tcrossprod(Theta_best,A_best)
    PFf_hat <- compute_PFf(X_hat,result$tau)
    P_hat <- PFf_hat$P_all
    prob_diff <- sum((P_hat-P_true)^2)/(length(P_hat))
    BIC <- 2*result$value+(sum(A_hat!=0)+N*K)*log(N)

    return(c(prob_diff=prob_diff,tau_diff=tau_diff,Theta_diff=Theta_diff,
             sens=sens,spec=spec,load_diff=load_diff,rho_diff=rho_diff,
             card_diff=card_diff,BIC=BIC))
  } else {
    beta_diff <- sum((result$beta - Data$beta_true)^2) / (J*(C-1))

    eta_true <- tcrossprod(Theta_true, A_true)
    eta_hat  <- tcrossprod(Theta_best, A_best)
    P_true <- compute_P_gpcm(eta_true, Data$beta_true, C)
    P_hat  <- compute_P_gpcm(eta_hat,  result$beta,  C)
    prob_diff <- sum((P_hat - P_true)^2) / length(P_hat)
    BIC <- 2*result$value+(sum(A_hat!=0)+N*K)*log(N)

    return(c(prob_diff=prob_diff,beta_diff=beta_diff,Theta_diff=Theta_diff,
             sens=sens,spec=spec,load_diff=load_diff,rho_diff=rho_diff,
             card_diff=card_diff,BIC=BIC))
  }
}

calSelect <- function(result, Data, N, K) {
  nnz <- sum(result$A != 0)
  c(card_diff = abs(nnz - sum(Data$A_true != 0)),
    BIC       = 2*result$value + (nnz + N*K)*log(N))
}
