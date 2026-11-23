# poly_mirt_wrappers.R
#
# R wrappers that expose the same API as the original Basics.R functions
# (logl, compute_H) but dispatch to the Rcpp implementations in
# poly_mirt_rcpp.cpp.
#
# Usage:
#   library(Rcpp)
#   sourceCpp("R/poly_mirt_rcpp.cpp")
#   source("R/poly_mirt_wrappers.R")
#
#   # Drop-in replacement for logl() and compute_H()
#   val <- logl(Y, model = "GPCM", par = list(C = C, eta = eta, tau = tau))
#   H   <- compute_H(Y, model = "GRM", par = list(tau = tau, P_all = P_all))

## Validate par (same logic as Basics.R)

validate_par <- function(par, model) {
  model <- match.arg(model, c("GRM", "GPCM"))
  req <- switch(model,
                GRM  = c("tau", "X", "P_all"),
                GPCM = c("C", "eta", "tau"))
  missing <- setdiff(req, names(par))
  if (length(missing)) {
    stop("Missing in par for ", model, ": ", paste(missing, collapse = ", "))
  }
  par
}

## logl wrapper

logl <- function(Y, model = c("GRM", "GPCM"), par) {
  model <- match.arg(model, c("GRM", "GPCM"))
  par   <- validate_par(par, model)

  if (model == "GRM") {
    tau   <- par$tau
    X     <- par$X
    P_all <- par$P_all

    if (!is.null(P_all)) {
      return(logl_grm(Y, tau, NULL, P_all))
    } else if (!is.null(X)) {
      return(logl_grm(Y, tau, X, NULL))
    } else {
      stop("GRM logl requires either X or P_all in par")
    }
  } else {
    # GPCM
    C   <- par$C
    eta <- par$eta
    tau <- par$tau
    return(logl_gpcm(Y, eta, tau, C))
  }
}

## compute_H wrapper

compute_H <- function(Y, model = c("GRM", "GPCM"), par) {
  model <- match.arg(model, c("GRM", "GPCM"))
  par   <- validate_par(par, model)

  if (model == "GRM") {
    P_all <- par$P_all
    C     <- ncol(par$tau) - 1
    return(compute_H_grm(Y, P_all, C))
  } else {
    # GPCM
    C   <- par$C
    eta <- par$eta
    tau <- par$tau
    return(compute_H_gpcm(Y, eta, tau, C))
  }
}
