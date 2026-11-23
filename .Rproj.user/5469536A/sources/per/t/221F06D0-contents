// poly_mirt_rcpp.cpp
//
// Combined Rcpp port of compute_PFf(), compute_H() (GRM & GPCM branches),
// and logl() (GRM & GPCM branches) from Basics.R.
//
// Compile and load in R with:
//   library(Rcpp)
//   sourceCpp("R/poly_mirt_rcpp.cpp")

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

using namespace Rcpp;

// ============================================================================
// compute_PFf  (GRM branch)
// ============================================================================
// [[Rcpp::export]]
List compute_PFf(NumericMatrix X, NumericMatrix tau, bool need_f = true) {
    const int N = X.nrow();
    const int J = X.ncol();
    const int C = tau.ncol() - 1;
    const double eps = std::numeric_limits<double>::epsilon();

    NumericMatrix P_all(N, J * C);
    NumericMatrix f_all_mat(need_f ? N : 0, need_f ? J * (C - 1) : 0);

    const double* p_X = REAL(X);
    const double* p_tau = REAL(tau);
    double* p_P = REAL(P_all);
    double* p_f = need_f ? REAL(f_all_mat) : nullptr;

    std::vector<double> F_prev(N);

    for (int j = 0; j < J; ++j) {
        std::fill(F_prev.begin(), F_prev.end(), 0.0);
        const double* X_col = p_X + static_cast<std::size_t>(j) * N;

        for (int c = 0; c < C - 1; ++c) {
            const double tau_jc = p_tau[static_cast<std::size_t>(c + 1) * J + j];
            double* P_col = p_P + static_cast<std::size_t>(j * C + c) * N;

            if (need_f) {
                double* f_col = p_f + static_cast<std::size_t>(j * (C - 1) + c) * N;
                for (int i = 0; i < N; ++i) {
                    const double z = tau_jc - X_col[i];
                    const double F_curr = 1.0 / (1.0 + std::exp(-z));
                    double P = F_curr - F_prev[i];
                    if (P < eps) P = eps;
                    P_col[i] = P;
                    f_col[i] = F_curr * (1.0 - F_curr);
                    F_prev[i] = F_curr;
                }
            } else {
                for (int i = 0; i < N; ++i) {
                    const double z = tau_jc - X_col[i];
                    const double F_curr = 1.0 / (1.0 + std::exp(-z));
                    double P = F_curr - F_prev[i];
                    if (P < eps) P = eps;
                    P_col[i] = P;
                    F_prev[i] = F_curr;
                }
            }
        }

        double* P_col_last = p_P + static_cast<std::size_t>(j * C + (C - 1)) * N;
        for (int i = 0; i < N; ++i) {
            double P = 1.0 - F_prev[i];
            if (P < eps) P = eps;
            P_col_last[i] = P;
        }
    }

    if (need_f) {
        return List::create(_["P_all"] = P_all, _["f_all"] = f_all_mat);
    } else {
        return List::create(_["P_all"] = P_all, _["f_all"] = R_NilValue);
    }
}

// ============================================================================
// tau_stats_grm  (fused GRM tau-step statistics)
// ============================================================================
// One streaming pass that returns exactly what the tau update needs:
//
//   G        = grad_tau(Y, J, C, PFf, clip_value)        (J x (C-1))
//   diag_vec = fisher_tau(J, C, PFf)$diag_vec            (J*(C-1))
//   off_vec  = fisher_tau(J, C, PFf)$off_vec             (J*(C-2))
//
// Replaces compute_PFf(need_f = TRUE) + grad_tau() + fisher_tau(). The
// N x J*C and N x J*(C-1) arrays (and the ~6 full-width temporaries the R
// versions build to reduce them) are never materialized: per row only C
// probabilities and C-1 densities are held in scratch.
//
// Accumulators are long double to match R's colSums(), which reduces in
// LDOUBLE; expression shapes below mirror the R code term by term
// (f2/P_l + f2/P_r, not f2*(1/P_l + 1/P_r)) so results match bitwise.
// [[Rcpp::export]]
List tau_stats_grm(IntegerMatrix Y, NumericMatrix X, NumericMatrix tau,
                   double clip_value) {
    const int N = X.nrow();
    const int J = X.ncol();
    const int C = tau.ncol() - 1;
    const double eps = std::numeric_limits<double>::epsilon();

    if (Y.nrow() != N || Y.ncol() != J) stop("Y dimensions must match X");
    if (C < 2) stop("C must be >= 2");

    NumericMatrix G(J, C - 1);
    NumericVector diag_vec(static_cast<std::size_t>(J) * (C - 1));
    NumericVector off_vec(static_cast<std::size_t>(J) * (C - 2));

    const double* p_X   = REAL(X);
    const double* p_tau = REAL(tau);
    const int*    p_Y   = INTEGER(Y);

    std::vector<double> P(C), f(C - 1);
    std::vector<long double> g_acc(C - 1), d_acc(C - 1),
                             o_acc(C > 2 ? C - 2 : 0);

    for (int j = 0; j < J; ++j) {
        const double* X_col = p_X + static_cast<std::size_t>(j) * N;
        const int*    Y_col = p_Y + static_cast<std::size_t>(j) * N;

        std::fill(g_acc.begin(), g_acc.end(), 0.0L);
        std::fill(d_acc.begin(), d_acc.end(), 0.0L);
        std::fill(o_acc.begin(), o_acc.end(), 0.0L);

        for (int i = 0; i < N; ++i) {
            const double x = X_col[i];

            // P[0..C-1], f[0..C-2] for this (i, j) -- same recursion as
            // compute_PFf, including the eps floor on P.
            double F_prev = 0.0;
            for (int c = 0; c < C - 1; ++c) {
                const double z = p_tau[static_cast<std::size_t>(c + 1) * J + j] - x;
                const double F_curr = 1.0 / (1.0 + std::exp(-z));
                double Pv = F_curr - F_prev;
                if (Pv < eps) Pv = eps;
                P[c] = Pv;
                f[c] = F_curr * (1.0 - F_curr);
                F_prev = F_curr;
            }
            double P_last = 1.0 - F_prev;
            if (P_last < eps) P_last = eps;
            P[C - 1] = P_last;

            const int y = Y_col[i];

            for (int k = 0; k < C - 1; ++k) {
                const double fk = f[k];

                // grad_tau: (Y == k+1) * fk/P[k] - (Y == k+2) * fk/P[k+1]
                double gv = 0.0;
                if (y == k + 1) gv += fk / P[k];
                if (y == k + 2) gv -= fk / P[k + 1];
                g_acc[k] += gv;

                // fisher diagonal: f^2/P_left + f^2/P_right
                const double f2 = fk * fk;
                d_acc[k] += f2 / P[k] + f2 / P[k + 1];
            }

            // fisher off-diagonal: (f[m] * f[m+1]) / P[m+1]
            for (int m = 0; m < C - 2; ++m) {
                o_acc[m] += (f[m] * f[m + 1]) / P[m + 1];
            }
        }

        for (int k = 0; k < C - 1; ++k) {
            double gv = -static_cast<double>(g_acc[k]);          // G <- -G
            if (gv >  clip_value) gv =  clip_value;              // pmin/pmax
            if (gv < -clip_value) gv = -clip_value;
            G(j, k) = gv;
            diag_vec[static_cast<std::size_t>(j) * (C - 1) + k] =
                static_cast<double>(d_acc[k]);
        }
        for (int m = 0; m < C - 2; ++m) {
            off_vec[static_cast<std::size_t>(j) * (C - 2) + m] =
                -static_cast<double>(o_acc[m]);
        }
    }

    return List::create(_["G"] = G, _["diag_vec"] = diag_vec,
                        _["off_vec"] = off_vec);
}

// ============================================================================
// compute_H_gpcm  (GPCM branch)
// ============================================================================
// [[Rcpp::export]]
List compute_H_gpcm(IntegerMatrix Y,
                        NumericMatrix eta,
                        NumericVector tau,
                        int C) {
    const int N = Y.nrow();
    const int J = Y.ncol();

    NumericMatrix Hw(N, J);
    NumericVector col_sums(J * C);

    const int*    p_Y   = INTEGER(Y);
    const double* p_eta = REAL(eta);
    const double* p_tau = REAL(tau);
    double*       p_Hw  = REAL(Hw);
    double*       p_cs  = REAL(col_sums);

    std::vector<double> Xc(C);
    std::vector<double> Pc(C);

    for (int j = 0; j < J; ++j) {
        const double* eta_col = p_eta + static_cast<std::size_t>(j) * N;
        const int*    Y_col   = p_Y   + static_cast<std::size_t>(j) * N;
        const double* tau_j   = p_tau + static_cast<std::size_t>(j) * C;
        double*       Hw_col  = p_Hw  + static_cast<std::size_t>(j) * N;
        double*       cs_j    = p_cs  + static_cast<std::size_t>(j) * C;

        for (int i = 0; i < N; ++i) {
            const double eta_ij = eta_col[i];

            double M = tau_j[0];
            Xc[0] = M;
            for (int c = 1; c < C; ++c) {
                const double X = c * eta_ij + tau_j[c];
                Xc[c] = X;
                if (X > M) M = X;
            }

            double S = 0.0;
            for (int c = 0; c < C; ++c) {
                Pc[c] = std::exp(Xc[c] - M);
                S += Pc[c];
            }

            double Hw_ij = 0.0;
            for (int c = 0; c < C; ++c) {
                const double P = Pc[c] / S;
                cs_j[c] -= P;
                Hw_ij  -= c * P;
            }

            const int y = Y_col[i];
            cs_j[y - 1] += 1.0;
            Hw_ij       += (y - 1);

            Hw_col[i] = Hw_ij;
        }
    }

    return List::create(_["Hw"] = Hw, _["col_sums"] = col_sums);
}

// ============================================================================
// compute_H_grm  (GRM branch)
// ============================================================================
// [[Rcpp::export]]
NumericMatrix compute_H_grm(IntegerMatrix Y, NumericMatrix P_all, int C) {
    const int N = Y.nrow();
    const int J = Y.ncol();

    NumericMatrix H(N, J);

    const int*    p_Y = INTEGER(Y);
    const double* p_P = REAL(P_all);
    double*       p_H = REAL(H);

    for (int j = 0; j < J; ++j) {
        const int* Y_col = p_Y + static_cast<std::size_t>(j) * N;
        double*    H_col = p_H + static_cast<std::size_t>(j) * N;

        for (int c = 0; c < C; ++c) {
            const double* P_col = p_P + static_cast<std::size_t>(j * C + c) * N;
            const int c_one = c + 1;

            for (int i = 0; i < N; ++i) {
                const int y = Y_col[i];
                const int sgn = (c_one > y) - (c_one < y);
                H_col[i] += sgn * P_col[i];
            }
        }
    }

    return H;
}

// ============================================================================
// compute_Z_grm_cpp  (fused GRM working response)
// ============================================================================
// Returns the MM working response Z = X - 2H in one streaming pass, i.e.
//
//   compute_H_grm(Y, compute_PFf(X, tau, need_f = FALSE)$P_all, C)
//   Z <- X - 2 * H
//
// without ever allocating the N x J*C P_all array or the N x J H and 2*H
// temporaries. Per row only C probabilities are held in scratch.
//
// Arithmetic is performed in the same order as the two-step version (P
// recursion identical to compute_PFf including the eps floor; the signed
// accumulation runs over c = 0..C-1 exactly as compute_H_grm does; then
// x - 2.0 * h), so the result is bitwise identical.
// [[Rcpp::export]]
NumericMatrix compute_Z_grm_cpp(IntegerMatrix Y, NumericMatrix X,
                                NumericMatrix tau) {
    const int N = X.nrow();
    const int J = X.ncol();
    const int C = tau.ncol() - 1;
    const double eps = std::numeric_limits<double>::epsilon();

    if (Y.nrow() != N || Y.ncol() != J) stop("Y dimensions must match X");
    if (C < 2) stop("C must be >= 2");

    NumericMatrix Z(N, J);

    const double* p_X   = REAL(X);
    const double* p_tau = REAL(tau);
    const int*    p_Y   = INTEGER(Y);
    double*       p_Z   = REAL(Z);

    std::vector<double> P(C);

    for (int j = 0; j < J; ++j) {
        const double* X_col = p_X + static_cast<std::size_t>(j) * N;
        const int*    Y_col = p_Y + static_cast<std::size_t>(j) * N;
        double*       Z_col = p_Z + static_cast<std::size_t>(j) * N;

        for (int i = 0; i < N; ++i) {
            const double x = X_col[i];

            double F_prev = 0.0;
            for (int c = 0; c < C - 1; ++c) {
                const double z = p_tau[static_cast<std::size_t>(c + 1) * J + j] - x;
                const double F_curr = 1.0 / (1.0 + std::exp(-z));
                double Pv = F_curr - F_prev;
                if (Pv < eps) Pv = eps;
                P[c] = Pv;
                F_prev = F_curr;
            }
            double P_last = 1.0 - F_prev;
            if (P_last < eps) P_last = eps;
            P[C - 1] = P_last;

            const int y = Y_col[i];
            double h = 0.0;
            for (int c = 0; c < C; ++c) {
                const int c_one = c + 1;
                const int sgn = (c_one > y) - (c_one < y);
                h += sgn * P[c];
            }

            Z_col[i] = x - 2.0 * h;
        }
    }

    return Z;
}

// ============================================================================
// init_M_tilde  (fused initialization transform for paramsINIT_binary)
// ============================================================================
// Replaces this block of paramsINIT_binary():
//
//   M <- suppressWarnings(
//          ifelse(Q_tilde < -1+epsilon, qlogis(epsilon),
//          ifelse(Q_tilde >  1-epsilon, qlogis(1-epsilon),
//                 qlogis(0.5*(Q_tilde+1)))))
//   d       <- apply(M, 2, sum) / N
//   M_tilde <- t(t(M) - d)
//
// R's ifelse() evaluates both branches over the whole matrix and then
// materializes another full-size copy via rep(yes, length.out = ...), and
// it is nested here -- roughly eight N x J allocations, measured at +420 MB
// (J = 500, N = 12500) and by far the largest single allocation in the
// whole GRM pipeline. t(t(M) - d) adds three more.
//
// This kernel does the same work with ONE N x J allocation: pass 1 writes M
// and accumulates the column sums, pass 2 centers in place. Column sums use
// a long double accumulator to match R's sum(); qlogis(p) is expanded to
// its nmath definition log(p / (1 - p)) (location 0, scale 1), and the
// branch order and arithmetic shapes mirror the R code, so the result is
// bitwise identical.
//
// Returns: list(M_tilde = <N x J>, d = <length J>)
// [[Rcpp::export]]
List init_M_tilde(NumericMatrix Q_tilde, double epsilon) {
    const int N = Q_tilde.nrow();
    const int J = Q_tilde.ncol();

    NumericMatrix M(N, J);
    NumericVector d(J);

    const double* p_Q = REAL(Q_tilde);
    double*       p_M = REAL(M);
    double*       p_d = REAL(d);

    // Cut points and clipped values, formed exactly as the R code does.
    const double lo_cut = -1.0 + epsilon;
    const double hi_cut =  1.0 - epsilon;
    const double lo_val = std::log(epsilon / (1.0 - epsilon));
    const double hi_val = std::log(hi_cut  / (1.0 - hi_cut));

    for (int j = 0; j < J; ++j) {
        const double* Q_col = p_Q + static_cast<std::size_t>(j) * N;
        double*       M_col = p_M + static_cast<std::size_t>(j) * N;

        long double acc = 0.0L;
        for (int i = 0; i < N; ++i) {
            const double q = Q_col[i];
            double m;
            if (q < lo_cut) {
                m = lo_val;
            } else if (q > hi_cut) {
                m = hi_val;
            } else {
                const double u = 0.5 * (q + 1.0);
                m = std::log(u / (1.0 - u));
            }
            M_col[i] = m;
            acc += m;
        }

        // d <- apply(M, 2, sum) / N   (sum() reduces in long double)
        const double dj = static_cast<double>(acc) / N;
        p_d[j] = dj;

        // M_tilde <- t(t(M) - d), done in place: M[i, j] - d[j]
        for (int i = 0; i < N; ++i) M_col[i] -= dj;
    }

    return List::create(_["M_tilde"] = M, _["d"] = d);
}

// ============================================================================
// logl_grm  (GRM branch)
// ============================================================================
// [[Rcpp::export]]
double logl_grm(IntegerMatrix Y, NumericMatrix tau,
                    Nullable<NumericMatrix> X_ = R_NilValue,
                    Nullable<NumericMatrix> P_all_ = R_NilValue) {
    const int N = Y.nrow();
    const int J = Y.ncol();
    const int C = tau.ncol() - 1;

    const int*    p_Y   = INTEGER(Y);
    const double* p_tau = REAL(tau);

    double logl_value = 0.0;
    const double eps = std::numeric_limits<double>::epsilon();

    if (P_all_.isNotNull()) {
        NumericMatrix P_all(P_all_);
        const double* p_P = REAL(P_all);
        const int P_C = P_all.ncol() / J;
        if (P_C != C) stop("P_all ncol mismatch with tau");

        for (int j = 0; j < J; ++j) {
            const int* Y_col = p_Y + static_cast<std::size_t>(j) * N;
            for (int i = 0; i < N; ++i) {
                const int y = Y_col[i];
                const int col_idx = j * C + (y - 1);
                double P = p_P[static_cast<std::size_t>(col_idx) * N + i];
                if (P < eps) P = eps;
                logl_value -= std::log(P);
            }
        }
    } else if (X_.isNotNull()) {
        NumericMatrix X(X_);
        const double* p_X = REAL(X);

        for (int j = 0; j < J; ++j) {
            const double* X_col   = p_X   + static_cast<std::size_t>(j) * N;
            const int*    Y_col   = p_Y   + static_cast<std::size_t>(j) * N;

            for (int i = 0; i < N; ++i) {
                const int y = Y_col[i];          // 1-indexed category
                const double z_lo = p_tau[(y - 1) * J + j] - X_col[i];
                const double z_hi = p_tau[y * J + j]       - X_col[i];
                const double F_lo = 1.0 / (1.0 + std::exp(-z_lo));
                const double F_hi = 1.0 / (1.0 + std::exp(-z_hi));
                double P = F_hi - F_lo;
                if (P < eps) P = eps;
                logl_value -= std::log(P);
            }
        }
    } else {
        stop("Either X or P_all must be provided for GRM logl");
    }

    return logl_value;
}

// ============================================================================
// logl_gpcm  (GPCM branch)
// ============================================================================
// [[Rcpp::export]]
double logl_gpcm(IntegerMatrix Y, NumericMatrix eta,
                     NumericVector tau, int C) {
    const int N = Y.nrow();
    const int J = Y.ncol();

    const int*    p_Y   = INTEGER(Y);
    const double* p_eta = REAL(eta);
    const double* p_tau = REAL(tau);

    // log_num_sum = sum((Y - 1) * eta) + sum(tau[tau_idx])
    double log_num_sum = 0.0;
    for (int j = 0; j < J; ++j) {
        const int*    Y_col   = p_Y   + static_cast<std::size_t>(j) * N;
        const double* eta_col = p_eta + static_cast<std::size_t>(j) * N;
        for (int i = 0; i < N; ++i) {
            const int y = Y_col[i] - 1;  // 0-indexed category
            log_num_sum += y * eta_col[i];
            log_num_sum += p_tau[static_cast<std::size_t>(j) * C + y];
        }
    }

    // log-denominator via streaming log-sum-exp over the C axis
    double log_den = 0.0;
    for (int j = 0; j < J; ++j) {
        const double* eta_col = p_eta + static_cast<std::size_t>(j) * N;
        for (int i = 0; i < N; ++i) {
            double M = p_tau[static_cast<std::size_t>(j) * C + 0];
            double S = 1.0;
            for (int c = 1; c < C; ++c) {
                const double Xc = c * eta_col[i] +
                                  p_tau[static_cast<std::size_t>(j) * C + c];
                const double M_new = std::max(M, Xc);
                S = S * std::exp(M - M_new) + std::exp(Xc - M_new);
                M = M_new;
            }
            log_den += std::log(S) + M;
        }
    }

    return log_den - log_num_sum;
}
