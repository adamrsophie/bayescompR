# Linear algebra with reverse mode adjoints.
#
# Adjoint formulas follow Giles (2008), "Collected matrix derivative results
# for forward and reverse mode algorithmic differentiation", and Murray
# (2016), "Differentiation of the Cholesky decomposition". The adjoints are
# checked against central finite differences in tests/testthat/test-ad.R.

sym_check <- function(f, A, tol = 1e-8) {
  v <- val(A)
  if (nrow(v) != ncol(v)) stop_domain(sprintf("%s: matrix is not square (%d x %d)", f, nrow(v), ncol(v)))
  if (any(!is.finite(v))) stop_domain(sprintf("%s: matrix has non-finite elements", f))
  if (max(abs(v - t(v)), 0) > tol * max(1, max(abs(v)))) stop_domain(sprintf("%s: matrix is not symmetric", f))
}
square_check <- function(f, A) {
  v <- val(A)
  if (nrow(v) != ncol(v)) stop_domain(sprintf("%s: matrix is not square (%d x %d)", f, nrow(v), ncol(v)))
}

ad_cholesky <- function(A) {
  A <- as_dbl(A)
  sym_check("cholesky_decompose", A)
  vA <- val(A)
  R <- tryCatch(chol(vA), error = function(e) NULL)
  if (is.null(R)) stop_domain("cholesky_decompose: matrix is not positive definite")
  L <- t(R)
  if (!is_ad(A)) return(L)
  ad_node(L, list(A), function(Lbar) {
    Lbar[upper.tri(Lbar)] <- 0
    P <- crossprod(L, Lbar)
    P[upper.tri(P)] <- 0
    diag(P) <- diag(P) / 2
    Linv <- backsolve(R, diag(nrow(L)))  # R^{-1} = L^{-T}
    S <- Linv %*% P %*% t(Linv)
    list(0.5 * (S + t(S)))
  })
}
def("cholesky_decompose", "matrix", "matrix", ad_cholesky)

ad_inverse <- function(A) {
  A <- as_dbl(A)
  square_check("inverse", A)
  Y <- tryCatch(solve(val(A)), error = function(e) stop_domain("inverse: matrix is singular"))
  if (!is_ad(A)) return(Y)
  ad_node(Y, list(A), function(g) list(-t(Y) %*% g %*% t(Y)))
}
def("inverse", "matrix", "matrix", ad_inverse)
def("inverse_spd", "matrix", "matrix", function(A) {
  sym_check("inverse_spd", A)
  if (!is_pos_def(val(A))) stop_domain("inverse_spd: matrix is not positive definite")
  ad_inverse(A)
})
def("chol2inv", "matrix", "matrix", function(L) {
  L <- as_dbl(L)
  ad_inverse(ad_matmul(L, ad_transpose(L)))
})

ad_log_det <- function(A) {
  A <- as_dbl(A)
  square_check("log_determinant", A)
  vA <- val(A)
  d <- determinant(vA, logarithm = TRUE)
  y <- as.numeric(d$modulus)
  if (!is_ad(A)) return(y)
  ad_node(y, list(A), function(g) list(g * t(solve(vA))))
}
def("log_determinant", "matrix", "real", ad_log_det)
def("log_determinant_spd", "matrix", "real", function(A) {
  sym_check("log_determinant_spd", A)
  L <- ad_cholesky(A)
  2 * ad_sum(log(.lib$diagonal[[1L]]$impl(L)))
})
def("determinant", "matrix", "real", function(A) {
  A <- as_dbl(A)
  square_check("determinant", A)
  vA <- val(A)
  y <- det(vA)
  if (!is_ad(A)) return(y)
  ad_node(y, list(A), function(g) list(g * y * t(solve(vA))))
})

# solve A x = B
ad_solve <- function(A, B, f = "mdivide_left") {
  A <- as_dbl(A)
  B <- as_dbl(B)
  vA <- val(A)
  vB <- val(B)
  vec <- is.null(dim(vB))
  Bm <- if (vec) matrix(vB, ncol = 1L) else vB
  square_check(f, A)
  if (nrow(vA) != nrow(Bm)) stop_shape(sprintf("%s: rows of the matrix (%d) must match rows of the right side (%d)", f, nrow(vA), nrow(Bm)))
  X <- tryCatch(solve(vA, Bm), error = function(e) stop_domain(paste0(f, ": matrix is singular")))
  out <- if (vec) as.vector(X) else X
  if (!any_ad(A, B)) return(out)
  ad_node(out, list(A, B), function(g) {
    G <- if (vec) matrix(g, ncol = 1L) else g
    Bbar <- solve(t(vA), G)
    list(
      if (is_ad(A)) -Bbar %*% t(X),
      if (is_ad(B)) (if (vec) as.vector(Bbar) else Bbar))
  })
}
def("mdivide_left", c("matrix", "vector"), "vector", function(A, b) ad_solve(A, b))
def("mdivide_left", c("matrix", "matrix"), "matrix", function(A, B) ad_solve(A, B))
def("mdivide_left_spd", c("matrix", "vector"), "vector", function(A, b) {
  sym_check("mdivide_left_spd", A); ad_solve(A, b, "mdivide_left_spd")
})
def("mdivide_left_spd", c("matrix", "matrix"), "matrix", function(A, B) {
  sym_check("mdivide_left_spd", A); ad_solve(A, B, "mdivide_left_spd")
})
# B / A = B * inverse(A) = t(solve(t(A), t(B)))
mdivide_right_impl <- function(B, A, rowvec = FALSE) {
  B <- as_dbl(B)
  Bm <- if (rowvec) as_row(B) else B
  X <- ad_transpose(ad_solve(ad_transpose(A), ad_transpose(Bm), "mdivide_right"))
  if (rowvec) ad_reshape(X, NULL) else X
}
def("mdivide_right", c("row_vector", "matrix"), "row_vector", function(b, A) mdivide_right_impl(b, A, TRUE))
def("mdivide_right", c("matrix", "matrix"), "matrix", function(B, A) mdivide_right_impl(B, A))
def("mdivide_right_spd", c("row_vector", "matrix"), "row_vector", function(b, A) mdivide_right_impl(b, A, TRUE))
def("mdivide_right_spd", c("matrix", "matrix"), "matrix", function(B, A) mdivide_right_impl(B, A))

# triangular solves use only the lower triangle of L
ad_tri_solve <- function(L, B, f) {
  L <- as_dbl(L)
  B <- as_dbl(B)
  vL <- val(L)
  vL[upper.tri(vL)] <- 0
  vB <- val(B)
  vec <- is.null(dim(vB))
  Bm <- if (vec) matrix(vB, ncol = 1L) else vB
  square_check(f, L)
  if (nrow(vL) != nrow(Bm)) stop_shape(sprintf("%s: size mismatch", f))
  if (any(diag(vL) == 0)) stop_domain(paste0(f, ": matrix is singular"))
  X <- forwardsolve(vL, Bm)
  out <- if (vec) as.vector(X) else X
  if (!any_ad(L, B)) return(out)
  ad_node(out, list(L, B), function(g) {
    G <- if (vec) matrix(g, ncol = 1L) else g
    Bbar <- backsolve(t(vL), G)
    gL <- -Bbar %*% t(X)
    gL[upper.tri(gL)] <- 0
    list(if (is_ad(L)) gL, if (is_ad(B)) (if (vec) as.vector(Bbar) else Bbar))
  })
}
def("mdivide_left_tri_low", c("matrix", "vector"), "vector", function(L, b) ad_tri_solve(L, b, "mdivide_left_tri_low"))
def("mdivide_left_tri_low", c("matrix", "matrix"), "matrix", function(L, B) ad_tri_solve(L, B, "mdivide_left_tri_low"))
# B * inverse(tri(L)): X = t(solve(t(L), t(B))) with t(L) upper triangular
ad_tri_solve_right <- function(B, L, rowvec) {
  L <- as_dbl(L)
  B <- as_dbl(B)
  Bm <- if (rowvec) as_row(B) else B
  vL <- val(L)
  mask <- matrix(1, nrow(vL), ncol(vL))
  mask[upper.tri(mask)] <- 0
  Lt <- ad_transpose(L * mask)
  X <- ad_transpose(ad_solve(Lt, ad_transpose(Bm), "mdivide_right_tri_low"))
  if (rowvec) ad_reshape(X, NULL) else X
}
def("mdivide_right_tri_low", c("row_vector", "matrix"), "row_vector", function(b, L) ad_tri_solve_right(b, L, TRUE))
def("mdivide_right_tri_low", c("matrix", "matrix"), "matrix", function(B, L) ad_tri_solve_right(B, L, FALSE))

def("crossprod", "matrix", "matrix", function(A) ad_matmul(ad_transpose(as_dbl(A)), as_dbl(A)))
def("tcrossprod", "matrix", "matrix", function(A) ad_matmul(as_dbl(A), ad_transpose(as_dbl(A))))
def("multiply_lower_tri_self_transpose", "matrix", "matrix", function(A) {
  A <- as_dbl(A)
  v <- val(A)
  mask <- matrix(1, nrow(v), ncol(v))
  mask[upper.tri(mask)] <- 0
  L <- A * mask
  ad_matmul(L, ad_transpose(L))
})
def("diag_pre_multiply", c("vector", "matrix"), "matrix", function(v, M) diag_pre(v, M))
def("diag_pre_multiply", c("row_vector", "matrix"), "matrix", function(v, M) diag_pre(v, M))
def("diag_post_multiply", c("matrix", "vector"), "matrix", function(M, v) diag_post(M, v))
def("diag_post_multiply", c("matrix", "row_vector"), "matrix", function(M, v) diag_post(M, v))
diag_pre <- function(v, M) {
  d <- dim(val(M))
  if (nelem(v) != d[1L]) stop_shape("diag_pre_multiply: size mismatch")
  as_dbl(M) * ad_gather(as_dbl(v), rep(seq_len(d[1L]), d[2L]), d)
}
diag_post <- function(M, v) {
  d <- dim(val(M))
  if (nelem(v) != d[2L]) stop_shape("diag_post_multiply: size mismatch")
  as_dbl(M) * ad_gather(as_dbl(v), rep(seq_len(d[2L]), each = d[1L]), d)
}
def("quad_form", c("matrix", "matrix"), "matrix", function(A, B) {
  B <- as_dbl(B)
  ad_matmul(ad_matmul(ad_transpose(B), as_dbl(A)), B)
})
def("quad_form", c("matrix", "vector"), "real", function(A, v) {
  v <- as_dbl(v)
  ad_sum(v * op_mat_vec(A, v))
})
def("quad_form_sym", c("matrix", "matrix"), "matrix", function(A, B) {
  sym_check("quad_form_sym", A)
  B <- as_dbl(B)
  R <- ad_matmul(ad_matmul(ad_transpose(B), as_dbl(A)), B)
  0.5 * (R + ad_transpose(R))
})
def("quad_form_sym", c("matrix", "vector"), "real", function(A, v) {
  sym_check("quad_form_sym", A)
  v <- as_dbl(v)
  ad_sum(v * op_mat_vec(A, v))
})
def("quad_form_diag", c("matrix", "vector"), "matrix", function(M, v) diag_post(diag_pre(v, M), v))
def("quad_form_diag", c("matrix", "row_vector"), "matrix", function(M, v) diag_post(diag_pre(v, M), v))
def("trace", "matrix", "real", function(A) ad_sum(.lib$diagonal[[1L]]$impl(A)))
def("trace_quad_form", c("matrix", "matrix"), "real", function(A, B) {
  B <- as_dbl(B)
  ad_sum(B * ad_matmul(as_dbl(A), B))
})
def("trace_quad_form", c("matrix", "vector"), "real", function(A, b) {
  b <- as_dbl(b)
  if (nrow(val(A)) != ncol(val(A)) || ncol(val(A)) != nelem(b)) stop_shape("trace_quad_form: size mismatch")
  ad_sum(b * op_mat_vec(A, b))
})
# trace(A * B)
def("trace_dot", c("matrix", "matrix"), "real", function(A, B) {
  dA <- dim(val(A)); dB <- dim(val(B))
  if (dA[2L] != dB[1L] || dA[1L] != dB[2L]) stop_shape("trace_dot: size mismatch")
  ad_sum(as_dbl(A) * ad_transpose(as_dbl(B)))
})
def("symmetrize_from_lower_tri", "matrix", "matrix", function(A) {
  d <- dim(val(A))
  if (d[1L] != d[2L]) stop_shape("symmetrize_from_lower_tri: the matrix must be square")
  n <- d[1L]
  i <- rep(seq_len(n), n); j <- rep(seq_len(n), each = n)
  # entry (i, j) above the diagonal takes (j, i)
  src <- ifelse(i >= j, i + (j - 1L) * n, j + (i - 1L) * n)
  if (is_ad(A)) ad_gather(A, src, d) else { out <- val(A)[src]; dim(out) <- d; out }
})
def("trace_gen_quad_form", c("matrix", "matrix", "matrix"), "real", function(D, A, B) {
  B <- as_dbl(B)
  ad_sum(.lib$diagonal[[1L]]$impl(ad_matmul(ad_matmul(as_dbl(D), ad_transpose(B)), ad_matmul(as_dbl(A), B))))
})
def("matrix_power", c("matrix", "int"), "matrix", function(A, n) {
  square_check("matrix_power", A)
  if (n < 0L) stop_domain("matrix_power: exponent must be nonnegative")
  k <- nrow(val(A))
  R <- diag(1, k)
  if (n == 0L) return(R)
  R <- as_dbl(A)
  if (n > 1L) for (i in 2:n) R <- ad_matmul(R, as_dbl(A))
  R
})

ad_eigen_sym <- function(A, want) {
  A <- as_dbl(A)
  sym_check(if (want == "values") "eigenvalues_sym" else "eigenvectors_sym", A)
  vA <- val(A)
  e <- eigen((vA + t(vA)) / 2, symmetric = TRUE)
  o <- order(e$values)  # Stan returns ascending order
  lam <- e$values[o]
  V <- e$vectors[, o, drop = FALSE]
  if (want == "values") {
    if (!is_ad(A)) return(lam)
    return(ad_node(lam, list(A), function(g) list(V %*% (g * t(V)))))
  }
  if (!is_ad(A)) return(V)
  ad_node(V, list(A), function(G) {
    k <- length(lam)
    Fm <- outer(lam, lam, function(a, b) b - a)
    Fm <- ifelse(abs(Fm) > 0, 1 / Fm, 0)
    diag(Fm) <- 0
    S <- V %*% (Fm * (t(V) %*% G)) %*% t(V)
    list(0.5 * (S + t(S)))
  })
}
def("eigenvalues_sym", "matrix", "vector", function(A) ad_eigen_sym(A, "values"))
def("eigenvectors_sym", "matrix", "matrix", function(A) ad_eigen_sym(A, "vectors"))
def("eigenvalues", "matrix", "vector", function(A) stop_unsupported("eigenvalues of a nonsymmetric matrix return complex values; complex numbers are not supported"))

def("singular_values", "matrix", "vector", function(A) {
  A <- as_dbl(A)
  s <- svd(val(A))
  if (!is_ad(A)) return(s$d)
  ad_node(s$d, list(A), function(g) list(s$u %*% (g * t(s$v))))
})
def("svd_U", "matrix", "matrix", function(A) {
  A <- as_dbl(A)
  s <- svd(val(A))
  if (is_ad(A)) stop_unsupported("svd_U with a parameter dependent argument has no gradient in this version")
  s$u
})
def("svd_V", "matrix", "matrix", function(A) {
  A <- as_dbl(A)
  s <- svd(val(A))
  if (is_ad(A)) stop_unsupported("svd_V with a parameter dependent argument has no gradient in this version")
  s$v
})

# matrix exponential through Matrix::expm; the adjoint uses the block matrix
# identity expm([[X^T, G], [0, X^T]]) whose upper right block is the
# Frechet derivative.
mexp <- function(X) as.matrix(Matrix::expm(Matrix::Matrix(X, sparse = FALSE)))
ad_matrix_exp <- function(A) {
  A <- as_dbl(A)
  square_check("matrix_exp", A)
  vA <- val(A)
  E <- mexp(vA)
  if (!is_ad(A)) return(E)
  ad_node(E, list(A), function(G) {
    k <- nrow(vA)
    Bm <- rbind(cbind(t(vA), G), cbind(matrix(0, k, k), t(vA)))
    list(mexp(Bm)[seq_len(k), k + seq_len(k)])
  })
}
def("matrix_exp", "matrix", "matrix", ad_matrix_exp)
def("matrix_exp_multiply", c("matrix", "matrix"), "matrix", function(A, B) ad_matmul(ad_matrix_exp(A), as_dbl(B)))
def("scale_matrix_exp_multiply", c("real", "matrix", "matrix"), "matrix", function(t, A, B) ad_matmul(ad_matrix_exp(as_dbl(t) * as_dbl(A)), as_dbl(B)))

# thin QR by modified Gram Schmidt written with differentiable operations,
# with Stan's sign convention (R has a nonnegative diagonal)
qr_thin <- function(A) {
  A <- as_dbl(A)
  v <- val(A)
  m <- nrow(v)
  n <- ncol(v)
  k <- min(m, n)
  cols <- lapply(seq_len(n), function(j) ad_gather(A, (j - 1L) * m + seq_len(m)))
  Q <- vector("list", k)
  Rentries <- list()
  rr <- integer(0)
  cc <- integer(0)
  for (j in seq_len(n)) {
    w <- cols[[j]]
    for (i in seq_len(min(j - 1L, k))) {
      rij <- ad_sum(Q[[i]] * w)
      Rentries[[length(Rentries) + 1L]] <- rij
      rr <- c(rr, i); cc <- c(cc, j)
      w <- w - rij * Q[[i]]
    }
    if (j <= k) {
      rjj <- sqrt(ad_sum(w * w))
      if (val(rjj) == 0) stop_domain("qr: matrix is rank deficient")
      Rentries[[length(Rentries) + 1L]] <- rjj
      rr <- c(rr, j); cc <- c(cc, j)
      Q[[j]] <- w / rjj
    }
  }
  Qm <- ad_reshape(ad_c(Q), c(m, k))
  Rm <- ad_scatter(matrix(0, k, n), (cc - 1L) * k + rr, ad_c(Rentries))
  list(Q = Qm, R = Rm)
}
def("qr_thin_Q", "matrix", "matrix", function(A) qr_thin(A)$Q)
def("qr_thin_R", "matrix", "matrix", function(A) qr_thin(A)$R)
def("qr_Q", "matrix", "matrix", function(A) {
  if (is_ad(A)) {
    d <- dim(val(A))
    if (d[1L] > d[2L]) stop_unsupported("qr_Q of a tall parameter dependent matrix is not differentiable in this version; use qr_thin_Q")
    return(qr_thin(A)$Q)
  }
  q <- qr(val(A))
  Q <- qr.Q(q, complete = TRUE)
  R <- qr.R(q, complete = TRUE)
  s <- sign(diag(R))
  s[s == 0] <- 1
  Q[, seq_along(s)] <- sweep(Q[, seq_along(s), drop = FALSE], 2L, s, `*`)
  Q
})
def("qr_R", "matrix", "matrix", function(A) {
  if (is_ad(A)) {
    d <- dim(val(A))
    if (d[1L] > d[2L]) stop_unsupported("qr_R of a tall parameter dependent matrix is not differentiable in this version; use qr_thin_R")
    return(qr_thin(A)$R)
  }
  q <- qr(val(A))
  R <- qr.R(q, complete = TRUE)
  s <- sign(diag(R))
  s[s == 0] <- 1
  R[seq_along(s), ] <- R[seq_along(s), , drop = FALSE] * s
  R
})
def("generalized_inverse", "matrix", "matrix", function(A) {
  A <- as_dbl(A)
  d <- dim(val(A))
  if (d[1L] == d[2L]) {
    r <- tryCatch(ad_inverse(A), bayescomp_domain_error = function(e) NULL)
    if (!is.null(r)) return(r)
  }
  if (is_ad(A)) {
    if (d[1L] >= d[2L]) return(ad_matmul(ad_inverse(ad_matmul(ad_transpose(A), A)), ad_transpose(A)))
    return(ad_matmul(ad_transpose(A), ad_inverse(ad_matmul(A, ad_transpose(A)))))
  }
  s <- svd(val(A))
  tol <- max(d) * max(s$d) * .Machine$double.eps
  dinv <- ifelse(s$d > tol, 1 / s$d, 0)
  s$v %*% (dinv * t(s$u))
})

# ---- Gaussian process covariance functions ---------------------------------------

gp_points <- function(x) {
  # x is array[] real (numeric vector) or array[] vector (matrix N x D)
  v <- val(x)
  if (is.null(dim(v))) return(list(n = length(v), d = 1L, rows = function(i) ad_select(x, i)))
  n <- nrow(v)
  D <- ncol(v)
  list(n = n, d = D, rows = function(i) ad_gather(x, i + (seq_len(D) - 1L) * n))
}
gp_sqdist <- function(x1, x2) {
  a <- gp_points(as_dbl(x1))
  b <- gp_points(as_dbl(x2))
  if (a$d != b$d) stop_shape("gp covariance: points have different dimensions")
  va <- val(x1); vb <- val(x2)
  if (!any_ad(x1, x2)) {
    A <- if (is.null(dim(va))) matrix(va, ncol = 1L) else va
    B <- if (is.null(dim(vb))) matrix(vb, ncol = 1L) else vb
    return(pmax(outer(rowSums(A^2), rowSums(B^2), `+`) - 2 * A %*% t(B), 0))
  }
  A <- if (is.null(dim(va))) ad_reshape(as_dbl(x1), c(length(va), 1L)) else as_dbl(x1)
  B <- if (is.null(dim(vb))) ad_reshape(as_dbl(x2), c(length(vb), 1L)) else as_dbl(x2)
  na <- a$n; nb <- b$n; D <- a$d
  out <- 0
  for (k in seq_len(D)) {
    ak <- ad_gather(A, (k - 1L) * na + seq_len(na))
    bk <- ad_gather(B, (k - 1L) * nb + seq_len(nb))
    diffm <- ad_gather(ak, rep(seq_len(na), nb), c(na, nb)) - ad_gather(bk, rep(seq_len(nb), each = na), c(na, nb))
    out <- out + diffm * diffm
  }
  out
}
gp_kernel <- function(name, kfun) {
  for (xt in c("real[]", "vector[]")) {
    local({
      x_t <- xt
      def(name, c(x_t, "real", "real"), "matrix", function(x, s, l) {
        chk_positive(name, "magnitude", s); chk_positive(name, "length scale", l)
        kfun(gp_sqdist(x, x), as_dbl(s), as_dbl(l))
      })
      def(name, c(x_t, x_t, "real", "real"), "matrix", function(x1, x2, s, l) {
        chk_positive(name, "magnitude", s); chk_positive(name, "length scale", l)
        kfun(gp_sqdist(x1, x2), as_dbl(s), as_dbl(l))
      })
    })
  }
}
# distance from squared distance; the derivative at zero distance is taken as
# zero (the diagonal and coincident inputs) instead of infinity
gp_dist <- function(d2) {
  v <- val(d2)
  r <- sqrt(v)
  if (!is_ad(d2)) return(r)
  ad_node(r, list(d2), function(g) list(ifelse(r > 0, g * 0.5 / r, 0)))
}
k_sqexp <- function(d2, s, l) s * s * exp(-0.5 * d2 / (l * l))
k_exp <- function(d2, s, l) s * s * exp(-gp_dist(d2) / l)
k_m32 <- function(d2, s, l) {
  r <- gp_dist(d2) * sqrt(3) / l
  s * s * (1 + r) * exp(-r)
}
k_m52 <- function(d2, s, l) {
  r <- gp_dist(d2) * sqrt(5) / l
  s * s * (1 + r + r * r / 3) * exp(-r)
}
# automatic relevance determination: one length scale per input dimension,
# applied by scaling the inputs and using a unit length scale
gp_ard_scale <- function(name, x, l) {
  v <- val(x)
  chk_positive(name, "length scale", l)
  if (ncol(v) != nelem(l)) stop_shape(sprintf("%s: the number of length scales must equal the dimension of the inputs", name))
  as_dbl(x) / ad_gather(as_dbl(l), rep(seq_len(ncol(v)), each = nrow(v)), dim(v))
}
gp_kernel_ard <- function(name, kfun) {
  def(name, c("vector[]", "real", "real[]"), "matrix", function(x, s, l) {
    chk_positive(name, "magnitude", s)
    xs <- gp_ard_scale(name, x, l)
    kfun(gp_sqdist(xs, xs), as_dbl(s), 1)
  })
  def(name, c("vector[]", "vector[]", "real", "real[]"), "matrix", function(x1, x2, s, l) {
    chk_positive(name, "magnitude", s)
    kfun(gp_sqdist(gp_ard_scale(name, x1, l), gp_ard_scale(name, x2, l)), as_dbl(s), 1)
  })
}
gp_kernel_ard("gp_exp_quad_cov", k_sqexp)
gp_kernel_ard("gp_exponential_cov", k_exp)
gp_kernel_ard("gp_matern32_cov", k_m32)
gp_kernel_ard("gp_matern52_cov", k_m52)
gp_kernel("gp_exp_quad_cov", k_sqexp)
gp_kernel("gp_exponential_cov", k_exp)
gp_kernel("gp_matern32_cov", k_m32)
gp_kernel("gp_matern52_cov", k_m52)
k_periodic <- function(d2, s, l, p) {
  chk_positive("gp_periodic_cov", "signal standard deviation", s); chk_positive("gp_periodic_cov", "length scale", l)
  chk_positive("gp_periodic_cov", "period", p)
  d <- gp_dist(d2)
  s <- as_dbl(s); l <- as_dbl(l); p <- as_dbl(p)
  sn <- sin(pi * d / p)
  s * s * exp(-2 * sn * sn / (l * l))
}
# points as an n x d matrix
gp_matrix <- function(x) {
  v <- val(x)
  if (is.null(dim(v))) ad_reshape(as_dbl(x), c(length(v), 1L)) else as_dbl(x)
}
k_dot_prod <- function(x1, x2, sigma) {
  chk_nonnegative("gp_dot_prod_cov", "sigma", sigma)
  A <- gp_matrix(x1); B <- gp_matrix(x2)
  if (ncol(val(A)) != ncol(val(B))) stop_shape("gp_dot_prod_cov: points have different dimensions")
  as_dbl(sigma) * as_dbl(sigma) + ad_matmul(A, ad_transpose(B))
}
for (xt in c("real[]", "vector[]")) local({
  x_t <- xt
  def("gp_periodic_cov", c(x_t, "real", "real", "real"), "matrix", function(x, s, l, p) k_periodic(gp_sqdist(x, x), s, l, p))
  def("gp_periodic_cov", c(x_t, x_t, "real", "real", "real"), "matrix", function(x1, x2, s, l, p) k_periodic(gp_sqdist(x1, x2), s, l, p))
  def("gp_dot_prod_cov", c(x_t, "real"), "matrix", function(x, sigma) k_dot_prod(x, x, sigma))
  def("gp_dot_prod_cov", c(x_t, x_t, "real"), "matrix", function(x1, x2, sigma) k_dot_prod(x1, x2, sigma))
})

# ---- compressed sparse row (CSR) matrices -----------------------------------------
# Stan stores the nonzeros in row major order: w the values, v their 1 based
# column indices, u the 1 based start of each row in w (length rows + 1).
csr_nonzero <- function(a) {
  va <- val(a)
  m <- nrow(va)
  # row major positions of the nonzeros, as column major indices into a
  idx <- which(t(va) != 0)
  cols <- (idx - 1L) %% ncol(va) + 1L
  rows <- (idx - 1L) %/% ncol(va) + 1L
  list(idx = (cols - 1L) * m + rows, v = as.integer(cols), u = as.integer(c(1L, cumsum(tabulate(rows, m)) + 1L)))
}
def("csr_extract_w", "matrix", "vector", function(a) {
  nz <- csr_nonzero(a)
  if (!length(nz$idx)) return(numeric(0))
  ad_gather(flat(as_dbl(a)), nz$idx)
})
def("csr_extract_v", "matrix", "int[]", function(a) csr_nonzero(a)$v)
def("csr_extract_u", "matrix", "int[]", function(a) csr_nonzero(a)$u)
csr_check <- function(f, m, n, w, v, u) {
  if (m < 0L || n < 0L) stop_domain(sprintf("%s: dimensions must be nonnegative", f))
  if (length(u) != m + 1L) stop_shape(sprintf("%s: u must have size m + 1 = %d, got %d", f, m + 1L, length(u)))
  if (length(v) != nelem(w)) stop_shape(sprintf("%s: w and v must have the same size", f))
  if (u[1L] != 1L || any(diff(u) < 0L) || u[m + 1L] - 1L != length(v))
    stop_domain(sprintf("%s: u is not a valid row start array for %d nonzeros", f, length(v)))
  if (length(v) && (any(v < 1L) || any(v > n))) stop_domain(sprintf("%s: column indices must be in 1 to %d", f, n))
  rep.int(seq_len(m), diff(u))
}
def("csr_to_dense_matrix", c("int", "int", "vector", "int[]", "int[]"), "matrix", function(m, n, w, v, u) {
  rows <- csr_check("csr_to_dense_matrix", m, n, w, v, u)
  out <- matrix(0, m, n)
  if (!length(rows)) return(out)
  pos <- (v - 1L) * m + rows
  if (anyDuplicated(pos)) stop_domain("csr_to_dense_matrix: duplicate entries")
  w <- flat(as_dbl(w))
  vw <- val(w)
  out[pos] <- vw
  if (!is_ad(w)) return(out)
  ad_node(out, list(w), function(g) list(g[pos]))
})
def("csr_matrix_times_vector", c("int", "int", "vector", "int[]", "int[]", "vector"), "vector", function(m, n, w, v, u, b) {
  f <- "csr_matrix_times_vector"
  rows <- csr_check(f, m, n, w, v, u)
  if (nelem(b) != n) stop_shape(sprintf("%s: b must have size n = %d, got %d", f, n, nelem(b)))
  w <- flat(as_dbl(w)); b <- flat(as_dbl(b))
  vw <- val(w); vb <- val(b)
  y <- if (length(rows)) tabulate_sum(rows, vw * vb[v], m) else numeric(m)
  if (!any_ad(w, b)) return(y)
  ad_node(y, list(w, b), function(g) {
    gr <- g[rows]
    list(if (is_ad(w)) gr * vb[v],
         if (is_ad(b)) if (length(v)) tabulate_sum(v, gr * vw, n) else numeric(n))
  })
})
tabulate_sum <- function(idx, x, n) {
  out <- numeric(n)
  s <- rowsum(x, idx)
  out[as.integer(rownames(s))] <- s[, 1L]
  out
}
