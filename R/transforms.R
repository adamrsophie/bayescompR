# Constraining transforms for parameters.
#
# Every constrained parameter type maps an unconstrained real vector `u` to a
# value `x` in the constrained space, together with the log absolute
# determinant of the Jacobian of that map. The definitions follow the Stan
# Reference Manual, version 2.40, chapter "Constraint Transforms". All
# functions here accept differentiable values and only use operations that
# are recorded on the AD tape.
#
# Each transform is a list with
#   free_size(shape)       number of unconstrained reals for one element
#   constrain(u, shape, a) list(x = value, lj = log Jacobian)
#   unconstrain(x, shape, a) numeric unconstrained vector
#   check(x, shape, a)     NULL when valid, otherwise a message
# where `a` holds bound or offset arguments already evaluated.

tf_identity <- list(
  free_size = function(shape) prod(shape),
  constrain = function(u, shape, a) list(x = u, lj = 0),
  unconstrain = function(x, shape, a) as.numeric(x),
  check = function(x, shape, a) NULL
)

# b + sign * exp(u) as one node
ad_exp_shift <- function(u, b, sign) {
  vu <- val(u)
  vb <- val(b)
  e <- exp(vu)
  y <- vb + sign * e
  if (!is.null(dim(vu))) dim(y) <- dim(vu)
  if (!is_ad(u) && !is_ad(b)) return(y)
  ad_node(y, list(u, b), function(g) list(if (is_ad(u)) sign * g * e, if (is_ad(b)) unbroadcast(g, vb)))
}

# lower, upper, lower and upper, offset and multiplier: elementwise
tf_bounds <- list(
  free_size = function(shape) prod(shape),
  constrain = function(u, shape, a) {
    lo <- a$lower
    up <- a$upper
    if (!is.null(a$offset) || !is.null(a$multiplier)) {
      o <- if (is.null(a$offset)) 0 else a$offset
      m <- if (is.null(a$multiplier)) 1 else a$multiplier
      if (any(!is.finite(val(m)) | val(m) <= 0)) stop_domain("multiplier must be positive and finite")
      if (any(!is.finite(val(o)))) stop_domain("offset must be finite")
      x <- o + m * u
      lm <- log(m)
      lj <- if (length(val(lm)) == 1L) lm * length(val(u)) else ad_sum(lm)
      return(list(x = x, lj = lj))
    }
    n_u <- length(val(u))
    # only lower = -Inf and upper = +Inf mean unbounded
    if (!is.null(lo) && any(val(lo) == Inf | is.nan(val(lo)))) stop_domain("lower bound must not be positive infinity or NaN")
    if (!is.null(up) && any(val(up) == -Inf | is.nan(val(up)))) stop_domain("upper bound must not be negative infinity or NaN")
    flo <- if (is.null(lo)) rep(FALSE, n_u) else rep_len(is.finite(val(lo)), n_u)
    fup <- if (is.null(up)) rep(FALSE, n_u) else rep_len(is.finite(val(up)), n_u)
    grp <- flo + 2L * fup
    if (n_u > 1L && length(unique(grp)) > 1L) {
      # finite and infinite bounds mixed elementwise: transform each group of
      # elements with the same kind of bounds, then restore the order
      pick <- function(b, idx) if (is.null(b) || length(val(b)) == 1L) b else ad_gather(b, idx)
      parts <- list(); lj <- 0; ord <- integer(0)
      for (k in unique(grp)) {
        idx <- which(grp == k)
        r <- tf_bounds$constrain(ad_gather(u, idx), length(idx),
                                 list(lower = if (k %% 2L == 1L) pick(lo, idx), upper = if (k >= 2L) pick(up, idx)))
        parts[[length(parts) + 1L]] <- r$x
        lj <- lj + r$lj
        ord <- c(ord, idx)
      }
      return(list(x = ad_gather(ad_c(parts), order(ord)), lj = lj))
    }
    has_lo <- !is.null(lo)
    has_up <- !is.null(up)
    if (has_lo && has_up) {
      vlo <- val(lo)
      vup <- val(up)
      if (all(is.finite(vlo)) && all(is.finite(vup))) {
        w <- up - lo
        x <- lo + w * ad_inv_logit(u)
        lj <- ad_sum(log(w) + ad_log_inv_logit(u) + ad_log1m_inv_logit(u))
        if (length(val(w)) == 1L && length(val(u)) > 1L)
          lj <- ad_sum(ad_log_inv_logit(u) + ad_log1m_inv_logit(u)) + length(val(u)) * log(w)
        return(list(x = x, lj = lj))
      }
      if (all(is.infinite(vlo)) && all(is.infinite(vup))) return(list(x = u, lj = 0))
      if (all(is.infinite(vup))) has_up <- FALSE
      if (all(is.infinite(vlo))) has_lo <- FALSE
      if (has_lo && has_up)
        stop_domain("bounds that are finite for some elements and infinite for others are not supported")
    }
    if (has_lo) {
      if (all(is.infinite(val(lo)))) return(list(x = u, lj = 0))
      return(list(x = ad_exp_shift(u, lo, 1), lj = ad_sum(u)))
    }
    if (has_up) {
      if (all(is.infinite(val(up)))) return(list(x = u, lj = 0))
      return(list(x = ad_exp_shift(u, up, -1), lj = ad_sum(u)))
    }
    list(x = u, lj = 0)
  },
  unconstrain = function(x, shape, a) {
    x <- as.numeric(x)
    if (!is.null(a$offset) || !is.null(a$multiplier)) {
      o <- if (is.null(a$offset)) 0 else as.numeric(a$offset)
      m <- if (is.null(a$multiplier)) 1 else as.numeric(a$multiplier)
      return((x - o) / m)
    }
    lo <- if (is.null(a$lower)) -Inf else as.numeric(a$lower)
    up <- if (is.null(a$upper)) Inf else as.numeric(a$upper)
    lo <- rep_len(lo, length(x))
    up <- rep_len(up, length(x))
    out <- x
    both <- is.finite(lo) & is.finite(up)
    onlylo <- is.finite(lo) & !is.finite(up)
    onlyup <- !is.finite(lo) & is.finite(up)
    out[both] <- stats::qlogis((x[both] - lo[both]) / (up[both] - lo[both]))
    out[onlylo] <- log(x[onlylo] - lo[onlylo])
    out[onlyup] <- log(up[onlyup] - x[onlyup])
    out
  },
  check = function(x, shape, a) {
    x <- as.numeric(x)
    if (!is.null(a$lower) && any(x < as.numeric(a$lower) | is.na(x)))
      return(sprintf("is %s, but must be greater than or equal to %s",
                     fmt_num(x[which(x < as.numeric(a$lower) | is.na(x))[1]]),
                     fmt_num(rep_len(as.numeric(a$lower), length(x))[which(x < as.numeric(a$lower) | is.na(x))[1]])))
    if (!is.null(a$upper) && any(x > as.numeric(a$upper) | is.na(x)))
      return(sprintf("is %s, but must be less than or equal to %s",
                     fmt_num(x[which(x > as.numeric(a$upper) | is.na(x))[1]]),
                     fmt_num(rep_len(as.numeric(a$upper), length(x))[which(x > as.numeric(a$upper) | is.na(x))[1]])))
    NULL
  }
)

tf_ordered <- list(
  free_size = function(shape) shape[1L],
  constrain = function(u, shape, a) {
    K <- shape[1L]
    if (K == 0L) return(list(x = u, lj = 0))
    if (K == 1L) return(list(x = u, lj = 0))
    tailu <- ad_gather(u, 2:K)
    x <- cumsum(ad_c(list(ad_gather(u, 1L), exp(tailu))))
    list(x = x, lj = ad_sum(tailu))
  },
  unconstrain = function(x, shape, a) {
    x <- as.numeric(x)
    if (length(x) <= 1L) return(x)
    c(x[1L], log(diff(x)))
  },
  check = function(x, shape, a) {
    x <- as.numeric(x)
    if (length(x) > 1L && any(!(diff(x) > 0)))
      return("is not a valid ordered vector: elements must be strictly increasing")
    NULL
  }
)

tf_positive_ordered <- list(
  free_size = function(shape) shape[1L],
  constrain = function(u, shape, a) {
    if (shape[1L] == 0L) return(list(x = u, lj = 0))
    list(x = cumsum(exp(u)), lj = ad_sum(u))
  },
  unconstrain = function(x, shape, a) {
    x <- as.numeric(x)
    if (length(x) == 0L) return(x)
    log(c(x[1L], diff(x)))
  },
  check = function(x, shape, a) {
    x <- as.numeric(x)
    if (length(x) > 0L && (x[1L] <= 0 || (length(x) > 1L && any(!(diff(x) > 0)))))
      return("is not a valid positive_ordered vector: elements must be positive and strictly increasing")
    NULL
  }
)

# Simplex transform (Stan 2.37 and later, Reference Manual 2.40, "Unit
# simplex"): the softmax of an isometric sum to zero vector. Stan Math adds
# sum(log(x)) + 0.5 * log(K) to the log density; the same constant is used
# here so that lp__ matches.
simplex_constrain <- function(u, K) {
  if (K == 1L) return(list(x = 1, lj = 0))
  z <- stz_apply(u, K)
  lse <- ad_log_sum_exp(z)
  lx <- z - lse
  list(x = exp(lx), lj = ad_sum(lx) + 0.5 * log(K))
}

simplex_free <- function(x) {
  K <- length(x)
  if (K <= 1L) return(numeric(0))
  stz_apply_t(log(x), K)
}

simplex_check <- function(x, tol = 1e-8) {
  x <- as.numeric(x)
  if (length(x) == 0L) return("is a simplex of size zero")
  if (any(is.na(x)) || any(x < 0)) return("is not a valid simplex: elements must be nonnegative")
  if (abs(sum(x) - 1) > tol) return(sprintf("is not a valid simplex: elements sum to %s, not 1", fmt_num(sum(x))))
  NULL
}

tf_simplex <- list(
  free_size = function(shape) shape[1L] - 1L,
  constrain = function(u, shape, a) simplex_constrain(u, shape[1L]),
  unconstrain = function(x, shape, a) simplex_free(as.numeric(x)),
  check = function(x, shape, a) simplex_check(x)
)

tf_unit_vector <- list(
  free_size = function(shape) shape[1L],
  constrain = function(u, shape, a) {
    ss <- ad_sum(u * u)
    if (val(ss) == 0) stop_domain("unit_vector parameter has zero norm on the unconstrained scale")
    list(x = u / sqrt(ss), lj = -0.5 * ss)
  },
  unconstrain = function(x, shape, a) as.numeric(x),
  check = function(x, shape, a) {
    if (abs(sum(as.numeric(x)^2) - 1) > 1e-8) return("is not a valid unit vector: its Euclidean norm is not 1")
    NULL
  }
)

# Isometric sum to zero transform (Stan 2.36 and later). The map is linear with
# orthonormal columns, so the log Jacobian is constant (zero). It is applied
# in linear time as Stan Math's recursion: with w_i = y_i / sqrt(i (i + 1)),
# z_k = sum_{i >= k} w_i - (k - 1) w_{k - 1}.
stz_apply <- function(u, K) {
  N <- K - 1L
  vu <- as.vector(val(u))
  i <- seq_len(N)
  sc <- 1 / sqrt(i * (i + 1))
  w <- vu * sc
  z <- c(rev(cumsum(rev(w))), 0) - c(0, w * i)
  if (!inherits(u, "bcad")) return(z)
  ad_node(z, list(u), function(g) list(stz_apply_t(g, K)))
}
# transpose of the map (also its inverse on sum to zero vectors)
stz_apply_t <- function(x, K) {
  N <- K - 1L
  i <- seq_len(N)
  x <- as.vector(x)
  (cumsum(x[i]) - i * x[i + 1L]) / sqrt(i * (i + 1))
}

.stz_cache <- new.env(parent = emptyenv())

sum_to_zero_basis <- function(K) {
  key <- as.character(K)
  if (!is.null(.stz_cache[[key]])) return(.stz_cache[[key]])
  .stz_cache[[key]] <- sum_to_zero_basis_compute(K)
}

sum_to_zero_basis_compute <- function(K) {
  N <- K - 1L
  V <- matrix(0, K, N)
  for (col in seq_len(N)) {
    y <- numeric(N)
    y[col] <- 1
    z <- numeric(K)
    sum_w <- 0
    for (i in N:1) {
      n <- i
      w <- y[i] / sqrt(n * (n + 1))
      sum_w <- sum_w + w
      z[i] <- z[i] + sum_w
      z[i + 1L] <- z[i + 1L] - w * n
    }
    V[, col] <- z
  }
  V
}

tf_sum_to_zero <- list(
  free_size = function(shape) shape[1L] - 1L,
  constrain = function(u, shape, a) {
    K <- shape[1L]
    if (K <= 1L) return(list(x = numeric(K), lj = 0))
    list(x = stz_apply(u, K), lj = 0)
  },
  unconstrain = function(x, shape, a) {
    K <- shape[1L]
    if (K <= 1L) return(numeric(0))
    stz_apply_t(x, K)
  },
  check = function(x, shape, a) {
    if (abs(sum(as.numeric(x))) > 1e-8) return("is not a valid sum_to_zero_vector: elements do not sum to zero")
    NULL
  }
)

# Sum to zero matrix (Stan 2.37): Z = B_N X t(B_M) with X the (N - 1) x
# (M - 1) free matrix. Stan Math builds the matrix with its own recursion,
# whose orthonormal basis B_K differs from the one of the vector transform:
# column c of B_K is e_c mapped by that recursion with a single column.
stz_matrix_basis <- function(K) {
  key <- paste0("m", K)
  if (!is.null(.stz_cache[[key]])) return(.stz_cache[[key]])
  N <- K - 1L
  B <- matrix(0, K, N)
  a <- 1 / sqrt((seq_len(N)) * (seq_len(N) + 1))
  b <- seq_len(N) * a
  for (c in seq_len(N)) {
    x <- numeric(N); x[c] <- 1
    z <- numeric(K)
    axp <- 0
    for (i in N:1) {
      z[i] <- b[i] * x[i] - axp
      z[K] <- z[K] - z[i]
      axp <- axp + a[i] * x[i]
    }
    B[, c] <- z
  }
  .stz_cache[[key]] <- B
}
tf_sum_to_zero_matrix <- list(
  free_size = function(shape) (shape[1L] - 1L) * (shape[2L] - 1L),
  constrain = function(u, shape, a) {
    N <- shape[1L]
    M <- shape[2L]
    if (N <= 1L || M <= 1L) return(list(x = matrix(0, N, M), lj = 0))
    X <- ad_reshape(u, c(N - 1L, M - 1L))
    Z <- ad_matmul(ad_matmul(stz_matrix_basis(N), X), t(stz_matrix_basis(M)))
    list(x = Z, lj = 0)
  },
  unconstrain = function(x, shape, a) {
    N <- shape[1L]
    M <- shape[2L]
    if (N <= 1L || M <= 1L) return(numeric(0))
    as.vector(crossprod(stz_matrix_basis(N), as.matrix(x)) %*% stz_matrix_basis(M))
  },
  check = function(x, shape, a) {
    X <- as.matrix(x)
    if (max(abs(rowSums(X)), abs(colSums(X)), 0) > 1e-8)
      return("is not a valid sum_to_zero_matrix: rows and columns must sum to zero")
    NULL
  }
)

# assemble a K x K lower triangular matrix from a list of entries indexed by
# (row, col), each possibly differentiable
assemble_matrix <- function(entries, rows, cols, K, M = K) {
  flat <- (cols - 1L) * M + rows
  vals <- ad_c(entries)
  ad_scatter(matrix(0, M, K), flat, vals)
}

# Cholesky factor of a correlation matrix (row wise canonical partial
# correlations, as in Stan's cholesky_corr_constrain).
tf_cholesky_corr <- list(
  free_size = function(shape) shape[1L] * (shape[1L] - 1L) / 2L,
  constrain = function(u, shape, a) {
    K <- shape[1L]
    if (K == 0L) return(list(x = matrix(0, 0, 0), lj = 0))
    if (K == 1L) return(list(x = matrix(1, 1, 1), lj = 0))
    z <- tanh(u)
    lj <- ad_sum(ad_log1m(z * z))
    entries <- list(1)
    rows <- 1L
    cols <- 1L
    k <- 0L
    for (i in 2:K) {
      k <- k + 1L
      zi <- ad_select(z, k)
      entries[[length(entries) + 1L]] <- zi
      rows <- c(rows, i)
      cols <- c(cols, 1L)
      sum_sqs <- zi * zi
      if (i > 2L) {
        for (j in 2:(i - 1L)) {
          lj <- lj + 0.5 * ad_log1m(sum_sqs)
          k <- k + 1L
          lij <- ad_select(z, k) * sqrt(1 - sum_sqs)
          entries[[length(entries) + 1L]] <- lij
          rows <- c(rows, i)
          cols <- c(cols, j)
          sum_sqs <- sum_sqs + lij * lij
        }
      }
      entries[[length(entries) + 1L]] <- sqrt(1 - sum_sqs)
      rows <- c(rows, i)
      cols <- c(cols, i)
    }
    list(x = assemble_matrix(entries, rows, cols, K), lj = lj)
  },
  unconstrain = function(x, shape, a) {
    L <- as.matrix(x)
    K <- nrow(L)
    if (K <= 1L) return(numeric(0))
    out <- numeric(K * (K - 1L) / 2L)
    k <- 0L
    for (i in 2:K) {
      k <- k + 1L
      out[k] <- atanh(L[i, 1L])
      sum_sqs <- L[i, 1L]^2
      if (i > 2L) {
        for (j in 2:(i - 1L)) {
          k <- k + 1L
          out[k] <- atanh(L[i, j] / sqrt(1 - sum_sqs))
          sum_sqs <- sum_sqs + L[i, j]^2
        }
      }
    }
    out
  },
  check = function(x, shape, a) {
    L <- as.matrix(x)
    if (nrow(L) != ncol(L)) return("is not a square matrix")
    if (any(abs(L[upper.tri(L)]) > 0)) return("is not lower triangular")
    if (any(diag(L) <= 0)) return("is not a valid Cholesky factor: diagonal elements must be positive")
    if (any(abs(rowSums(L^2) - 1) > 1e-8))
      return("is not a valid Cholesky factor of a correlation matrix: rows must have unit length")
    NULL
  }
)

# Correlation matrix through column wise canonical partial correlations
# (Stan's read_corr_L with the Lewandowski, Kurowicka and Joe Jacobian).
corr_L_from_cpcs <- function(z, K) {
  entries <- list(1)
  rows <- 1L
  cols <- 1L
  # acc[i]: remaining squared length for row i
  acc <- vector("list", K)
  for (i in 2:K) acc[[i]] <- 1
  pos <- 0L
  for (j in 1:(K - 1L)) {
    if (j > 1L) {
      entries[[length(entries) + 1L]] <- sqrt(acc[[j]])
      rows <- c(rows, j)
      cols <- c(cols, j)
    }
    for (i in (j + 1L):K) {
      pos <- pos + 1L
      cpc <- ad_select(z, pos)
      lij <- if (j == 1L) cpc else cpc * sqrt(acc[[i]])
      entries[[length(entries) + 1L]] <- lij
      rows <- c(rows, i)
      cols <- c(cols, j)
      acc[[i]] <- acc[[i]] * (1 - cpc * cpc)
    }
  }
  entries[[length(entries) + 1L]] <- sqrt(acc[[K]])
  rows <- c(rows, K)
  cols <- c(cols, K)
  assemble_matrix(entries, rows, cols, K)
}

tf_corr_matrix <- list(
  free_size = function(shape) shape[1L] * (shape[1L] - 1L) / 2L,
  constrain = function(u, shape, a) {
    K <- shape[1L]
    if (K == 0L) return(list(x = matrix(0, 0, 0), lj = 0))
    if (K == 1L) return(list(x = matrix(1, 1, 1), lj = 0))
    z <- tanh(u)
    lj <- ad_sum(ad_log1m(z * z))
    w <- numeric(length(val(z)))
    pos <- 0L
    for (k in 1:(K - 1L)) for (i in (k + 1L):K) {
      pos <- pos + 1L
      w[pos] <- K - k - 1L
    }
    lj <- lj + 0.5 * ad_sum(w * ad_log1m(z * z))
    L <- corr_L_from_cpcs(z, K)
    list(x = ad_matmul(L, ad_transpose(L)), lj = lj)
  },
  unconstrain = function(x, shape, a) {
    S <- as.matrix(x)
    K <- nrow(S)
    if (K <= 1L) return(numeric(0))
    L <- t(chol(S))
    out <- numeric(K * (K - 1L) / 2L)
    pos <- 0L
    for (j in 1:(K - 1L)) for (i in (j + 1L):K) {
      pos <- pos + 1L
      rem <- if (j == 1L) 1 else 1 - sum(L[i, seq_len(j - 1L)]^2)
      out[pos] <- atanh(L[i, j] / sqrt(rem))
    }
    out
  },
  check = function(x, shape, a) {
    S <- as.matrix(x)
    if (nrow(S) != ncol(S)) return("is not a square matrix")
    if (max(abs(S - t(S))) > 1e-8) return("is not symmetric")
    if (any(abs(diag(S) - 1) > 1e-8)) return("is not a valid correlation matrix: diagonal elements must be 1")
    if (!is_pos_def(S)) return("is not positive definite")
    NULL
  }
)

tf_cov_matrix <- list(
  free_size = function(shape) shape[1L] + shape[1L] * (shape[1L] - 1L) / 2L,
  constrain = function(u, shape, a) {
    K <- shape[1L]
    if (K == 0L) return(list(x = matrix(0, 0, 0), lj = 0))
    rows <- integer(0)
    cols <- integer(0)
    diag_pos <- integer(0)
    pos <- 0L
    for (m in 1:K) {
      if (m > 1L) for (j in 1:(m - 1L)) {
        pos <- pos + 1L
        rows <- c(rows, m)
        cols <- c(cols, j)
      }
      pos <- pos + 1L
      rows <- c(rows, m)
      cols <- c(cols, m)
      diag_pos <- c(diag_pos, pos)
    }
    is_diag <- seq_len(pos) %in% diag_pos
    vals <- ad_ifelse(is_diag, exp(u), u)
    L <- ad_scatter(matrix(0, K, K), (cols - 1L) * K + rows, vals)
    lj <- K * log(2) + ad_sum((K - seq_len(K) + 2) * ad_gather(u, diag_pos))
    list(x = ad_matmul(L, ad_transpose(L)), lj = lj)
  },
  unconstrain = function(x, shape, a) {
    S <- as.matrix(x)
    K <- nrow(S)
    L <- t(chol(S))
    out <- numeric(0)
    for (m in 1:K) {
      if (m > 1L) out <- c(out, L[m, seq_len(m - 1L)])
      out <- c(out, log(L[m, m]))
    }
    out
  },
  check = function(x, shape, a) {
    S <- as.matrix(x)
    if (nrow(S) != ncol(S)) return("is not a square matrix")
    if (max(abs(S - t(S))) > 1e-8 * max(1, max(abs(S)))) return("is not symmetric")
    if (!is_pos_def(S)) return("is not positive definite")
    NULL
  }
)

tf_cholesky_cov <- list(
  free_size = function(shape) {
    M <- shape[1L]
    N <- shape[2L]
    N * (N + 1L) / 2L + (M - N) * N
  },
  constrain = function(u, shape, a) {
    M <- shape[1L]
    N <- shape[2L]
    rows <- integer(0)
    cols <- integer(0)
    diag_pos <- integer(0)
    pos <- 0L
    for (m in seq_len(N)) {
      if (m > 1L) for (j in 1:(m - 1L)) {
        pos <- pos + 1L
        rows <- c(rows, m)
        cols <- c(cols, j)
      }
      pos <- pos + 1L
      rows <- c(rows, m)
      cols <- c(cols, m)
      diag_pos <- c(diag_pos, pos)
    }
    if (M > N) for (m in (N + 1L):M) for (j in seq_len(N)) {
      pos <- pos + 1L
      rows <- c(rows, m)
      cols <- c(cols, j)
    }
    if (pos == 0L) return(list(x = matrix(0, M, N), lj = 0))
    is_diag <- seq_len(pos) %in% diag_pos
    vals <- ad_ifelse(is_diag, exp(u), u)
    L <- ad_scatter(matrix(0, M, N), (cols - 1L) * M + rows, vals)
    list(x = L, lj = ad_sum(ad_gather(u, diag_pos)))
  },
  unconstrain = function(x, shape, a) {
    L <- as.matrix(x)
    M <- nrow(L)
    N <- ncol(L)
    out <- numeric(0)
    for (m in seq_len(N)) {
      if (m > 1L) out <- c(out, L[m, seq_len(m - 1L)])
      out <- c(out, log(L[m, m]))
    }
    if (M > N) for (m in (N + 1L):M) out <- c(out, L[m, ])
    out
  },
  check = function(x, shape, a) {
    L <- as.matrix(x)
    if (nrow(L) < ncol(L)) return("has fewer rows than columns")
    if (any(abs(L[upper.tri(L)]) > 0)) return("is not lower triangular")
    if (any(diag(L)[seq_len(ncol(L))] <= 0)) return("is not a valid Cholesky factor: diagonal elements must be positive")
    NULL
  }
)

tf_column_stochastic <- list(
  free_size = function(shape) (shape[1L] - 1L) * shape[2L],
  constrain = function(u, shape, a) {
    N <- shape[1L]
    M <- shape[2L]
    cols <- vector("list", M)
    lj <- 0
    for (j in seq_len(M)) {
      s <- simplex_constrain(ad_gather(u, (j - 1L) * (N - 1L) + seq_len(N - 1L)), N)
      cols[[j]] <- s$x
      lj <- lj + s$lj
    }
    list(x = ad_reshape(ad_c(cols), c(N, M)), lj = lj)
  },
  unconstrain = function(x, shape, a) {
    X <- as.matrix(x)
    unlist(lapply(seq_len(ncol(X)), function(j) simplex_free(X[, j])))
  },
  check = function(x, shape, a) {
    X <- as.matrix(x)
    for (j in seq_len(ncol(X))) {
      m <- simplex_check(X[, j])
      if (!is.null(m)) return(paste("has a column that", m))
    }
    NULL
  }
)

tf_row_stochastic <- list(
  free_size = function(shape) shape[1L] * (shape[2L] - 1L),
  constrain = function(u, shape, a) {
    N <- shape[1L]
    M <- shape[2L]
    # free values form an N x (M - 1) matrix in column major order
    rows <- vector("list", N)
    lj <- 0
    for (i in seq_len(N)) {
      s <- simplex_constrain(ad_gather(u, i + N * (seq_len(M - 1L) - 1L)), M)
      rows[[i]] <- s$x
      lj <- lj + s$lj
    }
    X <- ad_transpose(ad_reshape(ad_c(rows), c(M, N)))
    list(x = X, lj = lj)
  },
  unconstrain = function(x, shape, a) {
    X <- as.matrix(x)
    Y <- t(vapply(seq_len(nrow(X)), function(i) simplex_free(X[i, ]), numeric(ncol(X) - 1L)))
    if (ncol(X) - 1L == 1L) Y <- matrix(Y, ncol = 1L)
    as.vector(Y)
  },
  check = function(x, shape, a) {
    X <- as.matrix(x)
    for (i in seq_len(nrow(X))) {
      m <- simplex_check(X[i, ])
      if (!is.null(m)) return(paste("has a row that", m))
    }
    NULL
  }
)

transform_for <- function(kind) {
  switch(kind,
    none = tf_identity,
    bounds = tf_bounds,
    ordered = tf_ordered,
    positive_ordered = tf_positive_ordered,
    simplex = tf_simplex,
    unit_vector = tf_unit_vector,
    sum_to_zero_vector = tf_sum_to_zero,
    sum_to_zero_matrix = tf_sum_to_zero_matrix,
    cholesky_factor_corr = tf_cholesky_corr,
    cholesky_factor_cov = tf_cholesky_cov,
    corr_matrix = tf_corr_matrix,
    cov_matrix = tf_cov_matrix,
    column_stochastic_matrix = tf_column_stochastic,
    row_stochastic_matrix = tf_row_stochastic,
    stop("internal error: unknown transform ", kind, call. = FALSE))
}

is_pos_def <- function(S) {
  if (any(!is.finite(S))) return(FALSE)
  r <- tryCatch(chol(S), error = function(e) NULL)
  !is.null(r)
}

fmt_num <- function(x) format(signif(x, 6), scientific = NA)

# ---------------------------------------------------------------------------
# Array layout helpers.
#
# A declared variable with array dimensions `adims` and element shape `eshape`
# is stored as an R array with dim c(adims, eshape) (column major). Stan reads
# unconstrained values element by element, array elements in row major order,
# so the concatenated element layout is an array with dim c(eshape, rev(adims)).

layout_perm <- function(adims, eshape) {
  k <- length(adims)
  m <- length(eshape)
  total <- prod(c(adims, eshape))
  if (k == 0L || total == 0L) return(NULL)
  idx <- array(seq_len(total), dim = c(eshape, rev(adims), 1L)[seq_len(m + k)])
  if (m + k == 1L) return(seq_len(total))
  as.vector(aperm(idx, c((m + k):(m + 1L), seq_len(m))))
}

# from storage layout to concatenated element layout
to_element_order <- function(x, adims, eshape) {
  p <- layout_perm(adims, eshape)
  x <- as.vector(x)
  if (is.null(p)) return(x)
  out <- numeric(length(x))
  out[p] <- x
  out
}

full_shape <- function(adims, eshape) {
  d <- c(adims, eshape)
  if (length(d) <= 1L) NULL else d
}

# Constrain one declared variable.
#   u      differentiable unconstrained vector for the whole variable
#   decl   list(kind, adims, eshape)
#   args   evaluated constraint arguments (lower, upper, offset, multiplier)
constrain_variable <- function(u, decl, args) {
  kind <- decl$kind
  adims <- decl$adims
  eshape <- decl$eshape
  # fast path: an unconstrained variable whose storage order is the
  # element order (a scalar, a plain vector, matrix or array of reals)
  if (kind == "none" && !length(args) && (!length(adims) || !length(eshape))) {
    fs <- full_shape(adims, eshape)
    return(list(x = if (is.null(fs) && is.null(dim(val(u)))) u else ad_reshape(u, fs), lj = 0))
  }
  if (identical(kind, "tuple")) return(tuple_constrain(u, decl, args))
  tf <- transform_for(kind)
  n_elem <- if (length(adims)) prod(adims) else 1L
  fs <- tf$free_size(eshape)
  if (decl$kind %in% c("none", "bounds")) {
    # elementwise transform on the whole variable in element order
    a <- lapply(args, function(b) bound_to_element_order(b, adims, eshape))
    r <- tf$constrain(u, eshape, a)
    x <- r$x
    p <- layout_perm(adims, eshape)
    if (!is.null(p)) x <- ad_gather(x, p, full_shape(adims, eshape))
    else x <- ad_reshape(x, full_shape(adims, eshape))
    return(list(x = x, lj = r$lj))
  }
  if (n_elem == 1L && length(adims) == 0L) {
    r <- tf$constrain(u, eshape, args)
    return(list(x = ad_reshape(r$x, full_shape(adims, eshape)), lj = r$lj))
  }
  parts <- vector("list", n_elem)
  lj <- 0
  for (e in seq_len(n_elem)) {
    ue <- if (fs > 0L) ad_gather(u, (e - 1L) * fs + seq_len(fs)) else numeric(0)
    r <- tf$constrain(ue, eshape, args)
    parts[[e]] <- r$x
    lj <- lj + r$lj
  }
  x <- ad_c(parts)
  p <- layout_perm(adims, eshape)
  x <- ad_gather(x, p, full_shape(adims, eshape))
  list(x = x, lj = lj)
}

bound_to_element_order <- function(b, adims, eshape) {
  if (is.null(b)) return(NULL)
  vb <- val(b)
  if (length(vb) == 1L) return(b)
  total <- prod(c(adims, eshape))
  ne <- prod(eshape)
  if (length(vb) == total) {
    p <- layout_perm(adims, eshape)
    if (is.null(p)) return(ad_reshape(b, NULL))
    inv <- order(p)
    return(ad_gather(b, inv))
  }
  if (length(vb) == ne && length(adims)) {
    return(ad_gather(b, rep(seq_len(ne), prod(adims))))
  }
  ad_reshape(b, NULL)
}

unconstrain_variable <- function(x, decl, args) {
  if (identical(decl$kind, "tuple")) {
    msg <- tuple_check(x, decl, args)
    if (!is.null(msg)) stop_data(sprintf("parameter '%s' %s", decl$name %||% "", msg))
    return(tuple_unconstrain(x, decl, args))
  }
  msg <- check_variable(x, decl, args)
  if (!is.null(msg)) stop_data(sprintf("parameter '%s' %s", decl$name %||% "", msg))
  tf <- transform_for(decl$kind)
  adims <- decl$adims
  eshape <- decl$eshape
  if (decl$kind %in% c("none", "bounds")) {
    a <- lapply(args, function(b) if (is.null(b)) NULL else val(bound_to_element_order(b, adims, eshape)))
    xe <- to_element_order(x, adims, eshape)
    return(tf$unconstrain(xe, eshape, a))
  }
  n_elem <- if (length(adims)) prod(adims) else 1L
  if (length(adims) == 0L) return(tf$unconstrain(x, eshape, args))
  xe <- to_element_order(x, adims, eshape)
  ne <- prod(eshape)
  unlist(lapply(seq_len(n_elem), function(e) {
    xi <- xe[(e - 1L) * ne + seq_len(ne)]
    if (length(eshape) == 2L) dim(xi) <- eshape
    tf$unconstrain(xi, eshape, args)
  }))
}

check_variable <- function(x, decl, args) {
  if (identical(decl$kind, "tuple")) return(tuple_check(x, decl, args))
  tf <- transform_for(decl$kind)
  adims <- decl$adims
  eshape <- decl$eshape
  if (decl$kind %in% c("none", "bounds")) {
    a <- lapply(args, function(b) if (is.null(b)) NULL else val(bound_to_element_order(b, adims, eshape)))
    return(tf$check(to_element_order(x, adims, eshape), eshape, a))
  }
  if (length(adims) == 0L) return(tf$check(x, eshape, args))
  xe <- to_element_order(x, adims, eshape)
  ne <- prod(eshape)
  for (e in seq_len(prod(adims))) {
    xi <- xe[(e - 1L) * ne + seq_len(ne)]
    if (length(eshape) == 2L) dim(xi) <- eshape
    m <- tf$check(xi, eshape, args)
    if (!is.null(m)) return(paste0("element ", e, " ", m))
  }
  NULL
}

free_size_variable <- function(decl) {
  if (identical(decl$kind, "tuple")) return(decl$free_size)
  tf <- transform_for(decl$kind)
  n_elem <- if (length(decl$adims)) prod(decl$adims) else 1L
  as.integer(n_elem * tf$free_size(decl$eshape))
}
