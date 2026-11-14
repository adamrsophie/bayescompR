# Higher order functions: the ODE, DAE and algebraic solvers, the 1D
# integrator, map_rect and reduce_sum. Their first argument names a user
# defined function. Each entry gives the types of the fixed arguments after
# the function (with the ones that must be data), the argument types the
# function must accept given the call's argument types, the return types and
# the runtime implementation.

hof_entry <- function(fixed, data, fargs, fret, ret, impl, variadic = TRUE, opts = list()) {
  list(fixed = fixed, data = data, fargs = fargs, fret = fret, ret = ret, impl = impl, variadic = variadic, opts = opts)
}

T_ra <- ty("real", 1L); T_ia <- ty("int", 1L); T_va <- ty("vector", 1L)
T_ra2 <- ty("real", 2L); T_ia2 <- ty("int", 2L)

ode_fargs <- function(fixed, va) c(list(T_real, T_vec), va)
hof_specs <- local({
  sp <- list()
  for (m in c("rk45", "ckrk", "adams", "bdf")) local({
    meth <- m
    sp[[paste0("ode_", meth)]] <<- hof_entry(list(T_vec, T_real, T_ra), c(FALSE, FALSE, FALSE), ode_fargs, T_vec, T_va,
                                             function(...) rt_ode(...), opts = list(method = meth, tol = FALSE))
    sp[[paste0("ode_", meth, "_tol")]] <<- hof_entry(list(T_vec, T_real, T_ra, T_real, T_real, T_int), c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
                                                     ode_fargs, T_vec, T_va, function(...) rt_ode(...), opts = list(method = meth, tol = TRUE))
  })
  sp$ode_adjoint_tol_ctl <- hof_entry(list(T_vec, T_real, T_ra, T_real, T_vec, T_real, T_vec, T_real, T_real, T_int, T_int, T_int, T_int, T_int),
                                      c(FALSE, FALSE, FALSE, rep(TRUE, 11L)), ode_fargs, T_vec, T_va,
                                      function(...) rt_ode(...), opts = list(method = "adjoint", tol = TRUE))
  # DAE
  dae_f <- function(fixed, va) c(list(T_real, T_vec, T_vec), va)
  sp$dae <- hof_entry(list(T_vec, T_vec, T_real, T_ra), c(FALSE, FALSE, FALSE, FALSE), dae_f, T_vec, T_va,
                      function(...) rt_dae(...), opts = list(tol = FALSE))
  sp$dae_tol <- hof_entry(list(T_vec, T_vec, T_real, T_ra, T_real, T_real, T_int), c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
                          dae_f, T_vec, T_va, function(...) rt_dae(...), opts = list(tol = TRUE))
  # algebraic solvers
  alg_f <- function(fixed, va) c(list(T_vec), va)
  for (m in c("newton", "powell")) local({
    meth <- m
    sp[[paste0("solve_", meth)]] <<- hof_entry(list(T_vec), FALSE, alg_f, T_vec, T_vec, function(...) rt_solve(...),
                                               opts = list(method = meth, tol = FALSE))
    sp[[paste0("solve_", meth, "_tol")]] <<- hof_entry(list(T_vec, T_real, T_real, T_int), c(FALSE, TRUE, TRUE, TRUE), alg_f, T_vec, T_vec,
                                                       function(...) rt_solve(...), opts = list(method = meth, tol = TRUE))
  })
  legacy_alg <- function(fixed, va) list(T_vec, T_vec, T_ra, T_ia)
  for (nm in c("algebra_solver", "algebra_solver_newton")) local({
    meth <- if (nm == "algebra_solver") "powell" else "newton"
    sp[[nm]] <<- hof_entry(list(T_vec, T_vec, T_ra, T_ia), c(FALSE, FALSE, TRUE, TRUE), legacy_alg, T_vec, T_vec,
                           function(...) rt_solve_legacy(...), variadic = FALSE,
                           opts = list(method = meth, alt = list(T_vec, T_vec, T_ra, T_ia, T_real, T_real, T_real)))
  })
  # 1D integration
  sp$integrate_1d <- hof_entry(list(T_real, T_real, T_ra, T_ra, T_ia), c(FALSE, FALSE, FALSE, TRUE, TRUE),
                               function(fixed, va) list(T_real, T_real, T_ra, T_ra, T_ia), T_real, T_real,
                               function(...) rt_integrate_1d(...), variadic = FALSE,
                               opts = list(alt = list(T_real, T_real, T_ra, T_ra, T_ia, T_real)))
  # map_rect
  sp$map_rect <- hof_entry(list(T_vec, T_va, T_ra2, T_ia2), c(FALSE, FALSE, TRUE, TRUE),
                           function(fixed, va) list(T_vec, T_vec, T_ra, T_ia), T_vec, T_vec,
                           function(...) rt_map_rect(...), variadic = FALSE)
  # reduce_sum: the sliced argument may be an array of any type
  for (nm in c("reduce_sum", "reduce_sum_static")) {
    sp[[nm]] <- hof_entry(list(NULL, T_int), c(FALSE, FALSE),
                          function(fixed, va) c(list(fixed[[1L]], T_int, T_int), va), T_real, T_real,
                          function(...) rt_reduce_sum(...))
  }
  sp
})

# does the user function signature sg accept the expected argument types
hof_fun_ok <- function(sg, expected) {
  length(sg$args) == length(expected) && all(mapply(ty_assignable, sg$args, expected))
}

sem_hof <- function(s, e) {
  nm <- e$name
  spec <- hof_specs[[nm]]
  if (length(e$args) < 1L || e$args[[1L]]$kind != "var")
    stop_semantic(sprintf("the first argument of '%s' must be the name of a user defined function", nm), e$loc)
  fname <- e$args[[1L]]$name
  if (grepl("_lup(d|m)f$", fname)) fname <- sub("_lup", "_lp", fname)
  cands <- s$ufuns[[fname]]
  if (is.null(cands)) stop_semantic(sprintf("'%s' is not a user defined function", e$args[[1L]]$name), e$args[[1L]]$loc)
  args <- lapply(e$args[-1L], function(a) sem_expr(s, a))
  types <- lapply(args, `[[`, "type")
  fixed <- spec$fixed
  data_req <- spec$data
  nf <- length(fixed)
  # legacy forms with optional control arguments
  if (!spec$variadic && !is.null(spec$opts$alt) && length(types) == length(spec$opts$alt)) {
    fixed <- spec$opts$alt
    data_req <- c(data_req, rep(TRUE, length(fixed) - nf))
    nf <- length(fixed)
  }
  if (length(types) < nf || (!spec$variadic && length(types) != nf))
    stop_semantic(sprintf("'%s' takes %d arguments after the function, found %d", nm, nf, length(types)), e$loc)
  for (i in seq_len(nf)) {
    want <- fixed[[i]]
    if (is.null(want)) {
      # reduce_sum: an array to slice
      if (types[[i]]$nd < 1L) stop_semantic(sprintf("argument %d of '%s' must be an array", i + 1L, nm), args[[i]]$loc)
      fixed[[i]] <- types[[i]]
      next
    }
    if (!ty_assignable(want, types[[i]]))
      stop_semantic(sprintf("argument %d of '%s' must be %s, found %s", i + 1L, nm, ty_str(want), ty_str(types[[i]])), args[[i]]$loc)
    if (data_req[i] && !isTRUE(args[[i]]$data_only))
      stop_semantic(sprintf("argument %d of '%s' must be data only", i + 1L, nm), args[[i]]$loc)
  }
  va <- types[-seq_len(nf)]
  expected <- spec$fargs(fixed, va)
  sg <- NULL
  for (cand in cands) if (hof_fun_ok(cand, expected)) { sg <- cand; break }
  if (is.null(sg))
    stop_semantic(sprintf("function '%s' passed to '%s' must take (%s)", fname, nm,
                          paste(vapply(expected, ty_str, ""), collapse = ", ")), e$loc)
  if (!ty_assignable(spec$fret, sg$ret))
    stop_semantic(sprintf("function '%s' passed to '%s' must return %s", fname, nm, ty_str(spec$fret)), e$loc)
  # arguments the function declares data must be data
  for (k in seq_along(va)) {
    pos <- length(expected) - length(va) + k
    if (isTRUE(sg$data_only[pos]) && !isTRUE(args[[nf + k]]$data_only))
      stop_semantic(sprintf("argument %d of '%s' must be data only, as '%s' declares it data", nf + k + 1L, nm, fname), args[[nf + k]]$loc)
  }
  e$args <- args
  e$hof <- list(name = nm, fun = sg, fixed = fixed, nf = nf, expected = expected, opts = spec$opts, impl = spec$impl)
  e$type <- spec$ret
  e$data_only <- all(vapply(args, function(a) isTRUE(a$data_only), TRUE))
  e
}

cg_hof <- function(e) {
  h <- e$hof
  codes <- lapply(e$args, cg_expr)
  nf <- h$nf
  fixed <- Map(function(code, a, t) promote_to(code, a$type, t), codes[seq_len(nf)], e$args[seq_len(nf)], h$fixed)
  vexp <- h$fun$args[-seq_len(length(h$fun$args) - (length(e$args) - nf))]
  va <- if (length(e$args) > nf)
    Map(function(code, a, t) promote_to(code, a$type, t), codes[-seq_len(nf)], e$args[-seq_len(nf)], vexp) else list()
  cl(h$impl, as.name(h$fun$mangled), cl_list(as.name("list"), unname(fixed)), cl_list(as.name("list"), unname(va)),
     CTX, c(h$opts, list(name = h$name, nargs = nf)))
}
