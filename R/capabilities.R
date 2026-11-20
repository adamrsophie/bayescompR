#' Algorithms behind Stan's solver and integrator interfaces
#'
#' Lists, for each higher order function the language accepts, the
#' algorithm that runs, the control arguments that change it, how its
#' derivatives are computed and whether it runs in parallel. Several Stan
#' names share one algorithm here, so selecting a name does not always give
#' the behavior or scaling of the algorithm Stan uses for it.
#'
#' @return A data frame with columns `functions`, `algorithm`, `controls`
#'   (the controls that are honored), `ignored` (controls that are checked
#'   but have no effect), `derivatives` and `parallel`.
#' @examples
#' bayescomp_capabilities()[, c("functions", "algorithm")]
#' @export
bayescomp_capabilities <- function() {
  row <- function(functions, algorithm, controls, ignored, derivatives, parallel = "no")
    data.frame(functions = functions, algorithm = algorithm, controls = controls, ignored = ignored,
               derivatives = derivatives, parallel = parallel)
  fwd <- "forward sensitivities integrated with the states"
  rbind(
    row("ode_rk45, ode_ckrk, ode_adams (and _tol), integrate_ode_rk45, integrate_ode_adams, integrate_ode",
        "Dormand-Prince 5(4), adaptive step", "relative_tolerance, absolute_tolerance, max_num_steps", "", fwd),
    row("ode_bdf (and _tol), integrate_ode_bdf",
        "Rosenbrock RODAS3 (stiffly accurate, order 3), adaptive step", "relative_tolerance, absolute_tolerance, max_num_steps", "",
        paste(fwd, "(stages solved with one state block factorization)")),
    row("ode_adjoint_tol_ctl", "Dormand-Prince 5(4), adaptive step", "relative_tolerance_forward, absolute_tolerance_forward, max_num_steps",
        "backward and quadrature tolerances, num_steps_between_checkpoints, interpolation_polynomial, solver_forward, solver_backward",
        paste(fwd, "(no adjoint pass)")),
    row("dae, dae_tol", "Radau IIA, 3 stages (order 5), index 1", "relative_tolerance, absolute_tolerance, max_num_steps", "", fwd),
    row("solve_newton, solve_powell (and _tol), algebra_solver, algebra_solver_newton",
        "Newton iteration with backtracking (one solver for all names)", "function_tolerance, max_num_steps",
        "scaling_step_size, relative_tolerance (checked only)", "implicit function theorem"),
    row("integrate_1d, integrate_1d_double_exponential, integrate_1d_gauss_kronrod (and _tol)",
        "adaptive Gauss-Kronrod quadrature (stats::integrate)", "relative_tolerance, absolute_tolerance",
        "max_refinements", "quadrature of the integrand gradient, one adaptive quadrature per component"),
    row("map_rect", "sequential loop over jobs", "", "", "reverse mode through each job"),
    row("reduce_sum, reduce_sum_static", "one call over the whole slice", "", "grainsize", "reverse mode"))
}
