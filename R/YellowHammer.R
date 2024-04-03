#' Nonlinear optimization with constraints
#'
#' Augmented Lagrangian Minimization Algorithm for optimizing smooth nonlinear objective functions with constraints. Linear or nonlinear inequality constraints are allowed.
#'
#' @param par Starting vector of parameter values.  Any initial vector, even those violating inequality constraints, may be specified.  This is in contrast to \code{constrOptim.nl} which requires  "feasible" initial values with respect to inequality constraints
#' @param fn Nonlinear objective function that is to be optimized. A scalar function that takes a real vector as argument and returns a scalar that is the value of the function at that point (see details).
#' @param gr The gradient of the objective function \code{fn} evaluated at the argument.  This is a vector-function that takes a real vector as argument and returns a real vector of the same length. It defaults to "NULL", which means that gradient is evaluated numerically.  Computations are dramatically faster in high-dimensional problems when the exact gradient is provided.
#' @param hin a vector function specifying inequality constraints such that hin[j] > 0 for all j
#' @param hin.jac Jacobian of \code{hin}.  If unspecified, it will be computed using finite-difference, but computations will be faster if specified.
#' @param control.outer A list of control parameters to be used by the outer loop in \code{constrOptim.nl}.  See *Details* for more information.
#' @param control.optim A list of control parameters to be used by the unconstrained optimization algorithm in the inner loop. Identical to that used in \code{marqLevAlg}. See \code{marqLevAlg}'s documentation for details.
#' @param ... Additional arguments passed to \code{fn}, \code{gr}, \code{hin}, \code{heq}.  All of them must accept any specified arguments, either explicitly or by having a \dots argument, but they do not need to use them all.
#'
#' @details
#' Argument \code{control.outer} is a list specifying any changes to default values of algorithm control parameters for the outer loop.  Note that the names of these must be specified completely.  Partial matching will not work. The list items are as follows:
#' \code{lam0}: Initial value for the Lagrangian parameter.
#' \code{sig0}: A scaling parameter for penalty term that is augmented to the Lagrangian.
#' \code{eps}: Tolerance for convergence of outer iterations of the barrier and/or augmented lagrangian algorithm
#' \code{itmax}: Maximum number of outer iterations.
#' \code{ilack.max}: Maximum number of outer iterations where no change in parameters is tolerated.
#' \code{trace}: A logical variable indicating whether information on outer iterations should be printed out.  If TRUE, at each outer iteration information is displayed on: (i) how well the inequality and equalities are satisfied, (ii) current parameter values, and (iii) current objective function value.
#' \code{NMinit}: A logical variable indicating whether "Nelder-Mead" algorithm should be used for the first outer iteration.
#' \code{i.scale}: A vector of length equal to number of inequalities that may be used to scale the inequalities or it can be a scalar in which case all the inequalities are scaled by the same value.
#' \code{kkt2.check}: A logical variable (TRUE/FALSE) indicating whether the second-order KKT condition should be checked.  Deafult is TRUE.  It may be set to FALSE in problems where the Hessian computation can b etime consuming.
#'
#' @return A list with the following components:
#'  \item{par}{Parameters that optimize the nonlinear objective function, satisfying constraints, if convergence is successful.}
#'  \item{value}{The value of the objective function at termination.}
#'  \item{counts}{A vector of length 2 denoting the number of times the objective \code{fn} and the \code{gr} were evaluated, respectively.}
#'  \item{convergence}{An integer code indicating type of convergence.  \code{0} indicates successful convergence. Positive integer codes indicate failure to converge.}
#'  \item{outer.iterations}{Number of outer iterations}
#'  \item{lambda}{Values of the Lagrangian parameter.  This is a vector of same length as the total number of inequalities and equalities.  It must be zero for inactive inequalities; non-negative for active inequalities; and can have any sign for equalities.}
#'  \item{sigma}{Value of augmented penalty parameter for the quadratic term}
#'  \item{gradient}{Gradient of the augmented Lagrangian function at convergence. It should be small.}
#'  \item{hessian}{Hessian of the augmented Lagrangian function at convergence. It should be positive (negative) definite for minimization (maximization)}
#'  \item{ineq}{Values of inequlaity constraints at convergence. All of them must be non-negative}
#'  \item{equal}{Values of equlaity constraints at convergence. All of them must be close to zero.}
#'  \item{kkt1}{A logical variable indicating whether or not the first-order KKT conditions were satisfied.}
#'  \item{kkt2}{A logical variable indicating whether or not the second-order KKT conditions were satisfied.}
#' @export
yellowhammer <- function(
    par, fn, gr, hin, hin.jac, control.outer = list(), control.optim = list(),
    ...) {
  # Sanity Checks
  if (missing(hin)) stop("This is an unconstrained optimization problem - you should use `optim' \n")

  # Set/Modify Options
  control.outer.default <- list(
    lam0 = 10, sig0 = 100, eps = 1e-07,
    itmax = 50, trace = TRUE, NMinit = FALSE, ilack.max = 6,
    i.scale = 1, kkt2.check = TRUE
  )
  control.optim.default <- list(
    maxiter = 500, epsa = 1e-04, epsb = 1e-04, epsd = 1e-04,
    partialH = NULL, digits = 8, print.info = FALSE, blinding = TRUE,
    multipleTry = 25, nproc = 1, clustertype = NULL, file = "",
    .packages = NULL, minimize = TRUE)
  control.outer <- modifyList(control.outer.default, control.outer)
  control.inner <- modifyList(control.optim.default, control.optim)
  i.scale <- control.outer$i.scale

  # Gradients, scaling and Jacobian matrices
  if (missing(gr)) {
    gr <- function(par, ...) {
      numDeriv::grad(func = fn, x = par, method = "simple", ...)
    }
  }
  hin.scaled <- function(par, ...) {
    hin(par, ...) / i.scale
  }
  hin.jac.scaled <- if (missing(hin.jac)) {
    function(par, lower, upper, ...) {
      numDeriv::jacobian(func = hin.scaled, x = par, method = "simple", ...)
    }
  } else {
    function(par, lower, upper, ...) {
      hin.jac(par, ...) / i.scale
    }
  }

  ans <- .yellowhammer(par, fn, gr, hin.scaled, hin.jac.scaled,
    control.outer = control.outer, control.optim = control.inner, ...
  )

  ans$ineq <- ans$ineq * i.scale

  return(ans)
}

##################################################################
.yellowhammer <- function(par, fn, gr = NULL,
                    hin = NULL, hin.jac = NULL,
                    control.outer = list(), control.optim = list(),
                    ...) {
  sig <- control.outer$sig0
  lam0 <- control.outer$lam0
  trace <- control.outer$trace
  eps <- control.outer$eps
  itmax <- control.outer$itmax
  ilack.max <- control.outer$ilack.max
  NMinit <- control.outer$NMinit
  kkt2.check <- control.outer$kkt2.check
  pfact <- if (!is.null(control.optim$fnscale) && control.optim$fnscale <
    0) {
    -1
  } else {
    1
  }

  fun <- function(par, ...) {
    h0 <- hin(par, ...)
    d0 <- h0
    inactive <- (1:length(h0))[(h0 > lam[1:length(h0)] / sig)]
    d0[inactive] <- lam[inactive] / sig
    fn(par, ...) - pfact * sum(lam * d0) +
      pfact * sig / 2 * sum(d0 * d0)
  }

  gradient <- function(par, ...) {
    h0 <- hin(par, ...)
    d0 <- h0
    active <- (1:length(h0))[(h0 <= lam[1:length(h0)] / sig)]
    ij <- hin.jac(par, ...)[active, , drop = FALSE]
    gr(par, ...) - pfact * colSums(lam[active] *
      ij) + pfact * sig * drop(crossprod(ij, d0[active]))
  }

  h0 <- hin(par, ...)
  d0 <- h0
  lam <- rep(lam0, length(d0))
  inactive <- (1:length(h0))[(h0 > lam[1:length(h0)] / sig)]
  d0[inactive] <- lam[inactive] / sig
  dmax <- max(abs(d0))

  obj <- fn(par, ...)
  r <- obj
  feval <- 0
  geval <- 0
  ilack <- 0
  Kprev <- dmax
  sig0 <- sig / Kprev
  if (is.infinite(sig0)) {
    sig0 <- 1
  }
  sig <- sig0

  K <- Inf
  if (trace) cat("Min(hin): ", min(h0), "\n")
  for (i in 1:itmax) {
    if (trace) {
      cat("Outer iteration: ", i, "\n")
      cat("Min(hin): ", min(h0), "\n")
      cat("par: ", signif(par, 6), "\n")
      cat("fval =  ", signif(obj, 4), "\n \n")
    }
    par.old <- par
    obj.old <- obj
    r.old <- r
    if (NMinit & i == 1) {
      a <- optim(par = par, fn = fun, method = "Nelder-Mead", ...)
    } else {
      a <- marqLevAlg::marqLevAlg(
        b = par, fn = fun, gr = gradient,
        maxiter = control.optim$maxiter,
        epsa = control.optim$epsa,
        epsb = control.optim$epsb,
        epsd = control.optim$epsd,
        partialH = control.optim$partialH,
        digits = control.optim$digits,
        print.info = TRUE, #control.optim$print.info,
        blinding = TRUE, #control.optim$blinding,
        multipleTry = control.optim$multipleTry,
        nproc = control.optim$nproc,
        clustertype = control.optim$clustertype,
        file = control.optim$file,
        .packages = control.optim$.packages,
        minimize = control.optim$minimize, ...)
      r <- a$value
      feval <- feval + a$ni
      geval <- geval + a$ni
    }
    par <- a$b
    h0 <- hin(par, ...)
    d0 <- h0
    inactive <- (1:length(h0))[(h0 > lam[1:length(h0)] / sig)]
    d0[inactive] <- lam[inactive] / sig
    K <- max(abs(d0))
    if (K <= Kprev / 4) {
      lam <- lam - d0 * sig
      Kprev <- K
    } else {
      sig <- 10 * sig
    }
    obj <- fn(par, ...)

    pconv <- max(abs(par - par.old))
    if (pconv < eps) {
      ilack <- ilack + 1
    } else {
      ilack <- 0
    }

    if ((is.finite(r) && is.finite(r.old) && abs(r - r.old) <
      eps && K < eps) | ilack >= ilack.max) {
      break
    }
  }

  if (i == itmax) {
    a$convergence <- 7
    a$message <- "ALABaMA ran out of iterations and did not converge"
  } else if (K > eps) {
    a$convergence <- 9
    a$message <- "Convergence due to lack of progress in parameter updates"
  }
  a$outer.iterations <- i
  a$lambda <- lam
  a$sigma <- sig
  a$value <- fn(a$par, ...)
  a$gradient <- gradient(a$par, ...)
  a$ineq <- hin(a$par, ...)
  a$equal <- NA
  a$counts <- c(feval, geval)
  a$kkt1 <- max(abs(a$gradient)) <= 0.01 * (1 + abs(a$value))
  a$kkt2 <- NA
  if (kkt2.check) {
    a$hessian <- jacobian(x = a$par, func = gradient, ...)
    evs <- eigen(a$hessian, symmetric = TRUE, only.values = TRUE)$value
    a$kkt2 <- if (any(abs(Im(evs)) > 1.e-14)) FALSE else if (all(Re(evs) * control.optim$fnscale > 0)) TRUE else FALSE
  }
  a
}
