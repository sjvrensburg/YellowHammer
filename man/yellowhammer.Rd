% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/YellowHammer.R
\name{yellowhammer}
\alias{yellowhammer}
\title{Nonlinear optimization with constraints}
\usage{
yellowhammer(
  par,
  fn,
  gr,
  hin,
  hin.jac,
  control.outer = list(),
  control.optim = list(),
  ...
)
}
\arguments{
\item{par}{Starting vector of parameter values.  Any initial vector, even those violating inequality constraints, may be specified.  This is in contrast to \code{constrOptim.nl} which requires  "feasible" initial values with respect to inequality constraints}

\item{fn}{Nonlinear objective function that is to be optimized. A scalar function that takes a real vector as argument and returns a scalar that is the value of the function at that point (see details).}

\item{gr}{The gradient of the objective function \code{fn} evaluated at the argument.  This is a vector-function that takes a real vector as argument and returns a real vector of the same length. It defaults to "NULL", which means that gradient is evaluated numerically.  Computations are dramatically faster in high-dimensional problems when the exact gradient is provided.}

\item{hin}{a vector function specifying inequality constraints such that hin\link{j} > 0 for all j}

\item{hin.jac}{Jacobian of \code{hin}.  If unspecified, it will be computed using finite-difference, but computations will be faster if specified.}

\item{control.outer}{A list of control parameters to be used by the outer loop in \code{constrOptim.nl}.  See \emph{Details} for more information.}

\item{control.optim}{A list of control parameters to be used by the unconstrained optimization algorithm in the inner loop. Identical to that used in \code{marqLevAlg}. See \code{marqLevAlg}'s documentation for details.}

\item{...}{Additional arguments passed to \code{fn}, \code{gr}, \code{hin}, \code{heq}.  All of them must accept any specified arguments, either explicitly or by having a \dots argument, but they do not need to use them all.}
}
\value{
A list with the following components:
\item{par}{Parameters that optimize the nonlinear objective function, satisfying constraints, if convergence is successful.}
\item{value}{The value of the objective function at termination.}
\item{counts}{A vector of length 2 denoting the number of times the objective \code{fn} and the \code{gr} were evaluated, respectively.}
\item{convergence}{An integer code indicating type of convergence.  \code{0} indicates successful convergence. Positive integer codes indicate failure to converge.}
\item{outer.iterations}{Number of outer iterations}
\item{lambda}{Values of the Lagrangian parameter.  This is a vector of same length as the total number of inequalities and equalities.  It must be zero for inactive inequalities; non-negative for active inequalities; and can have any sign for equalities.}
\item{sigma}{Value of augmented penalty parameter for the quadratic term}
\item{gradient}{Gradient of the augmented Lagrangian function at convergence. It should be small.}
\item{hessian}{Hessian of the augmented Lagrangian function at convergence. It should be positive (negative) definite for minimization (maximization)}
\item{ineq}{Values of inequlaity constraints at convergence. All of them must be non-negative}
\item{equal}{Values of equlaity constraints at convergence. All of them must be close to zero.}
\item{kkt1}{A logical variable indicating whether or not the first-order KKT conditions were satisfied.}
\item{kkt2}{A logical variable indicating whether or not the second-order KKT conditions were satisfied.}
}
\description{
Augmented Lagrangian Minimization Algorithm for optimizing smooth nonlinear objective functions with constraints. Linear or nonlinear inequality constraints are allowed.
}
\details{
Argument \code{control.outer} is a list specifying any changes to default values of algorithm control parameters for the outer loop.  Note that the names of these must be specified completely.  Partial matching will not work. The list items are as follows:
\code{lam0}: Initial value for the Lagrangian parameter.
\code{sig0}: A scaling parameter for penalty term that is augmented to the Lagrangian.
\code{eps}: Tolerance for convergence of outer iterations of the barrier and/or augmented lagrangian algorithm
\code{itmax}: Maximum number of outer iterations.
\code{ilack.max}: Maximum number of outer iterations where no change in parameters is tolerated.
\code{trace}: A logical variable indicating whether information on outer iterations should be printed out.  If TRUE, at each outer iteration information is displayed on: (i) how well the inequality and equalities are satisfied, (ii) current parameter values, and (iii) current objective function value.
\code{NMinit}: A logical variable indicating whether "Nelder-Mead" algorithm should be used for the first outer iteration.
\code{i.scale}: A vector of length equal to number of inequalities that may be used to scale the inequalities or it can be a scalar in which case all the inequalities are scaled by the same value.
\code{kkt2.check}: A logical variable (TRUE/FALSE) indicating whether the second-order KKT condition should be checked.  Deafult is TRUE.  It may be set to FALSE in problems where the Hessian computation can b etime consuming.
}
