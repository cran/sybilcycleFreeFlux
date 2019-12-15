\name{llFBA}
\alias{llFBA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
looplessFBA
}
\description{
implements loopless FBA as described in Schellenberger et al 2011 algorithm.
}
\usage{
llFBA(model, lpdir = SYBIL_SETTINGS("OPT_DIRECTION"), solver = SYBIL_SETTINGS("SOLVER"),
   method = SYBIL_SETTINGS("METHOD"), solverParm = data.frame(CPX_PARAM_EPRHS = 1e-07),
  verboseMode = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
 An object of class \code{\link{modelorg}}.
}
   \item{lpdir}{
  Character value, direction of optimisation.  Can be set to \code{"min"} or
    \code{"max"}.\cr
    Default: \code{SYBIL_SETTINGS("OPT_DIRECTION")}.
}
  \item{solver}{
Single character value.  The solver to use.  See
    \code{\link{SYBIL_SETTINGS}} for possible values.\cr
    Default: \code{SYBIL_SETTINGS("SOLVER")}.
}
   \item{method}{
 Single character value.  The optimization algorithm to use.  Possible values
    depend on the setting in \code{solver}.  See \code{\link{SYBIL_SETTINGS}} for
    possible values.\cr
    Default: \code{LP_METHOD(SYBIL_SETTINGS)}.
}
  \item{solverParm}{
additional parameters to the solver
}
  \item{verboseMode}{
level of output messages
}
}

\value{
status returned from the solver, objective value, and the fluxes returned from solver.
}
\references{
Schellenberger J, Lewis NE, Palsson BO (2011) Elimination of thermodynamically infeasible
 loops in steadystate metabolic models. Biophysical journal 100: 544 53
}

\author{Abdelmoneim Amer Desouki}  


\seealso{
\code{\link{modelorg}}
\code{\link{cfFBA}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
library(sybilcycleFreeFlux)
data(Ec_core)
model=Ec_core
llFBA(model,solver="glpkAPI",verbose=2)

## The function is currently defined as
"llFBA"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
