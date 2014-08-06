\name{cfFBA}
\alias{cfFBA}
\alias{lrFBA}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cycle free flux
}
\description{
finds a cycle free flux distribution given a flux distribution and a network model.
}
\usage{
cfFBA(model, wtflux, objVal = NA, fixExchRxn = TRUE, excReactPos = NULL, 
lpdir = SYBIL_SETTINGS("OPT_DIRECTION"), solver = SYBIL_SETTINGS("SOLVER"),
 method = SYBIL_SETTINGS("METHOD"), verboseMode = 2, retOptSol=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
 An object of class \code{\link{modelorg}}.
}
  \item{wtflux}{
initial flux distribution that may contain loops.
}
  \item{objVal}{
value of the objective function.
}
  \item{fixExchRxn}{
a logical value default is true, which indicates that the exchange reactions should be 
fixed or not.
}
  \item{excReactPos}{
list of the exchange reactions that can be sent to avoid recalculation
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
 
  \item{verboseMode}{
level of displaying messages.
}
  \item{retOptSol}{
If set to \code{TRUE}, the function returns an object of class \code{\linkS4class{optsol_optimizeProb}}, otherwise a list.
}


}

\value{
return a list that contains the status of the solution , the objective value and the 
fluxes that are free of cycles. 
}

\author{Abdelmoneim Amer Desouki}

\seealso{
\code{\link{modelorg}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
data(iAF1260)
model=iAF1260
fba=optimizeProb(model)
cfopt=cfFBA(model,wtflux=getFluxDist(fba),objVal=lp_obj(fba), retOptSol=FALSE)
llflx=cfopt$fluxes
flx=getFluxDist(fba)
\dontrun{
layout(matrix(c(1,2,3,1,2,3), 2, 3, byrow = TRUE))
hist(log10(abs(flx[abs(llflx-flx)>1e-3])),main="loop flux",col="lightblue")
hist(log10(abs(llflx[abs(llflx-flx)>1e-3])),main="after removing loops",col="orange")
hist(log10(abs(flx[abs(llflx-flx)<1e-3])),main="fluxes not in loops",col="lightgreen")
}
## The function is currently defined as
"cfFBA"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cycle free flux }
\keyword{ loopless FBA }% __ONLY ONE__ keyword per line
