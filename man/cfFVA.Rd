\name{cfFVA}
\alias{cfFVA}
\alias{lrFVA}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cycle free flux variability
}
\description{
finds flux variability without loops. Maximize individual reactions and test if loop exists using \code{\link{cfFBA}}.
Then breaks the loop by setting the fluxes that goes to zero in loopless flux to zero.
}
\usage{
cfFVA(model, rxnList, solver = SYBIL_SETTINGS("SOLVER"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
An object of class \code{\link{modelorg}}.
}
  \item{rxnList}{
a character vector of reaction ID's to find their flux variability
  }
  \item{solver}{
Single character value.  The solver to use.  See
    \code{\link{SYBIL_SETTINGS}} for possible values.\cr
    Default: \code{SYBIL_SETTINGS("SOLVER")}.

}
}

\value{
return two lists res and maxFlx. The first list contains reactions with the computed values. 
The second list contains the details of calculations. 
If a reaction is involved in a loop it will appear at least twice in the second list.
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

fva=cfFVA(model,react_id(model))
write.csv(file="cfFVA1_cnd_res.csv",fva[[1]]);
write.csv(file="cfFVA1_cnd_det.csv",fva[[2]]);
#plot fluxes

## The function is currently defined as
"cfFVA"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Flux variability }
\keyword{ cycle free flux }% __ONLY ONE__ keyword per line
