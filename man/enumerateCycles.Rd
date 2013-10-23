\name{enumerateCycles}
\alias{enumerateCycles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to enumerate cycles
}
\description{
Uses cycleFreeFlux to enumerate cycles in a given metabolic network. 
 }
\usage{
enumerateCycles(model, rxnList, solver = SYBIL_SETTINGS("SOLVER"))
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
return a list of unique loops found in the model.
}


\author{Abdelmoneim Amer Desouki}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{modelorg}}
\code{\link{lrFBA}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
data(Ec_core)
model=Ec_core

cycles=enumerateCycles(model,react_id(model))

## The function is currently defined as
"enumerateCycles"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
