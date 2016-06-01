\name{getModel_WW}
\alias{getModel_WW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
 based on the algorithm given in Wright & Wagner 2008
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
get a model that contains only reactions in loops.
}
\usage{
getModel_WW(model, solver = SYBIL_SETTINGS("SOLVER"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
An object of class \code{\link{modelorg}}.
}
  \item{solver}{
Single character value.  The solver to use.  See
    \code{\link{SYBIL_SETTINGS}} for possible values.\cr
    Default: \code{SYBIL_SETTINGS("SOLVER")}.
}
}

\value{
return an object of class \code{\link{modelorg}} containing subset of 
reactions of the original model that participate in at least one loop.
}
\references{
Wright, J. and Wagner, A. (2008). Exhaustive identification of steady state cycles in
large stoichiometric networks. BMC systems biology, 2, 61.
}

\author{Abdelmoneim Amer Desouki}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{modelorg}}
\code{\link{cfFBA}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
\dontrun{
 data(iAF1260)
 loopiAF=getModel_WW(iAF1260)
 }
 }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
