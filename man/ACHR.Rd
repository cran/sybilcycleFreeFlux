\name{ACHR}
\alias{ACHR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random Sampling of Solution Space
}
\description{
implements sampling algorithm
}
\usage{
ACHR(model, W = 2000, nPoints = 5000, stepsPerPoint = 10, 
solver = SYBIL_SETTINGS("SOLVER"), method = SYBIL_SETTINGS("METHOD"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
An object of class \code{\link{modelorg}}.
}
  \item{W}{
Number of warmup points. It should be more than double the number of reactions of 
the model.
}
  \item{nPoints}{
Number of points to be generated
}
  \item{stepsPerPoint}{
number of steps per point, default is 10 steps.
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
 
}
\details{
Starts by calculating warm up points
}

\author{Abdelmoneim Amer Desouki}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{modelorg}}
\code{\link{cfFBA}}
}
\examples{
\dontrun{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
library(sybilcycleFreeFlux)
data(Ec_core);
model=Ec_core;
solver="cplexAPI"
W=500
nPnts=1000
s1=ACHR(model,W,nPoints=nPnts,stepsPerPoint=10)


sFVA=fluxVar(model,solver=solver) 
fva_min=sFVA@lp_obj[(c(1:length(react_id(model))))];
fva_max=sFVA@lp_obj[(c((length(react_id(model))+1):(2*length(react_id(model)))) )];
table(lp_stat(sFVA))

pnts=s1$Points
fvamin=apply(pnts,1,min)
fvamax=apply(pnts,1,max)

write.csv(file="fva.csv",cbind(fva_min,fvamin,fva_max,fvamax,lb=lowbnd(model),
ub=uppbnd(model)))
#####Plot samples
bmrxn=which(obj_coef(model)==1)
bmrow=S(model)[bmrxn,]

objvals=NULL
solver="glpkAPI"
nRxns=react_num(model);
llpnts=  matrix(rep(0,nRxns*nPnts),ncol=nPnts);
for(i in 1:nPnts){
	objvals=rbind(objvals,obj= pnts[bmrxn,i])
	lrf=lrFBA(model,wtflux=pnts[,i],solver=solver,objVal= pnts[bmrxn,i])
	llpnts[,i]=lrf$fluxes;
#Sys.time()
	print(sprintf("point \%d:\%f",i,objvals[i]))
}
llfvamin=apply(llpnts,1,min)
llfvamax=apply(llpnts,1,max)

 write.csv(file="objv.csv",objvals)
write.csv(file="llfva.csv",cbind(fva_min,llmin=llfvamin,fva_max,llmax=llfvamax,fvamin,
fvamax,lb=lowbnd(model),ub=uppbnd(model)))
nloopflux=NULL
loopflxll=NULL
loopflxlp=NULL

for(i in (1:length(react_id(model))))
	for(j in (1:nPnts)){
	#print(c(i,j))
		if(abs(pnts[i,j]-llpnts[i,j])<1e-7){
			nloopflux=c(nloopflux,pnts[i,j])
		}else{
		loopflxll=c(loopflxll,llpnts[i,j])
		loopflxlp=c(loopflxlp,pnts[i,j])
		}
	}
layout(matrix(c(1,2,3,1,2,3), 2, 3, byrow = TRUE))
hist(log10(abs(loopflxlp)),col="lightblue",main="a-loop fluxes",xlim=c(-3,3),
xlab="Log10(flux)")
hist(log10(abs(loopflxll)),col="orange",main="b-using cycleFreeFlux",
xlim=c(-3,3),xlab="Log10(flux)")
hist(log10(abs(nloopflux)),col="lightgreen",main="c-non-loop fluxes",
xlim=c(-3,3),xlab="Log10(flux)")

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Random Sampling }
\keyword{ loopless Sampling }% __ONLY ONE__ keyword per line
