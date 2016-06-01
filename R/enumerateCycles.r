enumerateCycles <- function(model,rxnList,solver = SYBIL_SETTINGS("SOLVER") ){
fva=cfFVA(model,rxnList,solver)
ldet=fva[[3]]
lp=ldet[ldet[,"llen"]>0,,drop=FALSE]
ulp=unique(lp[,c("lp","llen")])
return (ulp)
}