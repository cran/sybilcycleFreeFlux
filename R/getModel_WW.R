getModel_WW <- function(model,solver = SYBIL_SETTINGS("SOLVER") ){
solver='sybilGUROBI'
## This preprocessing part is based on the algorithm given in Wright & Wagner 2008
#preprocessing to get
loopicious=NULL;

t1=proc.time()
	mirv=mod2irrev(model)
t2=proc.time()
	excr=findExchReact(mirv)
	isExch=(react_id(mirv) %in% react_id(excr))
	#intRxn=which(!isExch)
	#remove Exchange rxns
	mirv=rmReact(mirv,react=which(isExch));
	de_met=which(rowSums(S(mirv)!=0)==1);
	 blockedRxn=apply(S(mirv)[de_met,],1,function(x){return (which(x!=0))})
	intRxn=c(1:react_num(mirv))
	lowbnd(mirv)[c(1:react_num(mirv))]=0;
	
	print(sprintf("number of internal reactions(irrev): %d",length(intRxn)))
	sols=NULL
t3=proc.time()	
	for(r in intRxn){
		print(r)
		uppbnd(mirv)[c(1:react_num(mirv))]=1000;
		mets1=which(S(mirv)[,r]!=0)
		candRxns=intRxn[colSums(S(mirv)[mets1,intRxn]!=0)==length(mets1) & intRxn!=r & colSums(S(mirv)[,intRxn]!=0)==length(mets1)]
		if(length(candRxns)>0){
			lp2=logical(length(candRxns))
			for(i in 1:length(candRxns)){
				if(!all(S(mirv)[,candRxns[i]]==-S(mirv)[,r])) {
					lp2[i]=FALSE ;
				}else{
					lp2[i]=TRUE;
				}				
			}
			uppbnd(mirv)[candRxns[lp2]]=0;
		}
		obj_coef(mirv)[c(1:react_num(mirv))]=0;
		obj_coef(mirv)[r]=1;
		sol=optimizeProb(mirv,solver=solver);
		sols=rbind(sols,cbind(r,stat=lp_stat(sol),obj=lp_obj(sol),RRlen=length(candRxns),flx=getFluxDist(sol)[r],ub0=sum(uppbnd(mirv)==0)))
		#print(sol)
		if(lp_obj(sol)>1e-7) loopicious=c(loopicious,r);
	}	

t4=proc.time();	

#the number of reactions in nontrivial loops
print(sprintf("Number of reactions in new model:%d",length(loopicious)))

#write.csv(file="WW_loopicious.csv",cbind(sn=sols[sols[,3]!=0,1],rxn=react_id(mirv)[sols[sols[,3]!=0,1]],obj=sols[sols[,3]!=0,3]))

#write.csv(file="WW_preproc.csv",cbind(sn=sols[,1],rxn=react_id(mirv)[sols[,1]],obj=sols[,3],flx=sols[,5],RRlen=sols[,4]))



lpModel=rmReact(mirv,react=react_id(mirv)[-loopicious])
##cbind(react_id(lpModel),lowbnd(lpModel),uppbnd(lpModel))
 return(lpModel)
}