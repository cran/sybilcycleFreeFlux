cfFVA <- function(model,rxnList,solver = SYBIL_SETTINGS("SOLVER"),
                  pct_objective=100
				  ,solverParm=NA
				  ,verboseMode = 2,includeRxnEqn=TRUE
				  ,boundFlg=FALSE  # should be set to FALSE to enumrateCycles
				  ){
	#when a loop exists set rxn(s) of no essential flux to 0: except rxn to be maximized 
	# if more than one rxn can be a problem!!!
	## send obj value to cfFBA
#function to get reaction equation
	getRxnEqn <- function(rxn=1){
  	m=S(model)[,rxn]
  	# input ==> output
  	mcf=ifelse(abs(m)>1,paste("(",abs(m),")",sep=""),"")
  	eqn=paste(gsub(" "," + ",Reduce(paste,paste(mcf[m<0],met_id(model)[which(m<0)],sep=""))),
  		ifelse(react_rev(model)[rxn],"<==>","-->"),
  	    gsub(" "," + ",Reduce(paste,paste(mcf[m<0],met_id(model)[which(m>0)],sep=""))) ,sep=" ")
  	
  return(eqn)
}
	tol=1e-6;
	
	MAX_CYCLE_FACTOR= (10 ^ floor((log10(max(uppbnd(model)))-log10(tol))/2)) #:
	# cycle_factor: max(flx)/min(flx) through cycle, it should be fixed for all fluxes through the cycle
	FVA_TOL = 0.001
	MIN_LOOP_FLUX=max(uppbnd(model))/MAX_CYCLE_FACTOR;# # 1000/(F1sn-F2sn)  (F3-F2) > 1000/MAX_CYCLE_FACTOR
	# F1 will be max possible(F1-F2) > MAX(F3)/MAX_CYCLE_FACTOR
				# as r is max there must be at least one rxn in the loop at 1000 so (F3-F2)/1000 > tol
		# #WRONG:there is always a zero because the two ExPw should be diff in at least 1 rxn
		
	sf=optimizeProb(model,solver=solver,solverParm=solverParm);
	bmrxn=which(obj_coef(model)!=0);# should be !=0
	objVal=getFluxDist(sf)[bmrxn]*pct_objective/100;#15/3/2015 add pct_objective
	
#####-------------------------########----------------------#######
maxFlx=NULL;
lowbnd(model)[bmrxn]=objVal-tol;#better to set one bound, problem in infeasibility [4/6/2015]
#uppbnd(model)[bmrxn]=objVal;#assume lpdir=max [14/4/2015]

#excReact = findExchReact(model)[1];# 1 is position
#excReactPos=react_pos(excReact$exchange);
excReact = findExchReact(model);
excReactPos=which(react_id(model) %in% react_id(excReact));
	
#mod=sysBiolAlg(model,solver=solver);
#prob_r=problem(mod)
model_r=model
#numerical issues  TODO
#lb=lb-tol*abs(lb)
#ub=ub+tol*abs(ub)

nCols=react_num(model);
# for tracing^^
if(includeRxnEqn){
	rxneqns=sapply(c(1:nCols),function(x) getRxnEqn(x))
}else{
	rxneqns=rep(NA,nCols);
}
res=NULL;
 #model=model_r
 #mod=sysBiolAlg(model,solver=solver);
		#prob=problem(mod)
	
 loops=NULL;
for (rxn in rxnList){
	sn=which(react_id(model)==rxn);
	if(length(sn)==0){
		stop(sprintf("Bad reaction id: %s",rxn));
	}
	for(dirxn in c(1,-1)){
		stime=proc.time();
		optid=1
		opts=list()
		brklist=data.frame(parent_id=0,brkrxns="",brktree="",stringsAsFactors=F)#parent_id,brkrxns(model order),brktree(seq of ops):[Searchable table]
		#ocf=obj_coef(model)
		#backupProb(mod) to restore original copy, other solver may have time issues
		##lb=lowbnd(model);	ub=uppbnd(model);
		bestVal=bestVal.proc=ifelse(dirxn==1,lowbnd(model)[sn]-1,uppbnd(model)[sn]+1); # the worst values
		bestVal.optid=0
		opts[[1]]=list(lb=lowbnd(model_r),ub=uppbnd(model_r),level=1,brkrxns=NULL,NLFlx=bestVal)# restore from model_r
		boundedsteps=0;
		#changeColsBnds(prob,j=c(1:nCols),lb=lb,ub=ub)
		#for(i in 1:5){# max no of loops
		while(optid<=length(opts)){
			if(verboseMode>2)
				print(sprintf("rxn:%s iteration:%d direction:%d",rxn,optid,dirxn));	
			# 0-Check if result can be better than the bestVal acheived so far
			###  Bound step -----------------
			
			#1- Max flux through rxn
			obj_coef(model)[bmrxn]=0;		obj_coef(model)[sn]=dirxn;#min/max
			
			#ocf[bmrxn]=0; ocf[sn]=1;
			#solL=optimizeProb(mod,react=c(1:nCols),obj_coef=ocf,resetChanges =TRUE)
			#changeObjCoefs(prob,c(sn,bmrxn),c(dirxn,0))
			#solL=optimizeProb(mod);#,solver=solver:solver defined in opt_obj
			
			#start simulations			
			lowbnd(model)=opts[[optid]]$lb;     uppbnd(model)=opts[[optid]]$ub;
			
			solL=optimizeProb(model,solver=solver,solverParm=solverParm)
			if(lp_ok(solL) != 0 ||  length(checkSolStat(stat=lp_stat(solL),solver=solver))!=0 ) {#check if solution is feasible
				if(is.na(checkSolStat(stat=lp_stat(solL),solver=solver))) stop("couldn't check solution status");
				stop(sprintf("Solution is not optimal, solL, rxn: %d, status:%d solver:%s",sn,lp_stat(solL),solver))
			}
		
			F1=getFluxDist(solL)#$fluxes#[solL$fldind];
			if(optid==1){
				F1org=F1;  
			} else # the only possible bound step (as rxn is maximized(F2[sn]:is not the uppbnd of rxn)
				if( ifelse(dirxn==1,bestVal - F1[sn], F1[sn] - bestVal ) > 0 && boundFlg && F1[sn] != bestVal.proc){
				#must be signed diff [if the bestVal not processed
					if(verboseMode>2) print(sprintf("sn:%d dirxn: %d, Bound F1, bestVal:%f, F1=%f",sn,dirxn,bestVal,F1[sn]));
					boundedsteps=boundedsteps+1
					optid=optid+1;  # noneed to calc F2 as it will not give better result
					next;
				}
			#model=model_r;
			# 2-F1 without loops
			obj_coef(model)[sn]=0;		obj_coef(model)[bmrxn]=obj_coef(model_r)[bmrxn];
			#changeObjCoefs(prob,c(sn,bmrxn),c(0,1))
			#print(getObjDir(prob))
			#browser()
			
			if(abs(F1[sn])>tol){#1e-7
				if(!(sn %in% excReactPos)){#send xchng rxn list to avoid recalc.
						solNL=cfFBA(model,wtflux=F1,excReactPos=excReactPos,solver=solver,objVal=F1[bmrxn],retOptSol=FALSE,
						    tol=tol,solverParm=solverParm);	
				}else{#fix exchange except one rxn
						excReactPos1=excReactPos[-which(sn %in% excReactPos)]
						solNL=cfFBA(model,wtflux=F1,excReactPos=excReactPos1,fixExchRxn=TRUE,solver=solver,objVal=F1[bmrxn],
							tol=tol,solverParm=solverParm,retOptSol=FALSE);	
				}
				if( solNL$ok != 0 ||  length(checkSolStat(stat=solNL$stat,solver=solver))!=0 ) {#check if solution is feasible
					if(is.na(checkSolStat(stat=solNL$stat,solver=solver))) 
						stop("couldn't check solution status, F2, cfFBA");
				
				write.csv(file=sprintf("cfFVA_solNL_%d_%d_%s.csv",sn ,solNL$stat,solver),
							cbind(react_id(model),lb=lowbnd(model),ub=uppbnd(model),cflb=solNL$lb,cfub=solNL$ub,F1))
				
					print(sprintf("Solution is not optimal, solNL, rxn: %d, status:%d (%s), solver:%s",
						sn,solNL$stat,getMeanStatus(solNL$stat,solver),solver))
					browser()
				}
				F2=solNL$fluxes;#F2 represents F1 without loops
			}else{F2=F1}
				
			# browser()
			###  Bound step ----------------- F2 is not the max of r (use costcoef to select favor r) but it can be a value
			# if( ifelse(dirxn==1,bestVal - F2[sn], F2[sn] - bestVal ) >= 0 && optid > 1 && boundFlg){ #must be signed diff
				## boundFlg is used to disable boundFlg step(when enumCycles)
				# optid=optid+1;
				# if(verboseMode>2) print(sprintf("sn:%d dirxn: %d, Bound F2, bestVal:%f, F2=%f",sn,dirxn,bestVal,F2[sn]));
				# next;
			# }
			### -----
			bestVal=F2[sn]  # best NLFlx (not yet achieved)
			# bestVal.processed
			### ------------------------
			
			#3-remove all loops except the one containing [rxn](Identify loop)
			if( abs(F1[sn]-F2[sn]) > FVA_TOL){# if a loop exists
				if(verboseMode>2)
					print(c("diff:" , F1[sn],F2[sn])) ;
				#lowbnd(model)[sn]=F1[sn];	uppbnd(model)[sn]=F1[sn]; # Force rxn to go with max flux to identify loop
				if(!(sn %in% excReactPos)){#fix rxn sn
						 sol1L=cfFBA(model,wtflux=F1,fixExchRxn=TRUE,excReactPos=c(sn,excReactPos),solver=solver,objVal=F1[bmrxn],retOptSol=FALSE);
				}else{
						excReactPos1=excReactPos[-which(sn %in% excReactPos)]
						sol1L=cfFBA(model,F1,excReactPos=excReactPos1,solver=solver,objVal=F1[bmrxn],retOptSol=FALSE);
				}
				if( sol1L$ok != 0 ||  length(checkSolStat(stat=sol1L$stat,solver=solver))!=0 ) {#check if solution is feasible
					if(is.na(checkSolStat(stat=sol1L$stat,solver=solver))) 
						stop("couldn't check solution status, F3, cfFBA");
					
					write.csv(file=sprintf("cfFVA_sol1L_%d_%d_%s.csv",sn ,sol1L$stat,solver),
							cbind(react_id(model),lb=lowbnd(model),ub=uppbnd(model),cflb=sol1L$lb,cfub=sol1L$ub,F1,F2))
					
					print(sprintf("Solution is not optimal, sol1L, rxn: %d, status:%d (%s), solver:%s",
						sn,sol1L$stat,getMeanStatus(sol1L$stat,solver),solver))
					browser()
				}

				if(verboseMode>2)
					print(sprintf("sol F3: stat=%d, obj=%f",sol1L$stat,sol1L$obj));
				F3=sol1L$fluxes#[sol1L$fldind];Only one Loop
				#write.csv(file=sprintf("fl_iter%d_%s.csv",i,rxn),cbind(rxn=react_id(model),rxneqns,lb,ub,F1,F2,F3))
				
				lp=(abs(F3-F2) > MIN_LOOP_FLUX );## detected loop (same flx rot?)
				flxrot=abs(F1[sn]-F2[sn])#nonessential flux through rxn[sn]
				cycle_factor=max(abs(F3[lp]-F2[lp]))/min(abs(F3[lp]-F2[lp])) # not F1-F2: more than 1 loop
					if(cycle_factor > MAX_CYCLE_FACTOR ) print(sprintf("cycle factor: %f",cycle_factor));
				#browser()
				if(verboseMode>2){
					print("Rxns in loop: ");
					print(cbind(which(lp),react_id(model)[lp],F3=F3[lp],F2=F2[lp],sFVA=F1[lp]));
				}
				#4- break loop and maximize again
				# set rxns going to 0 to 0 and maximize again
				#model=model_r; undo last change in model but keep changes thru the crnt rxn iteration
				#lowbnd(model)[sn]=lowbnd(model_r)[sn]; uppbnd(model)[sn]=uppbnd(model_r)[sn];
				
				brk=abs(F3-F2) > MIN_LOOP_FLUX & (abs(F2) <= tol) ; # rxn in loop goes to zero
				
				if(brk[sn] && sum(brk)==1) {# crnt rxn the only rxn having Zero flux without loop
					maxFlx=rbind(maxFlx,cbind(rxn,iter=optid,dirxn,initflx=F1org[sn],LoopyMax=F1[sn],NLMax=F2[sn],lp=paste(react_id(model)[lp],collapse=",")
					,brkSeq=brklist[optid,"brktree"],#,brk=paste(react_id(model)[brk],collapse=",")
					loopyflx=paste(F1[lp],collapse=","),NLflx=paste(F2[lp],collapse=","),LPF3flx=paste(F3[lp],collapse=","),
					rvrs=paste(react_rev(model)[lp],collapse=","),flxrot=abs(F1[sn]-F2[sn]),cycle_factor=cycle_factor,llen=sum(lp),
					eqn=paste(rxneqns[lp],collapse=",") ));
					
					bestVal.proc=0
					#bestVal=0 # used in bound step 4/6/2015
					#break;
				}else{
					brk[sn]=FALSE;
					if(sum(brk)==0){
							write.csv(file=sprintf("Error_cfFVA_rxn%s_%s.csv",react_id(model)[sn],ifelse(dirxn==1,"max","min")),
							  cbind(lb=lowbnd(model),ub=uppbnd(model),lp,F1,F2,F3))
							 # tobe stop
							 print(sprintf("No rxn will go to zero in this Loop! rxn:%s, dirxn:%d,flxrot:%f, lp:%s",rxn,dirxn,abs(F1[sn]-F2[sn]),
						      paste(react_id(model)[lp],collapse=",")));
							   browser();
					}
				
					loops=rbind(loops,cbind(rxn,iter=optid,dirxn,initflx=F1org[sn],LoopyMax=F1[sn],NLMax=F2[sn],lp=paste(react_id(model)[lp],collapse=",")
					,newbrk=paste(react_id(model)[brk],collapse=","),crntbrk=brklist[optid,"brkrxns"],crntlevel=opts[[optid]]$level,
					loopyflx=paste(F1[lp],collapse=","),NLflx=paste(F2[lp],collapse=","),LPF3flx=paste(F3[lp],collapse=","),
					rvrs=paste(react_rev(model)[lp],collapse=","),flxrot=abs(F1[sn]-F2[sn]),cycle_factor,llen=sum(lp),eqn=paste(rxneqns[lp],collapse=",") ));
					
					
					opts_len=length(opts)+1;
					for(r in which(brk)){ #lb[brk]=ub[brk]=F2[brk]; #Add new optimizations to opts		
					   if(verboseMode>3) print(sprintf("Add opt to break lp through rxn:%s, value:%f",react_id(model)[r],F3[r]));
					   v_newbrkrxns=c(opts[[optid]]$brkrxns,r)
					   v_tmp=paste(react_id(model)[v_newbrkrxns],collapse=",")
					   if(v_tmp %in% brklist[,"brkrxns"]){#skip opt as it is repeated
							if(verboseMode>2) print(c("repeated list:",v_tmp,optid,react_id(model)[r]));
					   }else{
							opts[[opts_len]]=list(lb=opts[[optid]]$lb,ub=opts[[optid]]$ub,level=opts[[optid]]$level+1,
								brkrxns=v_newbrkrxns,NLFlx=F2[sn])
							#close only one direction of the reaction[CHK:can rxn generate loops if run in any dirxn??]
							#store sign(F3[r])
							if(F3[r]<0){
								opts[[opts_len]]$lb[r]=0
							}else{
								opts[[opts_len]]$ub[r]=0
							}
							brklist[opts_len,]=data.frame(parent_id=optid,brkrxns=v_tmp,brktree=paste(brklist[optid,"brktree"],react_id(model)[r],sep=",")
								,stringsAsFactors=FALSE)
							#browser()
							opts_len=opts_len+1
					   }
					}#end for
				}
			}else{		 # max value is ok (there is no loop)
				if (optid>1) {
					print(c(rxn,optid,F1[sn],F2[sn]));
				}
				maxFlx=rbind(maxFlx,cbind(rxn,iter=optid,dirxn,initflx=F1org[sn],LoopyMax=F1[sn],NLMax=F2[sn],lp="",
								brkSeq=brklist[optid,"brktree"],loopyflx=NA,NLflx=NA,
				             LPF3flx="",rvrs=react_rev(model)[sn],flxrot=0,cycle_factor=NA,llen=0,eqn=""));
							 
				if(ifelse(dirxn==1,bestVal.proc - F2[sn], F2[sn] - bestVal.proc ) >= 0) bestVal.proc=F2[sn]; # used in bound step
				#break;
			}
			
			optid=optid+1;
		}#iteration for loop
		#find max/min
		rxnFlx=maxFlx[maxFlx[,1]==rxn & maxFlx[,3]==dirxn,,drop=FALSE]
		if(length(rxnFlx)==0) { print("not processed"); browser()}
		if(dirxn==1){
			maxid=which.max(as.numeric(rxnFlx[,"NLMax"]));
			maxVal=rxnFlx[maxid,"NLMax"];
			sFVAmax=rxnFlx[maxid,"initflx"];
			brkSeqMax=rxnFlx[maxid,"brkSeq"];
			#brkFlxMax=maxFlx[maxid,"brklistflx"]
			maxiter=optid-1
			maxTime=(proc.time()-stime)["elapsed"];
			maxBndSteps=boundedsteps
		}else {
			minid=which.min(as.numeric(rxnFlx[,"NLMax"]));
			minVal=rxnFlx[minid,"NLMax"];
			sFVAmin=rxnFlx[minid,"initflx"];
			brkSeqMin=rxnFlx[minid,"brkSeq"];
			#brkFlxMin=maxFlx[minid,"brklistflx"]
			miniter=optid-1
			minTime=(proc.time()-stime)["elapsed"];
			minBndSteps=boundedsteps
		}
 }#dirxn
 res=rbind(res,cbind(rxn=rxn,maxiter,miniter,minVal,maxVal,sFVAmin,sFVAmax,brkSeqMin,brkSeqMax,minTime,maxTime,minBndSteps,maxBndSteps))
 }#rxn
 
return(list(res,maxFlx,loops));
}