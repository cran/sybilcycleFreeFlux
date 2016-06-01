cfFBA <- function(model
,wtflux  # initial flux distribution when wtflux is NA don't include its constraint 
	,objVal = NA #min objval
	,fixExchRxn=TRUE
	,excReactPos=NULL
	 ,lpdir = SYBIL_SETTINGS("OPT_DIRECTION")
	 ,solver = SYBIL_SETTINGS("SOLVER")
         ,method = SYBIL_SETTINGS("METHOD")
		 ,tol=SYBIL_SETTINGS("TOLERANCE")
         ,solverParm=NA
		,verboseMode = 2
		,safeBounds=FALSE # can be used to deal with false infeasibility state returned by solvers
######### ADDED BY GABRIEL ###############
		,retOptSol = TRUE
##########################################
) {
# fix the direction of rxn, are not allowed to increase 
# minimize abs(flux)
# change lower and upper bounds only and obj fun
# forcing the same direction is a problem in some loops reversible rxns flux is c-1000
#N.B:23/10/2012, Keep exchange reactions fixed
#   considering cases of no objective or no exchange rxns
#	4/6/2015: sysBiolAlg made problem with 90% obj rintala, sim_all, cnd 1,'CEM1_1', max,
#   5/6/2015 numerical issues in solvers GUROBI and CPLEX, parameter feasibilityTol\CPX_PARAM_EPRHS=1e-9 must be set to very low value
#		parameter safeBounds = TRUE adjusts the values of bounds to deal with this issue

        if (! is(model, "modelorg")) {
            stop("needs an object of class modelorg!")
        }  
        
		nCols=react_num(model);
        # force objective constraint in model
		# mod_ocfs=getObjCoefs(prob,c(1:nCols))
		mod_ocfs=obj_coef(model)#getObjCoefs(prob,c(1:nCols))
		bmrxn=(mod_ocfs!=0); #Adjust objective coefficient to include values other than 1, 4/4/2015
	#it may be not important to fix obj value 
      if(is.na(objVal) && length(bmrxn)!=0){
      	sol=optimizeProb(model,solver=solver,solverParm=solverParm);#solver included in problem
      	objVal=sol$obj;
      }      
#####
mod_lb=lb=lowbnd(model);
mod_ub=ub=uppbnd(model);
ocf=rep(0,nCols);
# lb=max(min(wtflx,0),lb)    ub=min(max(0,wt),ub)
for (i in 1:nCols ){
	if(!is.na(wtflux[i]) ){  #6/4/2015: bounds less than tolerance return NO SOLUTION
		if(wtflux[i]<0){ # flux is -ve
			#bounds of model shouldn't be relaxed to the wtflux, correcting this will result in some infeasible solutions
			lb[i]=wtflux[i]; # old lb must be -ve also and |lb|>=|wtflx| 
			# lb[i]=max(lb[i],wtflux[i]);#17/4/2015 # gave infeasible solution because of solver numerical issues
			ub[i]=min(0,ub[i]); # old ub may be -ve, should not be relaxed
		}
		else #if(wtflux[i]>=0)
		{
			lb[i]=max(0,lb[i]);# made a problem in ATPM, should not be relaxed to 0
			ub[i]=wtflux[i];
			# ub[i]=min(ub[i],wtflux[i]);			
		}
		ocf[i]=sign(wtflux[i]);
	}
}        

        
if(fixExchRxn){  # Fixing exchange rxn from changing may preserve some rxns with biological evidence from changing
	if(is.null(excReactPos) && is(model, "modelorg")){
		excReact = findExchReact(model);
		if(!is.null(excReact)){
			excReactPos=react_pos(excReact );#excReactPos=which(react_id(model) %in% react_id(excReact));#excReact$exchange
		}else{
			excReactPos=NULL;
		}
	}#6/4/2015: fixed bug: NA value to bounds  | abs(wtflux[excReactPos]) < tol
 	if(length(excReactPos)>0){
		lb[excReactPos]=ifelse(is.na(wtflux[excReactPos]) ,lb[excReactPos],wtflux[excReactPos]) 
		ub[excReactPos]=ifelse(is.na(wtflux[excReactPos]) ,ub[excReactPos],wtflux[excReactPos])
	}
}
# Fix biomass: fix only one direction(not difference if optimal)#
#if(lpdir=="max"){lb else ub}

#Test not fixing objective value 5/6/2015 Needs testing
if(length(bmrxn)>0){
	lb[ bmrxn ] = objVal;
	ub[ bmrxn ] = objVal;
}

 # origObjDir=getObjDir(prob);	
 # setObjDir(prob,"min");
 # sol=optimizeProb(model,lpdir="min",react=c(1:nCols),lb=lb,ub=ub,obj_coef =ocf,retOptSol=FALSE)
 obj_coef(model)=ocf
 
 if(safeBounds){ # to deal with numerical rounding of fluxes that are calculated from high precision
	lb = lb - tol/2
	ub = ub + tol/2
}

 lowbnd(model)=lb 
 uppbnd(model)=ub 
 
 # write.csv(file="test_cfFBA.csv",cbind(react_id(model),mod_lb,mod_ub,lb,ub,ocf,mod_ocfs,wtflux))
 
 sol=optimizeProb(model,lpdir="min",solverParm=solverParm,solver=solver,retOptSol=FALSE)
 # setObjDir(prob,origObjDir); 

# 
 #--------------Output-------------------------

######### ADDED BY GABRIEL ###############
if (isTRUE(retOptSol)) {
            optsol <- new("optsol_optimizeProb",
                          mod_id       = mod_id(model),
                          mod_key      = mod_key(model),
                          solver       = solver,
                          method       = method,
                          algorithm    = "fba",#algorithm(mod),
                          num_of_prob  = 1L,
                          lp_dir       = lp_dir,
                          lp_num_rows  = met_num(model),
                          lp_num_cols  = nCols,
                          lp_ok        = as.integer(sol[["ok"]]),
                          lp_obj       = sol[["obj"]],
                          lp_stat      = as.integer(sol[["stat"]]),
                          obj_coef     = ocf,
                          obj_func     = printObjFunc(model),
                          fldind       = c(1:nCols),
                          fluxdist     = fluxDistribution(fluxes = sol[["fluxes"]],
                                                    nrow = length(sol[["fluxes"]]),
                                                    ncol = 1L),
                          alg_par      = solverParm #alg_par(mod) check!!
						  )
 
}
else{
        optsol <- list(ok = sol$ok,
	                   obj = sol$obj,
	                   stat =sol$stat,
	                   fluxes = sol$fluxes,
	                   wtflx=wtflux,
	                   lb=lb,
	                   ub=ub,
	                   ocf=ocf
	                  )	    

}
##########################################


    return(optsol)        
 }
