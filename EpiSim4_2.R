EpiSim4 <- function(
    nind,
    # For explicitly specifying group and index case status of each individual.
    # Requires a data.frame with 'nind' rows and columns for group (character),
    # donor/recipient ("D" or "R") and individual id
    group.and.index = NULL, 
    ngroups = 1,
    PairwiseBeta, RRgamma, 
    pop = NULL,
    offspring.only = FALSE,
    GRM, GRM.from.pedigree = FALSE,
    ntraits = 3, tmeans = rep(0, 3),
    tVA, trhoG, tVE, trhoE,
    traitnames = c("susceptibility", "infectivity", "recoverability"),
    trait.transform = c("lognormal", "invcloglog"),
    frailties = c("lognormal", "gamma"),
    Epi.type = c("SI", "SIS", "SIR")
) {
    
    gr.size <- nind / ngroups
    
    # Simulate trait values
    # call existing function
    # --------------------------------------------------------------------------
    
    if (is.null(pop)) {
        pop <- make.tvals.int(
            nind = nind,
            GRM = GRM,
            GRM.from.pedigree = GRM.from.pedigree,
            ntraits = ntraits,
            tmeans = tmeans,
            tVA = tVA,
            trhoG = trhoG,
            trhoE = trhoE,
            tVE = tVE,
            traitnames = traitnames,
            trait.transform = trait.transform
        )
    } else {
        pop <- pop
        
        pop[, frailties := frailties]
        
        if (offspring.only) {
            pop <- pop[dso == "offspring"]
        }
        
    }

    pop[, PairwiseBeta := PairwiseBeta]

    # Assign groups and index cases
    # Call existing function
    # --------------------------------------------------------------------------
    
    if (is.null(group.and.index)) {
        pop <- group.assign(pop = pop, ngroups = ngroups)
        pop[, frailties := frailties]
        pop[, Epi.type := Epi.type]
    } else {
        # Custom assignment of groups and index cases
        setkey(pop, id)
        setkey(group.and.index, individual_ForPedigree)
        pop = pop[group.and.index]
        setnames(pop, "Box_recoded", "group")
        pop[DR == "D", IndexCase_0 := 0]
        pop[is.na(IndexCase_0), IndexCase_0 := 1]
        pop[IndexCase_0 == 0, Tinf := 0]
        
        
        # setkey(pop, id)
        # setkey(group.and.index, individual_ForPedigree)
        # #DON'T WANT TO REORDER THE ROWS!!!
        # pop[, group := group.and.index[, Box_recoded]]
        # pop[id %in% group.and.index[DR == "D", individual_ForPedigree], IndexCase_0 := 0]
        # pop[is.na(IndexCase_0), IndexCase_0 := 1]
        # pop[IndexCase_0 == 0, Tinf := 0]
        
    } # will need to generalize this later
    
    # Generate epidemics
    # --------------------------------------------------------------------------
    
    # Store everything in a table indicating group and epidemic
    # Basically recreate POP and sirTS but with all runs in a single table
    
    # I = infected, R = recovered
    pop[Epi.type == "SIR", c("I", "R") := list(1 - IndexCase_0, 0)]
    # I = infected
    pop[Epi.type != "SIR", I := 1 - IndexCase_0]
    
    # Each individual has to have a contact matrix*******************************#DO THIS LATER#
    
    if (is.null(group.and.index)) {
        
        if (Epi.type == "SIR") {
            SIRts <- data.table(
                group = unique(pop[, group]),
                time = 0,
                S = gr.size - 1,
                I = 1,
                R = 0,
                totS = nind - ngroups,
                totI = ngroups,
                totR = 0
            )
            nextT = 0
            
        }#Endif
        
        if(Epi.type=="SI"){
            SIRts = data.table(
                group=unique(pop[,group]),
                time=0,
                S=gr.size-1,
                I=1,
                totS=nind-ngroups,
                totI=ngroups);
            nextT = 0;
        }#Endif
        
    }else{#end of; if(is.null(group.and.index))
        
        if(Epi.type=="SI"){
            SIRts = data.table(
                group=unique(pop[,group]),
                time=0);
            #*************************************
            
            S=pop[IndexCase_0==1,.N,by=group];setnames(S,"N","S");
            I=pop[IndexCase_0==0,.N,by=group];setnames(I,"N","I");
            totS=pop[IndexCase_0==1,.N];
            totI=pop[IndexCase_0==0,.N];
            
            setkey(SIRts,group);setkey(S,group);setkey(I,group);
            SIRts=SIRts[S][I];
            SIRts[,totS:=totS][,totI:=totI];
            
            nextT = 0;
        };#Endif SI
        
        if(Epi.type=="SIR"){
            SIRts = data.table(
                group=unique(pop[,group]),
                time=0);
            #*************************************
            
            S=pop[IndexCase_0==1,.N,by=group];setnames(S,"N","S");
            I=pop[IndexCase_0==0,.N,by=group];setnames(I,"N","I");
            #*NEW*#
            R=0;
            #*NEW*#
            
            totS=pop[IndexCase_0==1,.N];
            totI=pop[IndexCase_0==0,.N];
            #*NEW*#
            totR=0;
            #*NEW*#
            
            setkey(SIRts,group);setkey(S,group);setkey(I,group);
            SIRts=SIRts[S][I];
            SIRts[,R:=R][,totS:=totS][,totI:=totI][,totR:=totR];
            
            nextT = 0;
        };#Endif SIR
        
        
    };#End else
    
    
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    sumfs2=pop[I==1&frailties=="lognormal",sum(infectivity_ln),by=group];
    
    sumfs=data.table(group=unique(pop[,group]));
    setkey(sumfs,group);setkey(sumfs2,group);
    sumfs=sumfs2[sumfs];
    sumfs[is.na(V1),V1:=0];
    
    #pop[,V1:=NULL];
    setkey(pop,group);setkey(sumfs,group);
    pop = pop[sumfs];
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #while loop for SI model should stop once everyone is infected as well as if no-one is infected.
    #For SIS, while loop stops when no-one is infected (**could go on a long time!!**).
    #Current while loop should be able to accommodate SIS model, but need another stopping argument.
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if(Epi.type=="SI"){
        #SI loop#######################################################################################################
        while(SIRts[,last(totI)>0&last(totI)<nind]){#Start of epidemic simulation loop#
            
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PairwiseBeta%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
            pop[I==0&frailties=="lognormal",eRates:=susceptibility_ln*PairwiseBeta*V1];#V1 is the summed infectivity_ln.
            pop[I==1,eRates:=0];
            
            
            ###Taking the whole population together###
            nextT=nextT+rexp(1, rate = pop[,sum(eRates)]);
            IDnextT=sample(pop[,id],size=1,prob=pop[,eRates]/pop[,sum(eRates)]);
            groupnextT=pop[id==IDnextT,group];
            
            SIRtstemp=data.table(group=groupnextT,time=nextT,
                                 S=SIRts[group==groupnextT,last(S)],
                                 I=SIRts[group==groupnextT,last(I)]);
            
            SIRtstemp[pop[IDnextT==id,I!=1],
                      S:=SIRts[group==groupnextT,last(S)] - 1][pop[IDnextT==id,I!=1],
                                                               I:=SIRts[group==groupnextT,last(I)] + 1][pop[IDnextT==id,I!=1],
                                                                                                        totS:=SIRts[,last(totS)] - 1][pop[IDnextT==id,I!=1],
                                                                                                                                      totI:=SIRts[,last(totI)] + 1];
            
            SIRts=rbind(SIRts,SIRtstemp);
            
            pop[IDnextT==id&I!=1,
                Tinf:=nextT][IDnextT==id&I!=1,
                             I:=1];
            ###
            
            #RECALCULATE V1#
            sumfs2=pop[I==1&frailties=="lognormal",sum(infectivity_ln),by=group];
            
            sumfs=data.table(group=unique(pop[,group]));
            setkey(sumfs,group);setkey(sumfs2,group);
            sumfs=sumfs2[sumfs];
            sumfs[is.na(V1),V1:=0];
            
            pop[,V1:=NULL];
            setkey(pop,group);setkey(sumfs,group);
            pop = pop[sumfs];
            
        };#End of epidemic simulation loop#
        #SI loop above##################################################################################################
    }#End of SI if loop.
    
    
    if(Epi.type=="SIR"){
        #SIR loop#######################################################################################################
        while(SIRts[,last(totI)>0]){#&last(totI)<nind])#Start of epidemic simulation loop#I don't get why you would stop if all individuals are infected. They still need time to recover.
            
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PairwiseBeta%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
            pop[R==1,eRates:=0];
            pop[R!=1&I==0&frailties=="lognormal",eRates:=susceptibility_ln*PairwiseBeta*V1];#V1 is the summed infectivity_ln
            pop[R!=1&I==1&frailties=="lognormal",eRates:=recoverability_ln*RRgamma];
            
            ###Taking the whole population together###
            nextT=nextT+rexp(1, rate = pop[,sum(eRates)]);
            IDnextT=sample(pop[,id],size=1,prob=pop[,eRates]/pop[,sum(eRates)]);
            groupnextT=pop[id==IDnextT,group];
            
            SIRtstemp=data.table(group=groupnextT,time=nextT,
                                 S=SIRts[group==groupnextT,last(S)],
                                 I=SIRts[group==groupnextT,last(I)],
                                 R=SIRts[group==groupnextT,last(R)]);
            
            SIRtstemp[pop[IDnextT==id,I!=1], S:=SIRts[group==groupnextT,last(S)] - 1]
            SIRtstemp[pop[IDnextT==id,I!=1], I:=SIRts[group==groupnextT,last(I)] + 1]
            SIRtstemp[pop[IDnextT==id,I!=1], R:=SIRts[group==groupnextT,last(R)]]
            SIRtstemp[pop[IDnextT==id,I!=1], totS:=SIRts[,last(totS)] - 1]
            SIRtstemp[pop[IDnextT==id,I!=1], totI:=SIRts[,last(totI)] + 1]
            SIRtstemp[pop[IDnextT==id,I!=1], totR:=SIRts[,last(totR)]]
            SIRtstemp[pop[IDnextT==id,I==1], I:=SIRts[group==groupnextT,last(I)] - 1]
            SIRtstemp[pop[IDnextT==id,I==1], R:=SIRts[group==groupnextT,last(R)] + 1]
            SIRtstemp[pop[IDnextT==id,I==1], S:=SIRts[group==groupnextT,last(S)]]
            SIRtstemp[pop[IDnextT==id,I==1], totI:=SIRts[,last(totI)] - 1]
            SIRtstemp[pop[IDnextT==id,I==1], totR:=SIRts[,last(totR)] + 1]
            SIRtstemp[pop[IDnextT==id,I==1], totS:=SIRts[,last(totS)]];
            
            SIRts=rbind(SIRts,SIRtstemp);
            
            pop[IDnextT==id&I!=1,
                Tinf:=nextT][IDnextT==id&I==1,
                             Trec:=nextT][IDnextT==id&I==1,
                                          c("I","R"):=list(0,1)][IDnextT==id&I!=1&R!=1,
                                                                 I:=1];
            ###
            
            #RECALCULATE V1#
            sumfs2=pop[I==1&frailties=="lognormal",sum(infectivity_ln),by=group];
            
            sumfs=data.table(group=unique(pop[,group]));
            setkey(sumfs,group);setkey(sumfs2,group);
            sumfs=sumfs2[sumfs];
            sumfs[is.na(V1),V1:=0];
            
            pop[,V1:=NULL];
            setkey(pop,group);setkey(sumfs,group);
            pop = pop[sumfs];
            
        };#End of epidemic simulation loop#
        #SIR loop above##################################################################################################
    }#End of SIR if loop.
    
    
    ########################################################################################
    ########################################################################################
    ########################################################################################
    
    output = list();
    output$pop = pop;
    output$SIRts = SIRts;
    
    return(output);
    
}
