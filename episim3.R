EpiSim3 <- function(
    # No. of individuals
    nind,
    # For explicitly specifying group and index case status of each individual.
    # Requires a data.frame with 'nind' rows and columns for group (character),
    # donor/recipient ("D" or "R") and individual id
    group.and.index = NULL,
    # ngroups now has a default (17 July 2020)
    ngroups = 1,
    # Here I use a single scalar value, but this can also be a vector of 'fixed
    # effect' Beta values the same length as nind and assuming individuals are
    # in the same order as in the pedigree. Could be useful if some variable
    # e.g. body weight affects epidemiological trait values
    PairwiseBeta,
    RRgamma,
    # this is the traitdata previously generated
    pop = NULL,
    offspring.only = FALSE,
    # The genomic relationship matrix excluding parents, created above
    GRM,
    GRM.from.pedigree = FALSE,
    ntraits = 3,
    tmeans = rep(0, 3),
    tVA, trhoG, tVE, trhoE,
    traitnames = c("susceptibility", "infectivity", "recoverability"),
    trait.transform = c("lognormal", "invcloglog"),
    # Currently individual frailties are assumed lognormal, but this could be
    # made more flexible, e.g. allowing gamma-distributed epidemiological trait
    # values
    frailties = c("lognormal", "gamma"),
    # Whether this is an SI or SIR epidemic. More options could be added
    Epi.type = c("SI", "SIS", "SIR"),
    DEBUG = FALSE
) {
    gr.size <- nind / ngroups
    
    # Simulate trait values, call existing function
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
        if (offspring.only == TRUE) {
            pop <- pop[dso == "offspring"]
        }
    }
    
    pop[, PairwiseBeta := PairwiseBeta]
    
    
    # Assign groups and index cases, call existing function
    # --------------------------------------------------------------------------
    if (is.null(group.and.index)) {
        pop <- group.assign(pop = pop, ngroups = ngroups)
        pop[, frailties := frailties]
        pop[, Epi.type := Epi.type]
    } else {
        # Custom assignment of groups and index cases
        # setkey(pop,id)
        # setkey(group.and.index, individual_ForPedigree)
        # Don't want to reorder the rows!
        pop[, group := group.and.index[, Box_recoded]]
        pop[id %in% group.and.index[DR == "D", individual_ForPedigree], IndexCase_0 := 0]
        pop[is.na(IndexCase_0), IndexCase_0 := 1]
        pop[IndexCase_0 == 0, Tinf := 0]
    }
    # Will need to generalize this later
    
    # Generate epidemics
    #---------------------------------------------------------------------------
    
    # Store everything in a table indicating group and epidemic
    # Basically recreate `pop` and `SIRts` but with all runs in a single table
    # I = infected, R = recovered.
    pop[Epi.type == "SIR", c("I", "R") := list(1 - IndexCase_0, 0)]
    pop[Epi.type != "SIR", I := 1 - IndexCase_0]
    
    # Each individual has to have a contact matrix
    # Do this later
    # --------------------------------------------------------------------------
    
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
            nextT <- 0
        } else if (Epi.type == "SI") {
            SIRts <- data.table(
                group = unique(pop[, group]),
                time = 0,
                S = gr.size - 1,
                I = 1,
                totS = nind - ngroups,
                totI = ngroups)
            nextT <- 0
        }
    } else {
        # if group.and.index is not null
        
        if (Epi.type ==  "SI") {
            SIRts <- data.table(
                group = unique(pop[, group]),
                time = 0)
            
            S <- pop[IndexCase_0 == 1, .N, by = group]
            setnames(S, "N", "S")
            
            I <- pop[IndexCase_0 == 0, .N, by = group]
            setnames(I, "N", "I")
            
            totS <- pop[IndexCase_0 == 1, .N]
            totI <- pop[IndexCase_0 == 0, .N]
            
            
            setkey(SIRts, group)
            setkey(S, group)
            setkey(I, group)
            
            SIRts <- SIRts[S][I]
            SIRts[, totS := totS][, totI := totI]
            
            #,
            # S <- gr.size - 1,
            # I <- 1,
            # totS <- nind - ngroups,
            # totI <- ngroups);
            #
            
            # Deal with SIR later
            nextT <- 0
        }
    }
    
    # For each group in `pop`, create value V1 which is the sum of the log
    # infectivities. I'm guessing this gives an overall infectivity for the
    # group?
    sumfs <- pop[I == 1 & frailties == "lognormal", sum(infectivity_ln), by = group]
    setkey(pop, group)
    setkey(sumfs, group)
    pop <- pop[sumfs]
    
    
    
    # While loop for SI model
    # 
    # Should stop once everyone is infected as well as if no-one is infected.
    # For SIS, while loop stops when no-one is infected (**could go on a long
    # time!!**). Current while loop should be able to accommodate SIS model, but
    # need another stopping argument.
    # --------------------------------------------------------------------------
    
    
    # SI loop
    if (Epi.type == "SI") {
        # Start of epidemic simulation loop
        while (SIRts[, last(totI) > 0 & last(totI) < nind]) { 
            # PairwiseBeta
            # V1 is the summed infectivity_ln.
            pop[I == 0 & frailties == "lognormal", eRates := susceptibility_ln * PairwiseBeta * V1]
            pop[I == 1, eRates := 0]
            
            # fwrite(pop[, .(susceptibility_ln, PairwiseBeta, V1, eRates)], "poptest.csv")
            
            
            # Taking the whole population together
            nextT <- nextT + rexp(1, rate = pop[, sum(eRates)])
            IDnextT <- sample(pop[, id], size = 1, prob = pop[, eRates] / pop[, sum(eRates)])
            groupnextT <- pop[id == IDnextT, group]
            
            SIRts_tmp <- data.table(
                group = groupnextT,
                time = nextT, 
                S = SIRts[group == groupnextT, last(S)], 
                I = SIRts[group == groupnextT, last(I)]
            )
            
            SIRts_tmp[pop[IDnextT == id, I != 1], S := SIRts[group == groupnextT, last(S)] - 1]
            SIRts_tmp[pop[IDnextT == id, I != 1], I := SIRts[group == groupnextT, last(I)] + 1]
            SIRts_tmp[pop[IDnextT == id, I != 1], totS := SIRts[, last(totS)] - 1]
            SIRts_tmp[pop[IDnextT == id, I != 1], totI := SIRts[, last(totI)] + 1]
            
            SIRts <- rbind(SIRts, SIRts_tmp)
            
            pop[IDnextT == id & I != 1,
                Tinf := nextT][IDnextT == id & I != 1, I := 1]
            
            
            # print(head(pop))
            # print(SIRts)
            
            # Recalculate V1
            sumfs2 <- pop[I == 1 & frailties == "lognormal", sum(infectivity_ln), by = group]
            
            
            sumfs <- data.table(group = unique(pop[, group]))
            
            setkey(sumfs, group)
            setkey(sumfs2, group)
            
            sumfs <- sumfs2[sumfs]
            
            sumfs[is.na(V1), V1 := 0]
            
            print(list(sumfs, nrow(sumfs)))
            
            pop[,V1:=NULL]
            setkey(pop,group)
            setkey(sumfs,group)
            pop <- pop[sumfs]
            
            # print(head(pop))
            # print(SIRts)
            
        } # End of epidemic simulation loop
        # SI loop above
        cat("finished the SI model\n")
    } # End of SI loop
    
    
    # SIR loop
    if (Epi.type == "SIR") {
        
        # Start of epidemic simulation loop
        # I don't get why you would stop if all individuals are infected. They
        # still need time to recover.
        
        #&last(totI)<nind])
        while (SIRts[, last(totI) > 0]) {
            
            # PairwiseBeta
            pop[R == 1, eRates := 0]
            # V1 is the summed infectivity_ln
            pop[R != 1 & I == 0 & frailties == "lognormal", eRates := susceptibility_ln * PairwiseBeta * V1]
            pop[R != 1 & I == 1 & frailties == "lognormal", eRates := recoverability_ln * RRgamma]
            
            
            # Taking the whole population together
            nextT <- nextT + rexp(1, rate = pop[, sum(eRates)])
            IDnextT <- sample(pop[, id],  size = 1,
                              prob = pop[, eRates] / pop[, sum(eRates)])
            groupnextT <- pop[id == IDnextT, group]
            
            SIRts_tmp <- data.table(
                group = groupnextT,
                time = nextT,
                S = SIRts[group == groupnextT, last(S)],
                I = SIRts[group == groupnextT, last(I)],
                R = SIRts[group == groupnextT, last(R)]
            )
            
            SIRts_tmp[pop[IDnextT == id, I != 1], S := SIRts[group == groupnextT, last(S)] - 1]
            SIRts_tmp[pop[IDnextT == id, I != 1], I := SIRts[group == groupnextT, last(I)] + 1]
            SIRts_tmp[pop[IDnextT == id, I != 1], R := SIRts[group == groupnextT, last(R)]]
            SIRts_tmp[pop[IDnextT == id, I != 1], totS := SIRts[, last(totS)] - 1]
            SIRts_tmp[pop[IDnextT == id, I != 1], totI := SIRts[, last(totI)] + 1]
            SIRts_tmp[pop[IDnextT == id, I != 1], totR := SIRts[, last(totR)]]
            SIRts_tmp[pop[IDnextT == id, I == 1], I := SIRts[group == groupnextT, last(I)] - 1]
            SIRts_tmp[pop[IDnextT == id, I == 1], R := SIRts[group == groupnextT, last(R)] + 1]
            SIRts_tmp[pop[IDnextT == id, I == 1], S := SIRts[group == groupnextT, last(S)]]
            SIRts_tmp[pop[IDnextT == id, I == 1], totI := SIRts[, last(totI)] - 1]
            SIRts_tmp[pop[IDnextT == id, I == 1], totR := SIRts[, last(totR)] + 1]
            SIRts_tmp[pop[IDnextT == id, I == 1], totS := SIRts[, last(totS)]]
            
            SIRts <- rbind(SIRts, SIRts_tmp)
            
            pop[IDnextT == id & I != 1, Tinf := nextT]
            pop[IDnextT == id & I == 1, Trec := nextT]
            pop[IDnextT == id & I == 1, c("I", "R") := list(0, 1)]
            pop[IDnextT == id & I != 1 & R != 1, I := 1]
            
            
            if (DEBUG) {
                print(head(pop))
                print(SIRts)
            }
            
            # Recalculate V1
            sumfs2 <- pop[I == 1 & frailties == "lognormal", sum(infectivity_ln), by = group]
            
            
            sumfs <- data.table(group = unique(pop[, group]))
            setkey(sumfs, group)
            setkey(sumfs2, group)
            sumfs <- sumfs2[sumfs]
            sumfs[is.na(V1), V1 := 0]
            if (DEBUG) {
                print(list(sumfs, nrow(sumfs)))
            }
            
            pop[, V1 := NULL]
            setkey(pop, group)
            setkey(sumfs, group)
            pop <- pop[sumfs]
            
            
            if (DEBUG) {
                print(head(pop))
                print(SIRts)
            }
        } # End of epidemic simulation loop
        # SIR loop above
    } # End of SIR if loop.
    
    
    # return
    list(pop = pop, SIRts = SIRts)
}
