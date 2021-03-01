# Function for discretising time ready for GLMM
#-------------------------------------------------------------------------------

discretize.time <- function(
    # The population file
    pop,
    # No. of individuals
    nind,
    # No. of groups
    ngroups,
    # The duration of each discrete time interval. If not specified, a new
    # interval starts at each transition event
    interval.duration = "event", 
    # If TRUE, removes all time points after infection for each individual (none
    # removed for individuals that didn't become infected). Needed for GLMM
    # analysis. This should be FALSE for Chris's data.
    truncate = TRUE,
    # Epidemic type
    Epi.type = c("SI", "SIS", "SIR")
) {
    
    
    if (interval.duration == "event") {
        setkey(pop, NULL)
        setkey(pop, group, Tinf)
        pop[2:nind, Tinflast := pop[1:(nind - 1), Tinf]]
        pop[Tinf == 0, Tinflast := NA]
        pop[, int.duration := Tinf - Tinflast]
        breaks <- seq(0, (nind / ngroups) - 1)
        # FIXME: breaks at this line
        pop[, Tinf2 := breaks, by = group]
        
        if (Epi.type == "SIR") {
            pop[, Trec2 := ceiling(Trec / int.duration)]
        }
    } else {
        cat("not event")
        pop[, int.duration := interval.duration]
        pop[, Tinf2 := ceiling(Tinf / int.duration)]
        if (Epi.type == "SIR") {
            pop[, Trec2 := ceiling(Trec / int.duration)]
        }
    }
    # TODO: Still haven't sorted Trec2
    
    if (Epi.type == "SIR") {
        pop[is.na(Tinf), censored := 0]
        pop[Tinf <= max(na.omit(Trec2)), censored := 1]
    } else {
        pop[is.na(Tinf), censored := 0]
        pop[Tinf <= max(na.omit(Tinf2)), censored := 1]
    }
    
    if (Epi.type == "SIR") {
        dt1 <- data.table(expand.grid(pop[, id], seq(1, pop[, max(na.omit(Trec2))])))
    } else {
        dt1 <- data.table(expand.grid(pop[, id], seq(1, pop[, max(na.omit(Tinf2))])))
    }
    
    setnames(dt1, c("Var1", "Var2"), c("id", "Tinfind"))
    # FIXME: why set this to factor?
    #dt1[, id := as.factor(id)]
    
    setkey(dt1, id);
    setkey(pop, id);
    
    if (Epi.type == "SIR") {
        if (interval.duration == "event") {
            popvv <- pop[, .(id, group, Tinf2, Trec2, Tinf, Trec, censored, IndexCase_0)][dt1]
            pop[, Tinfind := Tinf2]
            pop[, Trecind := Trec2]
            setkey(pop, group, Tinfind)
            setkey(popvv, group, Tinfind)
            popvv <- pop[, .(group, Tinfind, int.duration)][popvv]
        } else {
            popvv <- pop[, .(id, group, Tinf2, Trec2, Tinf, Trec, censored, IndexCase_0, int.duration)][dt1]
        }
    } else {
        if (interval.duration == "event") {
            popvv <- pop[, .(id, group, Tinf2, Tinf, censored, IndexCase_0)][dt1]
            pop[, Tinfind := Tinf2]
            setkey(pop, group, Tinfind)
            setkey(popvv, group, Tinfind)
            popvv <- pop[, .(group, Tinfind, int.duration)][popvv]
        } else {
            popvv <- pop[, .(id, group, Tinf2, Tinf, censored, IndexCase_0, int.duration)][dt1]
        }
    }
    
    # if (Epi.type == "SIR") {
    #     # FIXME: this breaks if dt1$id is a factor
    #     popvv <- pop[, .(id, group, Tinf2, Trec2, Tinf, Trec, censored, IndexCase_0)][dt1]
    # }
    # 
    # if (Epi.type != "SIR" & interval.duration == "event") {
    #     popvv <- pop[, .(id, group, Tinf2, Tinf, censored, IndexCase_0)][dt1]
    #     pop[, Tinfind := Tinf2]
    #     setkey(pop, group, Tinfind)
    #     setkey(popvv, group, Tinfind)
    #     popvv <- pop[, .(group, Tinfind, int.duration)][popvv]
    # }
    # 
    # if (Epi.type != "SIR" & interval.duration != "event") {
    #     popvv <- pop[, .(id, group, Tinf2, Tinf, censored, IndexCase_0, int.duration)][dt1]
    # }
    
    popvv[Tinfind == Tinf2, pit := 1][is.na(pit), pit := 0]
    popvv[Tinfind == Tinf2 & censored == 1, pit := 1][is.na(pit), pit := 0]
    
    if (Epi.type == "SIR") {
        # (Already infected by current time interval - already dead by cti).
        # Died or became infected during current time interval don't count.
        I_t <- popvv[, sum(na.omit(Tinf2) < Tinfind) - sum(na.omit(Trec2) < Tinfind), by = c("group", "Tinfind")]
        
    } else {
        I_t <- popvv[, sum(na.omit(Tinf2) < Tinfind), by = c("group", "Tinfind")]
    }
    
    
    if (Epi.type == "SIR") {
        popvv[Tinf2 < Tinfind & Trec2 >= Tinfind, infected := 1]
    } else {
        popvv[Tinf2 < Tinfind, infected := 1]
    }
    popvv[is.na(infected), infected := 0]
    
    setnames(I_t, "V1", "I_t")
    setkey(I_t, group, Tinfind)
    setkey(popvv, group, Tinfind)
    popvv = I_t[popvv]
    
    
    # Need NAs added as these are always susceptible
    S_t <- popvv[, sum(na.omit(Tinf2) >= Tinfind) + sum(is.na(Tinf2)), by = c("group", "Tinfind")]
    N_t <- data.table(
        Tinfind = S_t[, Tinfind],
        group = S_t[, group],
        N_t = S_t[, V1] + I_t[, I_t],
        S_t = S_t[, V1]
    )
    
    setkey(N_t, group, Tinfind)
    setkey(popvv, group, Tinfind)
    popvv <- N_t[popvv]
    
    # This is the Biemens offset. Osvaldo doesn't divide by N in his eRates
    # calculation
    # FIXME: int.duration not found?
    popvv[, logoffset_RelFreq := log((I_t / N_t) * int.duration)]
    # Offset to match Osvaldo's simulation.
    popvv[, logoffset_deltat := log(I_t * int.duration)]
    
    # Removing the index cases.
    popvv <- popvv[id != popvv[Tinf == 0, unique(id)]]
    
    popvv[, group := as.factor(group)]
    popvv[, animal := id]
    popvv[, id2 := id]
    popvv[, pit_F := as.factor(pit)]
    
    
    # Removing all time points after infection for each individual (not removing
    # any for those that didn't become infected)
    if (truncate) {
        popvv <- popvv[Tinfind <= Tinf2 | is.na(Tinf2)]
    }
    
    # filtering individuals still S with no more I's left.
    popvv <- popvv[logoffset_deltat > -Inf]
    
    # return
    list(pop = pop, discpop = popvv)
}
