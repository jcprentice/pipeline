model_SIR <- function(pars, traitdata) {
    cat("Simulating an SIR epidemic ...\n")

    with(pars, { # handy for expanding parameters
        # ngroups <- pars$ngroups # 1
        # IFbeta <- pars$PairwiseBeta
        # RRgamma <- pars$RRgamma
        # pop <- traitdata
        # offspring_only = pars$offspring_only
        # ntraits <- pars$ntraits # 3
        # tmeans <- pars$tmeans # rep(0, 3),
        # tVA <- pars$tVA
        # trhoG <- pars$trhoG
        # tVE <- pars$tVE
        # trhoE <- pars$trhoE
        # traitnames <- pars$traitnames
        # trait.transform <- pars$trait_transform # c("lognormal", "invcloglog"),
        # # Currently individual frailties are assumed lognormal, but this could be
        # # made more flexible, e.g. allowing gamma-distributed epidemiological trait
        # # values
        # frailties <- pars$frailties # c("lognormal", "gamma"),
        # epi_type <- pars$epi_type
        # DEBUG <- pars$DEBUG


        pop <- traitdata

        if (offspring_only) {
            pop <- pop[dso == "offspring"]
        }

        nind <- nrow(pop)
        group_size <- ceiling(nind / ngroups)


        # Assign groups and index cases
        # --------------------------------------------------------------------------
        pop <- group_assign(pop, ngroups)
        setkey(pop, group)

        # Generate epidemics
        #---------------------------------------------------------------------------

        # Store everything in a table indicating group and epidemic
        # Basically recreate `pop` and `SIRts` but with all runs in a single table
        # S = Susceptible, E = Exposed, I = Infected, R = Recovered.
        pop[, status := "S"]
        pop[index_case == 1, c("status", "Tinf") := list("I", 0)]

        SIRts <- data.table(
            group = 1:ngroups,
            time = 0,
            S = group_size - 1,
            I = 1,
            R = 0,
            totS = nind - ngroups, # because 1 infected per group
            totI = ngroups,
            totR = 0
        )

        # Start of epidemic simulation loop

        # For each group in `pop`, create value `sum_inf` which is the sum of
        # the log infectivities. I'm guessing this gives an overall infectivity
        # for the group?

        # if ("lognormal" %in% frailties)
        sumfs <- pop[status == "I", sum(infectivity_ln), by = group]
        setnames(sumfs, "V1", "sum_inf")
        setkey(sumfs, group)
        pop <- pop[sumfs]


        epi_time <- 0
        while (SIRts[, last(totI) > 0]) {
            if (DEBUG) {
                cat("time =", epi_time)
            }

            # if S, infection at rate beta SI
            pop[status == "S", eRates := susceptibility_ln * r_beta * sum_inf]
            # if I, recover at rate gamma I
            pop[status == "I", eRates := recoverability_ln * r_gamma]
            # if individual is R, set event rate to 0 (it's done)
            pop[status == "R", eRates := 0]


            # generate random timestep
            dt <- rexp(1, rate = pop[, sum(eRates)])
            epi_time <- epi_time + dt
            # randomly select individual
            IDnextT <- sample(pop[, id],
                              size = 1,
                              prob = pop[, eRates] / pop[, sum(eRates)])
            # select group
            groupnextT <- pop[id == IDnextT, group]

            SIRts_tmp <- data.table(
                group = groupnextT,
                time = epi_time,
                S = SIRts[group == groupnextT, last(S)],
                I = SIRts[group == groupnextT, last(I)],
                R = SIRts[group == groupnextT, last(R)]
            )

            status = pop[IDnextT == id, status]
            if (status == "S") {
                SIRts_tmp[, S := SIRts[group == groupnextT, last(S)] - 1]
                SIRts_tmp[, I := SIRts[group == groupnextT, last(I)] + 1]
                SIRts_tmp[, R := SIRts[group == groupnextT, last(R)]]
                SIRts_tmp[, totS := SIRts[, last(totS)] - 1]
                SIRts_tmp[, totI := SIRts[, last(totI)] + 1]
                SIRts_tmp[, totR := SIRts[, last(totR)]]
            } else if (status == "I") {
                SIRts_tmp[, I := SIRts[group == groupnextT, last(I)] - 1]
                SIRts_tmp[, R := SIRts[group == groupnextT, last(R)] + 1]
                SIRts_tmp[, S := SIRts[group == groupnextT, last(S)]]
                SIRts_tmp[, totI := SIRts[, last(totI)] - 1]
                SIRts_tmp[, totR := SIRts[, last(totR)] + 1]
                SIRts_tmp[, totS := SIRts[, last(totS)]]
            } else {
                print("Individual is R, no event should have happened!")
            }

            SIRts <- rbind(SIRts, SIRts_tmp)

            if (status == "S") {
                pop[id == IDnextT, Tinf := epi_time]
                pop[id == IDnextT, status := "I"]
            } else {
                pop[id == IDnextT, Trec := epi_time]
                pop[id == IDnextT, status := "R"]
            }

            # Recalculate sum_inf
            sumfs2 <- pop[status == "I", sum(infectivity_ln), by = group]
            setkey(sumfs2, group)

            sumfs[, sum_inf := NULL]
            sumfs <- sumfs2[sumfs]
            setnames(sumfs, "V1", "sum_inf")
            sumfs[is.na(sum_inf), sum_inf := 0]

            pop[, sum_inf := NULL]
            pop <- pop[sumfs]


            if (DEBUG) {
                cat("\npop =\n")
                print(pop[, .(group, index_case, status, sum_inf, eRates)])
                cat("\nSIRts =\n")
                print(SIRts)
                cat("\nsumfs =\n")
                print(sumfs)
            }
        }


        # tidy up pop
        pop[, c("sum_inf", "eRates") := NULL]
        setkey(pop, id)

        # browser()

        # return
        return(list(pop = pop, SIRts = SIRts))

    }) # with()
}
