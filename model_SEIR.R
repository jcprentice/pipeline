get_infectives <- function(X) {
    return(X[, sum(status %in% c("E", "I"))])
}


simulate_SEIR <- function(pars, traitdata) {
    cat("Simulating an SEIR epidemic ...\n")

    with(pars, { # handy for expanding parameters
        pop <- traitdata

        if (offspring_only)
            pop <- pop[dso == "offspring"]

        # nind <- nrow(pop)
        # group_size <- ceiling(nind / ngroups)


        # Assign groups and index cases
        pop <- group_assign(pop, ngroups)
        setkey(pop, group)

        # Generate epidemics

        # Store everything in a table indicating group and epidemic
        # Basically recreate `pop` and `SIRts` but with all runs in a single table
        # S = Susceptible, E = Exposed, I = Infected, R = Recovered.
        pop[, status := factor("S", levels = c("S", "E", "I", "R"))]
        pop[, c("Tinf", "Trec") := as.numeric(NA)]
        pop[index_case == 1, c("status", "Tinf") := list("E", 0.0)]

        # ----------------------------------------------------------------------
        # Start epidemic simulation loop
        epi_time <- 0
        while (get_infectives(pop) > 0) {
            if (DEBUG) {
                cat("time =", epi_time)
            }

            # Calculate infection rates in each group
            # this is the sum of the log infectivitities
            # if ("lognormal" %in% frailties)
            pop[, group_inf := sum(infectivity_ln * (status == "I")), by = group]


            # if S, infection at rate beta SI
            pop[status == "S", event_rate := susceptibility_ln * r_beta * group_inf]
            # if E, advance to I at rate sigma
            pop[status == "E", event_rate := r_sigma]
            # if I, recover at rate gamma
            pop[status == "I", event_rate := recoverability_ln * r_gamma]
            # if individual is R, set event rate to 0 (it's done)
            pop[status == "R", event_rate := 0]


            # generate random timestep
            dt <- rexp(1, rate = pop[, sum(event_rate)])
            epi_time <- epi_time + dt
            # randomly select individual
            IDnextT <- sample(pop[, id],
                              size = 1,
                              prob = pop[, event_rate] / pop[, sum(event_rate)])
            # select group
            groupnextT <- pop[id == IDnextT, group]


            status = pop[IDnextT == id, status]

            # apparently this doesn't like factors?
            if (status == "S") {
                pop[id == IDnextT, Tinf := epi_time]
                pop[id == IDnextT, status := "I"]
            } else if (status == "E") {
                pop[id == IDnextT, Tinf := epi_time]
                pop[id == IDnextT, status := "I"]
            } else if (status == "I") {
                pop[id == IDnextT, Trec := epi_time]
                pop[id == IDnextT, status := "R"]
            } else {
                disp("unexpected event!")
            }


            if (DEBUG) {
                cat("\npop =\n")
                print(pop[, .(group, index_case, status, sum_inf, event_rate)])
                cat("\nSIRts =\n")
                print(SIRts)
                cat("\nsumfs =\n")
                print(sumfs)
            }
        }


        # tidy up pop
        pop[, c("group_inf", "event_rate") := NULL]
        setkey(pop, id)

        # browser()

        # return
        return(list(pop = pop, SIRts = NULL))

    }) # with()
}
