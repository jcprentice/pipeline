get_parameters <- function(epi_type = "SIR") {
    cat("running get_parameters() ...\n")
    # No. of sires, dams per sire, and offspring per pair, and groups
    nsires <- 3
    dpsire <- 2
    oppair <- 2
    ngroups <- 2

    # Derived population numbers
    ndams <- nsires * dpsire
    nparents <- nsires + nsires * dpsire
    nprogeny <- nsires * dpsire * oppair
    nind <- nprogeny
    ntotal <- nparents + nprogeny

    # TODO: I think it should be this, but we'll keep it as is for the moment
    # nind = nparents + nprogeny

    # TODO: what if this isn't an integer?
    group_size = nprogeny / ngroups

    # Traits
    if (epi_type == "SEIR") {
        traitnames <- c("susceptibility", "latency", "infectivity", "recoverability")
    } else
        traitnames <- c("susceptibility", "infectivity", "recoverability")
    }


    # No. of traits to simulate.
    ntraits <- length(traitnames)

    # Trait means
    tmeans <- rep(0, ntraits)

    # Trait additive genetic variance
    if (epi_type == "SEIR") {
        tVA <- c(0.5, 0.1, 0.1, 0.1)
    } else {
        tVA <- c(0.5, 0.1, 0.1)
    }

    # Trait non-additive (environmental) variance
    # with no association to the pedigree
    tVE <- rep(0.3, ntraits)
    # Trait non-additive correlation matrix
    trhoE <- diag(ntraits)
    # Trait additive genetic correlation matrix.
    # In this case correlations are zero.
    trhoG <- diag(ntraits)

    # also used:
    trait_transform = c("lognormal") # c("lognormal", "invcloglog")
    frailties = "lognormal"


    # Infection and recoverability, targeting R0 = 2.5?
    # epi_type <- "SIR"
    r_beta <- 0.5
    r_gamma <- 0.2
    r_sigma <- 1.0


    # MCMC settings
    nsample <- 5000
    burnin <- 500

    # Additional options
    trim_parents <- FALSE
    offspring_only <- TRUE
    interval_duration <- "event"
    truncate <- TRUE
    DEBUG <- FALSE


    # This will be a list for passing to functions
    pars <- list(nsires = nsires, ndams = ndams, nparents = nparents, nprogeny = nprogeny,
                 nind = nind, ntotal = ntotal, dpsire = dpsire, oppair = oppair,
                 ngroups = ngroups, traitnames = traitnames, ntraits = ntraits,
                 tmeans = tmeans, tVA = tVA, tVE = tVE, trhoE = trhoE,
                 trhoG = trhoG, trait_transform = trait_transform, frailties = frailties,
                 epi_type = epi_type, r_beta = r_beta, r_gamma = r_gamma, r_sigma = r_sigma,
                 nsample = nsample, burnin = burnin,
                 trim_parents = trim_parents, offspring_only = offspring_only,
                 interval_duration = interval_duration, truncate = truncate,
                 DEBUG = DEBUG)

    # Quickly summarise what we've got
    cat(paste0(
        "Generating data with:\n",
        nsires, " sires, ", ndams, " dams, ", nprogeny, " progeny, (", ntotal, " individuals), ", ngroups, " groups\n",
        "R0 = ", r_beta / r_gamma, "?\n"
    ))

    return(pars)
}
