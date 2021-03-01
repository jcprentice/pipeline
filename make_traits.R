make_traits <- function(pars, pedigree, GRM) {
    # require(data.table)
    # require(MASS)
    # require(corpcor)

    cat("Making trait values ...\n")

    # extract parameters
    ntotal <- pars$ntotal
    traitnames <- pars$traitnames
    ntraits <- pars$ntraits
    tmeans <- pars$tmeans
    tVA <- pars$tVA
    tVE <- pars$tVE
    trhoE <- pars$trhoE
    trhoG <- pars$trhoG
    DEBUG <- pars$DEBUG
    trait_transform <- pars$trait_transform


    # TODO: figure out what this is doing exactly
    tab1 <- data.table(
        expand.grid(trait1 = 1:ntraits, trait2 = 1:ntraits),
        corrval = c(trhoG) # as.numeric(combn(trhoG, 1))
    )
    tab1[, tvar1 := tVA[trait1]]
    tab1[, tvar2 := tVA[trait2]]
    tab1[, traitcov := corrval * sqrt(tvar1) * sqrt(tvar2)]

    # Environmental covariance
    tab1[, Ecorrval := c(trhoE)] # as.numeric(combn(trhoE, 1))
    tab1[, Evar1 := tVE[trait1]]
    tab1[, Evar2 := tVE[trait2]]
    tab1[, traitEcov := Ecorrval * sqrt(Evar1) * sqrt(Evar2)]

    # Creating the global multitrait individual pairwise covariance matrix

    # Create a similar pairwise matrix for the GRM.
    tab2 <- as.data.table(reshape2::melt(GRM))
    setnames(tab2, c("Var1", "Var2", "value"), c("ind1", "ind2", "pshared"))

    # GRM * covariance AND GRM * variance.
    tab3 <- as.data.table(tab2[, pshared] %*% t(tab1[, traitcov]))

    tab4 <- cbind(tab2, tab3)
    tab4 <- tab4[, pshared := NULL]

    GRM_Nt <- Matrix(0, # why NA?
                     nrow = ntraits * ntotal,
                     ncol = ntraits * ntotal)
    k <- 0
    for (i in 1:ntraits) { for (j in 1:ntraits) {
        # i = columns, j = Rows
        k <- k + 1

        # TODO: rewrite it to make it less evil

        # this is a horrible line and I have no idea what it does
        # GRM_Nt[(1 + ntotal * (j - 1)):(j * ntotal), (1 + ntotal * (i - 1)):(i * ntotal)] = as.matrix(dcast(melt(tab4, id = c("ind1", "ind2"))[, V := substr(variable, 2, 2)][V == k], ind1 ~ ind2, value.var = "value")[, ind1 := NULL])
        # irng = seq(1 + ntotal * (j - 1), j * ntotal)
        # jrng = seq(1 + ntotal * (i - 1), i * ntotal)

        irng <- (1:ntotal) + (j - 1) * ntotal
        jrng <- (1:ntotal) + (i - 1) * ntotal

        if (DEBUG) {
            # show the loop and what irng and jrng are doing
            cat("loop: (i_", i, ", j_", j, ")\n",
                "  dim GRM_Nt  = ", dim(GRM_Nt)[1], " x ", dim(GRM_Nt)[2], "\n",
                "  irng = [", irng[1], "..", tail(irng, n=1), "]\n",
                "  jrng = [", jrng[1], "..", tail(jrng, n=1), "]\n",
                sep = "")

            # figure out what the heck the following bit was doing
            tmp1 <- melt(tab4, id = c("ind1", "ind2"))
            tmp2 <- tmp1[, V := substr(variable, 2, 2)]
            tmp3 <- tmp2[V == k]
            tmp4 <- dcast(tmp3, ind1 ~ ind2, value.var = "value")
            tmp5 <- tmp4[, ind1 := NULL]
            GRM_Nt[irng, jrng] <- as.matrix(tmp5)
        } else {
            GRM_Nt[irng, jrng] <- as.matrix(
                dcast(melt(tab4, id = c("ind1", "ind2"))[, V := substr(variable, 2, 2)][V == k],
                      ind1 ~ ind2, value.var = "value")[, ind1 := NULL]
            )
        }
    }}
    # GRM_Nt for any number of traits is now done.

    # *** MAKE THIS AN OPTION ***

    summary(GRM_Nt)
    diag(GRM_Nt) <- diag(GRM_Nt) + 0.001

    # Creating breeding values for all traits
    #---------------------------------------------------------------------------

    # Requires GRM_Nt, tmeans, ntotal.
    # GRM_NtR <- mvrnorm(1, rep(tmeans, each = ntotal), GRM_Nt)
    GRM_NtR <- mvrnorm(n = 1,
                       mu = rep(tmeans, each = ntotal),
                       Sigma = GRM_Nt)

    # Dividing up the breeding values into the individual traits
    #---------------------------------------------------------------------------

    # Requires GRM_NtR, ntotal, ntraits.

    # NEW CODE
    # traitvals <- matrix(0, nrow = ntotal, ncol = ntraits)
    # for (i in 1:ntraits) {
    #     traitvals[, i] <- GRM_NtR[seq(ntotal) + (i - 1) * ntotal]
    # }
    traitvals <- matrix(GRM_NtR, nrow = ntotal, ncol = ntraits)

    # Adding environmental (co)variation
    #---------------------------------------------------------------------------
    sigmaE <- matrix(tab1[, traitEcov], nrow = ntraits)
    Etraitvals <- mvrnorm(n = ntotal,
                          mu = rep(0, ntraits),
                          Sigma = sigmaE)


    # Now the overall trait values

    # *** I HAVE TO CALCULATE THIS BEFORE THE TRANSFORMATIONS ***

    # Needs to allow any number of traits.
    # tottraitvals <- matrix(nrow = ntotal, ncol = ntraits)
    # for (i in 1:ntraits) {
    #     tottraitvals[, i] = traitvals[, i] + Etraitvals[, i]
    # }
    tottraitvals1 <- matrix(traitvals + Etraitvals,
                           nrow = ntotal, ncol = ntraits)

    if ("lognormal" %in% trait_transform) {
        cat("- adding lognormal\n")
        tottraitvals_ln <- exp(tottraitvals)
        tottraitvals_ln <- data.table(tottraitvals_ln)
        names(tottraitvals_ln) <- paste(traitnames, "ln", sep = "_")
    }

    if ("invprobit" %in% trait_transform) {
        cat("- adding invprobit\n")
        tottraitvals_ip <- pnorm(tottraitvals)
        tottraitvals_ip <- data.table(tottraitvals_ip)
        names(tottraitvals_ip) <- paste(traitnames, "ip", sep = "_")
    }

    if ("invcloglog" %in% trait_transform) {
        cat("- adding invcloglog\n")
        tottraitvals_icll <- invcloglog(tottraitvals)
        tottraitvals_icll <- data.table(tottraitvals_icll)
        names(tottraitvals_icll) <- paste(traitnames, "icll", sep = "_")
    }

    if ("invlogit" %in% trait_transform) {
        cat("- adding invlogit\n")
        tottraitvals_il <- plogis(tottraitvals)
        tottraitvals_il <- data.table(tottraitvals_il)
        names(tottraitvals_il) <- paste(traitnames, "il", sep = "_")
    }

    #---------------------------------------------------------------------------

    # Join together

    traitvals <- data.table(traitvals)
    names(traitvals) <- paste(traitnames, "BV", sep = "_")
    Etraitvals <- data.table(Etraitvals)
    names(Etraitvals) <- paste(traitnames ,"EV", sep = "_")
    tottraitvals <- data.table(tottraitvals)
    names(tottraitvals) <- traitnames

    traitvals <- cbind(traitvals,Etraitvals,tottraitvals)

    if ("lognormal" %in% trait_transform)
        traitvals <- cbind(traitvals, tottraitvals_ln)

    if ("invprobit" %in% trait_transform)
        traitvals <- cbind(traitvals, tottraitvals_ip)

    if ("invcloglog" %in% trait_transform)
        traitvals <- cbind(traitvals, tottraitvals_icll)

    if ("invlogit" %in% trait_transform)
        traitvals <- cbind(traitvals, tottraitvals_il)

    # id <- rownames(GRM) # seq(ntotal)
    traitvals <- cbind(pedigree[, c("id", "dso")], traitvals)
    #traitvals[, id := as.factor(id)] # FIXME: I don't think this is necessary
    # traitvals[, id := as.numeric(id)]

    # if (!is.null(pedigree))
        # traitvals[, dso := pedigree[, dso]]

    return(traitvals)
}
