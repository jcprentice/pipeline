# Bailey breeding value simulation function
# 21 April 2020
#-------------------------------------------------------------------------------

# Can be called from within the EpiSim3 function or run separately.

make.tvals.int <- function(
    # number of individuals
    nind,
    # genomic relationship matrix
    GRM,
    # can be useful for adding the "sire/dam/offspring" column
    pedigree = NULL,
    # no. of traits (3 for SIR)
    ntraits,
    # trait means
    tmeans = rep(0, ntraits),
    # trait variances and correlations
    tVA, trhoG, tVE, trhoE,
    traitnames = c("susceptible, infected, recovered"),
    trait.transform = c("lognormal", "invprobit", "invcloglog", "invlogit"),
    # This determines how the function ensures the matrix is positive definite.
    # A GRM can also be made from genotype data, which takes a different form
    # (it always has one zero eigenvalue)
    GRM.from.pedigree = FALSE,
    # turn on debug information
    DEBUG = FALSE
) {
    
    require(data.table)
    require(MASS)
    require(corpcor)
    
    # overwrite things here
    GRM[is.na(GRM)] <- 0
    nind <- nrow(GRM)
    
    # I'm creating a melted table of pairwise trait variances and covariances
    # Now I have all the trait variances and covariances. From this I should be
    # able to create GRM_crosst and GRM*tVA
    tab1 <- data.table(expand.grid(trait1 = 1:ntraits, trait2 = 1:ntraits),
                       corrval = as.numeric(combn(trhoG, 1)))
    tab1[, tvar1 := tVA[trait1]]
    tab1[, tvar2 := tVA[trait2]]
    tab1[, traitcov := corrval * sqrt(tvar1) * sqrt(tvar2)]
    
    # Environmental covariance
    tab1[, Ecorrval := as.numeric(combn(trhoE, 1))]
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
                     nrow = (ntraits * nind),
                     ncol = (ntraits * nind))
    k <- 0
    for (i in 1:ntraits) {
        for (j in 1:ntraits) {
            # i = columns, j = Rows
            k = k + 1
            # this is a horrible line and I have no idea what it does
                same_sire <- FALSE
                same_dam <- FALSE
            # TODO: rewrite it to make it less evil
            
            # GRM_Nt[(1 + nind * (j - 1)):(j * nind), (1 + nind * (i - 1)):(i * nind)] = as.matrix(dcast(melt(tab4, id = c("ind1", "ind2"))[, V := substr(variable, 2, 2)][V == k], ind1 ~ ind2, value.var = "value")[, ind1 := NULL])
            irng = (1 + nind * (j - 1)):(j * nind)
            jrng = (1 + nind * (i - 1)):(i * nind)
            
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
                GRM_Nt[irng, jrng] = as.matrix(tmp5)
            } else {
                GRM_Nt[irng, jrng] = as.matrix(
                    dcast(melt(tab4, id = c("ind1", "ind2"))[, V := substr(variable, 2, 2)][V == k],
                          ind1 ~ ind2, value.var = "value")[, ind1 := NULL]
                )
            }
        }
    }
    # GRM_Nt for any number of traits is now done.
    
    # *** MAKE THIS AN OPTION ***
    
    summary(GRM_Nt)
    if (GRM.from.pedigree) {
        # GRM_Nt = make.positive.definite(GRM_Nt)
    } else {
        diag(GRM_Nt) = diag(GRM_Nt) + 0.001
    }
    
    # Creating breeding values for all traits
    #---------------------------------------------------------------------------
    
    # Requires GRM_Nt, tmeans, nind.
    GRM_NtR = mvrnorm(1,
                      mu = rep(tmeans, each = nind),
                      Sigma = GRM_Nt)
    
    # Dividing up the breeding values into the individual traits
    #---------------------------------------------------------------------------
    
    # Requires GRM_NtR, nind, ntraits.
    
    # NEW CODE
    traitvals = matrix(nrow = nind, ncol = ntraits)
    for (i in 1:ntraits) {
        traitvals[, i] = GRM_NtR[(1 + nind * (i - 1)):(i * nind)]
    }
    
    # Adding environmental (co)variation
    #---------------------------------------------------------------------------
    sigmaE = matrix(tab1[, traitEcov], nrow = ntraits)
    Etraitvals = mvrnorm(nind,
                         mu = rep(0, ntraits),
                         Sigma = sigmaE)
    
    # Now the overall trait values
    
    # *** I HAVE TO CALCULATE THIS BEFORE THE TRANSFORMATIONS ***
    
    # Needs to allow any number of traits.
    tottraitvals = matrix(nrow = nind, ncol = ntraits)
    for (i in 1:ntraits) {
        tottraitvals[, i] = traitvals[, i] + Etraitvals[, i]
    }
    
    if ("lognormal" %in% trait.transform) {
        tottraitvals_ln = exp(tottraitvals)
        tottraitvals_ln = data.table(tottraitvals_ln)
        names(tottraitvals_ln) = paste(traitnames, "ln", sep = "_")
    }
    
    if ("invprobit" %in% trait.transform) {
        tottraitvals_ip = pnorm(tottraitvals)
        tottraitvals_ip = data.table(tottraitvals_ip)
        names(tottraitvals_ip) = paste(traitnames, "ip", sep = "_")
    }
    
    if ("invcloglog" %in% trait.transform) {
        tottraitvals_icll = invcloglog(tottraitvals)
        tottraitvals_icll = data.table(tottraitvals_icll)
        names(tottraitvals_icll) = paste(traitnames, "icll", sep = "_")
    }
    
    if ("invlogit" %in% trait.transform) {
        tottraitvals_il = plogis(tottraitvals)
        tottraitvals_il = data.table(tottraitvals_il)
        names(tottraitvals_il) = paste(traitnames, "il", sep = "_")
    }
    
    #---------------------------------------------------------------------------
    
    # Join together
    
    traitvals = data.table(traitvals)
    names(traitvals) = paste(traitnames, "BV", sep="_")
    Etraitvals = data.table(Etraitvals)
    names(Etraitvals) = paste(traitnames ,"EV", sep="_")
    tottraitvals = data.table(tottraitvals)
    names(tottraitvals) = traitnames
    
    traitvals = cbind(traitvals,Etraitvals,tottraitvals)
    
    if ("lognormal" %in% trait.transform) {
        traitvals = cbind(traitvals, tottraitvals_ln)
    }
    
    if ("invprobit" %in% trait.transform) {
        traitvals = cbind(traitvals, tottraitvals_ip)
    }
    
    if ("invcloglog" %in% trait.transform) {
        traitvals = cbind(traitvals, tottraitvals_icll)
    }
    
    if ("invlogit" %in% trait.transform) {
        traitvals = cbind(traitvals, tottraitvals_il)
    }
    
    # id = rownames(GRM) # seq(1, nind) *********************
    id = 1:nrow(GRM)
    traitvals = cbind(id, traitvals)
    #traitvals[, id := as.factor(id)] # FIXME: I don't think this is necessary
    traitvals[, id := as.numeric(id)]
    
    if (is.null(pedigree) == FALSE) {
        print("pedigree is not NULL")
        traitvals[id %in% pedigree[, unique(na.omit(sire))],
                  dso := "sire"]
        traitvals[id %in% pedigree[, unique(na.omit(dam))],
                  dso := "dam"]
        traitvals[id %in% pedigree[, unique(na.omit(sire))] == F &
                      id %in% pedigree[, unique(na.omit(dam))] == F,
                  dso := "offspring"]
    } else {
        print("pedigree is NULL")
    }
    
    
    return(traitvals)
}
