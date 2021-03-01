make_pedigree <- function(pars) {
    cat("Making pedigree ...\n")

    with (pars, {
        # create an empty data table
        ped <- data.table(id = 1:ntotal,
                          sire = as.numeric(NA), dam = as.numeric(NA),
                          dso = factor("", levels = c("sire", "dam", "offspring")))

        sire_ids   <- seq(nsires)
        dam_ids    <- seq(ndams) + nsires
        parent_ids <- seq(nparents)
        prog_ids   <- seq(nprogeny) + nparents

        # set type
        ped[sire_ids, dso := "sire"]
        ped[dam_ids,  dso := "dam"]
        ped[prog_ids, dso := "offspring"]
        ped[prog_ids, dso := "offspring"]

        # set sires and dams
        ped[prog_ids, sire := rep(sire_ids, each = dpsire * oppair)]
        ped[prog_ids, dam := rep(dam_ids, each = oppair)]

        return(ped)

        # tidy up
        # ped[, id := as.factor(id)]
        # ped[, sire := as.factor(sire)]
        # ped[, dam := as.factor(dam)]
        # ped[, dso := as.factor(dso)]

        # duplicate this but set parents' values to NA
        # ped2 <- copy(ped)
        # ped2[parent_ids, c("sire", "dam") := NA]

        # return(list(pedigree = ped, pedigree2 = ped2))
    })
}
