make_grm <- function(pedigree, pars) {
    cat("Making the GRM ...\n")

    tmp <- as.data.table(pedigree[, expand.grid(id, id)])
    setnames(tmp, c("Var1", "Var2"), c("id1", "id2"))

    tmp[, c("sire1", "dam1") := pedigree[tmp$id1, c("sire", "dam")]]
    tmp[, c("sire2", "dam2") := pedigree[tmp$id2, c("sire", "dam")]]

    tmp[, related := 0]
    # full sibs
    tmp[sire1 == sire2 & dam1 == dam2, related := 1/2]
    # half sibs
    tmp[sire1 == sire2 & dam1 != dam2 & !is.na(dam1) & !is.na(dam2), related := 1/4]
    tmp[sire1 != sire2 & !is.na(sire1) & !is.na(sire2) & dam1 == dam2, related := 1/4]
    # parent/child
    tmp[id1 == sire2 | id1 == dam2 | id2 == sire1 | id2 == dam1, related := 1/2]
    # same individual
    tmp[id1 == id2, related := 1]

    relatemat <- dcast(tmp, id1 ~ id2, value.var = "related")
    relatemat[, id1 := NULL]

    GRM <- as.matrix(relatedmat)
    colnames(GRM) <- pedigree[, id]
    rownames(GRM) <- pedigree[, id]

    return(GRM)
}


