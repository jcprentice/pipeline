# General function for converting a pedigree to a GRM
#-------------------------------------------------------------------------------

# Returns both a pedigree and a GRM as outputs

# The analysis package MCMCglmm uses a pedigree in standard format, but can also
# take a GRM from either a pedigree or genomic data

ped.to.grm <- function(
    # An existing pedigree
    pedigree = NULL,
    # No. of sires if no existing pedigree
    nsires = NULL,
    # Dams per sire if no existing pedigree
    dpsire = NULL,
    # Offspring per pair if no existing pedigree
    oppair = NULL,
    # Remove all individuals with no known parentage
    trim_parents = TRUE,
    # Whether to assume all those with no identified parents are unrelated. Sets
    # relatedness to 0 if TRUE, otherwise NA.
    noparents.unrelated = FALSE
) {
    # handy values
    ndams     <- nsires * dpsire
    nparents  <- nsires + ndams
    nprogeny <- nsires * dpsire * oppair
    opsire   <- dpsire * oppair
    
    if (is.null(pedigree)) {
        pedi_parents <- data.table(id = 1:nparents)
        pedi_offspring <- data.table(
            id = (1 + nparents):(nprogeny + nparents),
            dam = rep((1 + nsires):(ndams + nsires), each = oppair),
            sire = rep(1:nsires, each = opsire)
        )
        pedi <- rbind(pedi_parents, pedi_offspring, fill = T)
    } else {
        pedi <- as.data.table(pedigree)
        pedi[, id := as.integer(as.character(id))]
        pedi[, sire := as.integer(as.character(sire))]
        pedi[, dam := as.integer(as.character(dam))]
    }
    
    if (trim_parents) {
        pedi <- pedi[is.na(dam) == FALSE]
    } else if(noparents.unrelated) {
        pedi[is.na(dam), dam := 0]
        pedi[is.na(sire), sire := 0]
    }
    
    pedi3 <- as.data.table(pedi[, expand.grid(id, id)])
    setnames(pedi3, c("Var1", "Var2"), c("iding1", "iding2"))
    
    setkey(pedi, id)
    setkey(pedi3, iding1)
    pedi4 <- pedi[pedi3]
    setnames(pedi4, c("sire", "dam", "id"), c("sire1_1", "dam1_2", "iding1"))
    
    setkey(pedi, id)
    setkey(pedi4, iding2)
    pedi5 <- pedi[pedi4]
    setnames(pedi5, c("sire", "dam", "id"), c("sire2_1", "dam2_2", "iding2"))
    
    # Same individual (matrix diagonal)
    pedi5[sire1_1 == sire2_1 & dam1_2 == dam2_2 & iding1 == iding2,
          related := 1]
    # First parent-parent
    # Unrelated parents
    pedi5[dam2_2 == 0 & sire2_1 == 0 & dam1_2 == 0 & sire1_1 == 0 & is.na(related),
          related := 0]
    # Now add parent-offspring relationships
    pedi5[dam2_2  == iding1, related := 0.5]
    pedi5[sire2_1 == iding1, related := 0.5]
    pedi5[dam1_2  == iding2, related := 0.5]
    pedi5[sire1_1 == iding2, related := 0.5]
    
    # Full sibs
    pedi5[sire1_1 == sire2_1 & dam1_2 == dam2_2 & iding1 != iding2 & sire1_1 != 0 & sire2_1 != 0 & dam1_2 != 0 & dam2_2 != 0,
          related := 0.5]
    # Half sibs
    pedi5[sire1_1 == sire2_1 & dam1_2 != dam2_2 & sire1_1 != 0 & sire2_1 != 0 & dam1_2 != 0 & dam2_2 != 0,
          related := 0.25]
    # Half sibs
    pedi5[sire1_1 != sire2_1 & dam1_2 == dam2_2 & sire1_1 != 0 & sire2_1 != 0 & dam1_2 != 0 & dam2_2 != 0,
          related := 0.25]
    # Unrelated offspring
    pedi5[sire1_1 != sire2_1 & dam1_2 != dam2_2 & sire1_1 != 0 & sire2_1 != 0 & dam1_2 != 0 & dam2_2 != 0,
          related := 0]
    #  The remainder should be unrelated parents and offspring
    pedi5[is.na(related),
          related := 0]
    
    relatemat <- dcast(pedi5, iding1 ~ iding2, value.var = "related")
    relatemat[, iding1 := NULL]
    
    relatemat2 <- as.matrix(relatemat)
    colnames(relatemat2) = pedi[, id]
    rownames(relatemat2) = pedi[, id]
    
    pedi[id %in% pedigree[, unique(na.omit(sire))],
              dso := "sire"]
    pedi[id %in% pedigree[, unique(na.omit(dam))],
              dso := "dam"]
    pedi[id %in% pedigree[, unique(na.omit(sire))] == F &
                  id %in% pedigree[, unique(na.omit(dam))] == F,
              dso := "offspring"]
    
    
    # copy of pedi but with 0s as NAs.
    pedi2 <- copy(pedi)
    pedi2[dam == "0", dam := NA]
    pedi2[sire == "0", sire := NA]

    # return
    list(pedigree = pedi,
         pedigree2 = pedi2,
         GRM = relatemat2)
}


