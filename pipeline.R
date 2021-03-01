# Full Pipeline
#
# Send simulated data to SIRE 2.0
# Eventually intend to replace simulated data with Turbot data
#
# R functions derived from Richard Bailey's originals
# SIRE 2.0 modified from Chris Pooley's original

# ------------------------------------------------------------------------------

# load libraries and source files

{
    library(data.table) # used for all data table manipulations
    library(Matrix)     # allows Matrix
    library(xml2)       # used to generate XML file
    # library(corpcor)    # efficient estimation of Cov and (partial) Corr
    library(MASS)       # needed for mvnorm
    # library(MCMCglmm)
    library(codetools)  # check for global variables
    library(here)       # useful according to Hadley Wickham
}

{
    source("get_parameters.R")
    source("make_pedigree.R")
    source("make_grm.R")
    source("make_traits.R")
    source("group_assignment.R")
    # source("episim3.R")
    source("group_assignment.R")
    source("simulate_epidemic.R")
    source("discretise.R")
    source("data_table_to_tsv_string.R")
    source("generate_sire_xml.R")
    # source("sire_data.R")
}

# test functions for global variables
# intersect(findGlobals(model_SIR), ls(envir=.GlobalEnv))


pars      <- get_parameters("SIR")
# pars      <- get_parameters("SEIR")
pedigree  <- make_pedigree(pars)
GRM       <- make_grm(pedigree)
traitdata <- make_traits(pars, pedigree, GRM)

str(pedigree)
str(traitdata)


# Simulate an SIR epidemic
epidemic <- simulate_epidemic(pars, traitdata)
str(epidemic)

pop   <- epidemic$pop
SIRts <- epidemic$SIRts

# discretised version of pop (Tinf, Trec -> "[S,1][I,2][R,3]")
dpop <- discretise(pop, pars)


# Tidy up ready for saving data to XML file
# ------------------------------------------------------------------------------

# setorder(pop, id)
parent_traits <- traitdata[1:pars$nsires, .(id, susceptibility_BV, infectivity_BV, recoverability_BV)]
parent_traits[, c("group", "Tinf", "Trec") := "."]

progeny_traits <- pop[, .(id, group, Tinf, Trec, susceptibility_BV, infectivity_BV, recoverability_BV)]
# find end time for model
tmax <- ceiling(max(progeny_traits[, .(Tinf, Trec)], na.rm = TRUE))

# need to turn NAs in Tinf and Trec to "no"
progeny_traits[, Tinf := as.character(Tinf)]
progeny_traits[, Trec := as.character(Trec)]
progeny_traits[is.na(Tinf), Tinf := "no"]
progeny_traits[is.na(Trec), Trec := "no"]

# these are the columns we want for the data table
dt_generated <- rbind(
    parent_traits[, .(id, group, Tinf, Trec, susceptibility_BV, infectivity_BV, recoverability_BV)],
    progeny_traits[, .(id, group, Tinf, Trec, susceptibility_BV, infectivity_BV, recoverability_BV)],
    fill = TRUE
)

# id should be in the form of [Ind0, ..., Ind99]
dt_generated[, id := paste0("Ind", seq(pars$nsires + pars$nprogeny) - 1)]


# Generate the XML file for SIRE 2.0
# ------------------------------------------------------------------------------

generate_sire_xml(dt_generated, GRM, pars)

# copy the XML file to the SIRE 2.0 folder to run
system("cp foo.xml ../sire2/")

# Run SIRE 2.0 with foo.xml
cat("Running SIRE 2.0 ...\n")
system("../sire2/sire foo.xml 0")


# unmixed
# [1 1 1] [2 2 2] [3 3 3]

# random
# same sires for susceptibility
# different sires for infectivity
# [1 1 2] [2 2 3] [3 1 2]

# striped
# [1 2 3] [1 2 3] [1 2 3]
#
