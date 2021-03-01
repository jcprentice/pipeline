# Full Pipeline
#
# Send simulated data to SIRE 2.0
# Eventually intend to replace simulated data with Turbot data
#
# R functions based on Richard Bailey's originals
# SIRE 2.0 modified from Chris Pooley's original

# ------------------------------------------------------------------------------

# load libraries and source files
using Parameters

library(data.table) # used for all data table manipulations
library(xml2)       # used to generate XML file
library(corpcor)    # efficient estimation of Cov and (partial) Corr
library(MASS)
library(MCMCglmm)

setwd("~/Dropbox/Roslin/pipeline")

source("parameters.jl")
source("pedigree_grm-old.R")
source("pedigree_grm.R")
source("make_traits.R")
source("episim3.R")
source("episim3-new.R")
source("group_assignment.R")
source("discrete_time.R")
source("data_table_to_tsv_string.R")

sp = SirePars()

# Quickly print out what we've got
summarise_pars(sp)

# ------------------------------------------------------------------------------

# Create pedigree and GRM

cat("Running ped.to.grm() ...\n")
ped, ped2, sire = pedigree_sire(sp)

pedtest.fullhalfsib.old = ped.to.grm.old(
    nsire = nsire,
    dpsire = dpsire,
    oppair = oppair,
    # remove individuals with no known parentage
    trim_parents = FALSE,
    # assume those with no identified parents are unrelated.
    # sets relatedness to 0 if TRUE, otherwise NA.
    noparents.unrelated = TRUE
)

pedtest.fullhalfsib = ped.to.grm(
    nsire = nsire,
    dpsire = dpsire,
    oppair = oppair,
    # remove individuals with no known parentage
    trim_parents = FALSE,
    # assume those with no identified parents are unrelated.
    # sets relatedness to 0 if TRUE, otherwise NA.
    noparents.unrelated = TRUE
)

str(pedtest.fullhalfsib)

pedigree_old = pedtest.fullhalfsib.old$pedigree
pedigree_old
GRM1 = pedtest.fullhalfsib.old$GRM
Matrix(GRM1)
pedigree = pedtest.fullhalfsib$pedigree
pedigree
GRM = pedtest.fullhalfsib$GRM
Matrix(GRM)

# ------------------------------------------------------------------------------

# Generate trait values
# Note: can be skipped

cat("Running make.tvals.int() ...\n")
traitdata = make.tvals.int(
    nind = nind,
    # The genomic relationship matrix excluding parents, created above
    GRM = GRM,
    ntraits = ntraits,
    tmeans = tmeans,
    tVA = tVA,
    trhoG = trhoG,
    tVE = tVE,
    trhoE = trhoE,
    traitnames = traitnames,
    # What trait transformations to carry out. Options:
    # "lognormal","invprobit","invcloglog","invlogit"
    trait.transform = "lognormal",
    # This determines how the function ensures the matrix is positive
    # definite. A GRM can also be made from genotype data, which takes a
    # different form (it always has one zero eigenvalue)
    GRM.from.pedigree = TRUE,
    DEBUG = FALSE
)

# IMPLEMENT A TEST HERE TO CHECK THIS ARE AS WE EXPECT!!
# correlations vs GRM
# with a GLMM

# Look at the output (run the line twice if it doesn't show up first time):
traitdata


# ------------------------------------------------------------------------------

# Simulate an SIR epidemic

# If we've included parents in the GRM, we need to lose them for the moment in
# orger to run EpiSim3(). We'll recombine it with the sires afterwards.
if (nrow(GRM) == nparent + nprogeny) {
    subGRM = GRM[(nparent + 1):nrow(GRM), (nparent + 1):nrow(GRM)]
} else {
    subGRM = GRM
}

# this really should have been done sooner
traitdata[1:nsire, dso := "sire"]
traitdata[(nsire+1):(nsire+ndam), dso := "dam"]
traitdata[(nsire+ndam+1):(nsire+ndam+nprogeny), dso := "offspring"]

cat("Running EpiSim3() ...\n")
epitest_sir = EpiSim3(
    nind = nind,
    group.and.index = NULL,
    ngroups = ngroups,
    PairwiseBeta = PairwiseBeta,
    RRgamma = RRgamma,
    pop = traitdata,
    offspring.only = TRUE,
    GRM = subGRM,
    GRM.from.pedigree = TRUE,
    ntraits = ntraits,
    tmeans = tmeans,
    tVA = tVA,
    trhoG = trhoG,
    tVE = tVE,
    trhoE = trhoE,
    traitnames = traitnames,
    trait.transform = c("lognormal"),
    frailties = "lognormal",
    Epi.type = "SIR",
    DEBUG = FALSE
)

# epitest_sir = EpiSim3_new(
#     nind = nind,
#     group.and.index = NULL,
#     ngroups = ngroups,
#     PairwiseBeta = PairwiseBeta, RRgamma = RRgamma,
#     pop = traitdata,
#     offspring.only = FALSE,
#     GRM = subGRM,
#     GRM.from.pedigree = TRUE,
#     ntraits = ntraits,
#     tmeans = tmeans,
#     tVA = tVA,
#     trhoG = trhoG,
#     tVE = tVE,
#     trhoE = trhoE,
#     traitnames = traitnames,
#     trait.transform = "lognormal",
#     frailties = "lognormal",
#     Epi.type = "SIR",
#     DEBUG = FALSE
# )

pop   = epitest_sir$pop
SIRts = epitest_sir$SIRts

dpop = discretize.time(
    pop,
    nind,
    ngroups,
    interval.duration = "event",
    truncate = TRUE, # need this to be F for Chris's data
    Epi.type = "SIR" # should be SIR
)

# setorder(pop, id)
parent_traits = traitdata[1:nsire, .(id, susceptibility_BV, infectivity_BV, recoverability_BV)]
parent_traits[, group := "."]
parent_traits[, Tinf := "."]
parent_traits[, Trec := "."]

progeny_traits = pop[, .(id, group, Tinf, Trec, susceptibility_BV, infectivity_BV, recoverability_BV)]
# find end time for model
tmax = ceiling(max(progeny_traits[, .(Tinf, Trec)], na.rm = TRUE))

# need to turn NAs in Tinf and Trec to "no"
progeny_traits[, Tinf := as.character(Tinf)]
progeny_traits[, Trec := as.character(Trec)]
progeny_traits[is.na(Tinf), Tinf := "no"]
progeny_traits[is.na(Trec), Trec := "no"]

# these are the columns we want for the data table
dt_generated = rbind(
    parent_traits[, .(id, group, Tinf, Trec, susceptibility_BV, infectivity_BV, recoverability_BV)],
    progeny_traits[, .(id, group, Tinf, Trec, susceptibility_BV, infectivity_BV, recoverability_BV)],
    fill = TRUE
)

# id should be in the form of [Ind0, ..., Ind99]
dt_generated[, id := paste0("Ind", seq(0, nsire + nprogeny - 1))]


# Generate the XML file for SIRE 2.0
# ------------------------------------------------------------------------------

source("generate_sire_xml.R")

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
