# Richard Ian Bailey
# 12 August 2020
#
# Example code for simulating quantitative genetic epidemiological trait values
# and epidemics, and GLMM analysis
#-------------------------------------------------------------------------------

# Install these packages first before loading if they're not already installed

library(data.table) # I use this for all data table manipulations
library(corpcor)
library(MASS)
library(MCMCglmm)

# Set working directory (network Wilson group folder):
# Jamie set this to wherever you put the function code files.

# setwd("Z:/wilson_group/Fishboost/Richard_Simulations_2020/Report_10Aug2020")
setwd("~/Dropbox/Roslin/pipeline")

# source all relevant files

source("pedigree_grm.R")
source("episim3.R")
source("internal_bv.R")
source("group_assignment.R")
source("discrete_time.R")


# First look at a couple of existing data files produced by the function for
# simulating epidemics:
#-------------------------------------------------------------------------------

# Jamie you can remove the file path if these files are in your working directory

# pop = fread("Z:/wilson_group/Fishboost/Richard_Simulations_2020/Report_10Aug2020/Generic_simulated_data_files/fs_B01_h1_V01_gr10_1_pop.csv")
# SIRts= fread("Z:/wilson_group/Fishboost/Richard_Simulations_2020/Report_10Aug2020/Generic_simulated_data_files/fs_B01_h1_V01_gr10_1_SIRts.csv")

# Or just:

pop = fread("data/fs_B01_h1_V01_gr10_1_pop.csv")
SIRts = fread("data/fs_B01_h1_V01_gr10_1_SIRts.csv")

# This contains the trait values and time to infection for each individual.
pop

# susceptibility_ln & infectivity_ln are the lognormally distributed trait
# values. susceptibility_BV & infectivity_BV are the heritable component of the
# normally distributed trait values, i.e. the 'breeding values'.

# Contains the epidemic data, including the numbers of susceptible and infected
# in each group and overall at each infection event.
SIRts




# Open, select and run all the functions to load them into the workspace
#-------------------------------------------------------------------------------

# Creates a genomic relationship matrix from a pedigree, and also creates both a
# GRM and a pedigree from specified relationships.
ped.to.grm

# Simulates trait values for multiple traits from a multivariate normal
# distribution with the specified GRM and trait variances and correlations.
# Called internally by the EpiSim3 function for simulating epidemics, but can
# also be run separately if only the trait values are needed.
make.tvals.int

# Assigns individuals in a dataset randomly to the specified number of groups,
# and randomly assigns one index case per group. Called internally by the
# EpiSim3 function if group and index case status are not specified explicitly.
group.assign

# Assigns individuals to groups and as index cases, simulates quantitative
# genetic trait values and runs SI or SIR epidemics using the Gillespie
# algorithm.
EpiSim3

# Splits time units into discrete intervals in preparation for GLMM analysis.
discretize.time



# Create a pedigree and GRM with each family being a mixture of full- and
# half-sibs:
#-------------------------------------------------------------------------------

# Jamie, I think the plan is to use half-sib families only, in which case
# oppair=1.


# First exclude the parents from the pedigree. This is for the simulation of
# trait values:
cat("Running ped.to.grm() ...\n")
pedtest.fullhalfsib = ped.to.grm(
    nsires = 100, # No. of sires
    dpsire = 5, # Dams per sire
    oppair = 2, # Offspring per pair
    trim_parents = T # Whether to remove individuals with no known parentage
)

# This creates both a pedigree and GRM (genomic relationship matrix, with
# pairwise relatdness among individuals):
str(pedtest.fullhalfsib)

# The pedigree:
# NOTE: If nothing appears first time, run this line again. It's just a quirk of
# the data.table package.
pedtest.fullhalfsib$pedigree

# Make the same pedigree and GRM but including the parents this time. This is to
# make a standard format pedigree for GLMM analysis:
pedtest2.fullhalfsib = ped.to.grm(
    nsires = 100,
    dpsire = 20,
    oppair = 1,
    trim_parents = F # fix this s.t. diagonal is correct (1s on diag, 0s o/w)
)

pedtest2.fullhalfsib$pedigree


# Create trait values but don't run any epidemics
#------------------------------------------------

# The make.tvals.int function can be run internally within Episim3 below when
# running epidemics, or you can just create a table of trait values.

# For the additive genetic component ("BV" in the output) trait values are
# sampled from a multivariate normal distribution based on a covariance matrix
# with N individuals * N traits rows and columns, so this can get pretty big!
# 1000 individuals and 3 traits works on my 8Gb RAM laptop, I haven't tried much
# more.

cat("Running make.tvals.int() ...\n")
traitdata = make.tvals.int(
    # Number of individuals
    nind = 1000,
    
    # The genomic relationship matrix excluding parents, created above
    GRM = pedtest.fullhalfsib$GRM,
    
    # Number of traits to simulate. In this case I'm simulating
    # susceptibility, infectivity and recoverability
    ntraits = 3,
    
    # Trait means
    tmeans = rep(0, 3),
    
    # Trait additive genetic variances
    tVA = c(0.5, 0.1, 0.1),
    
    # Trait additive genetic correlation matrix. In this case correlations
    # are zero.
    trhoG = diag(3),
    
    # Trait non-additive (environmental) variance, with no association to
    # the pedigree
    tVE = c(0.3, 0.3, 0.3),
    
    # Trait non-additive correlation matrix
    trhoE = diag(3),
    
    # Names to be assigned to the traits
    traitnames = c("susceptibility", "infectivity", "recoverability"),
    
    # What trait transformations to carry out. Options:
    # "lognormal","invprobit","invcloglog","invlogit"
    trait.transform = "lognormal",
    
    # This determines how the function ensures the matrix is positive
    # definite. A GRM can also be made from genotype data, which takes a
    # different form (it always has one zero eigenvalue)
    GRM.from.pedigree = TRUE
)

# Look at the output (run the line twice if it doesn't show up first time):

traitdata


# The above can be used to make a table of individual trait values e.g. to be
# plugged into your software.


# Simulate data and epidemics
#----------------------------

# Here I'm making the trait values again as part of the epidemic simulation
# setup, independently of the above trait values (I only need susceptibility and
# infectivity variation at the moment).

cat("Running EpiSim3() with SI ...\n")
epitest = EpiSim3(
    # No. of individuals
    nind = 1000,
    
    # For explicitly specifying group and index case status of each
    # individual. Requires a data.frame with 'nind' rows and columns for
    # group (character), donor/recipient ("D" or "R") and individual id
    group.and.index = NULL,
    
    # Use ngroups if you want to randomly assign individuals to groups, with
    # a single randomly assigned index case per group. This will be ignored
    # if 'group.and.index' is set
    ngroups = 5,
    
    # Here I use a single scalar value, but this can also be a vector of
    # 'fixed effect' Beta values the same length as nind and assuming
    # individuals are in the same order as in the pedigree. Could be useful
    # if some variable e.g. body weight affects epidemiological trait values
    PairwiseBeta = 0.1,
    
    # The genomic relationship matrix excluding parents, created above
    GRM = pedtest.fullhalfsib$GRM,
    
    # This determines how the function ensures the matrix is positive
    # definite
    GRM.from.pedigree = TRUE,
    
    # No. of traits to simulate
    ntraits = 2,
    
    # Trait means
    tmeans = rep(0, 2),
    
    # Trait additive genetic variances (here I'm assuming no variation in
    # infectivity)
    tVA = c(0.5, 0),
    
    # Trait additive genetic correlation matrix. In this case correlations
    # are zero.
    trhoG = diag(2),
    
    # Trait non-additive variance
    tVE = c(0.5, 0),
    
    # Trait non-additive correlation matrix
    trhoE = diag(2),
    
    # Names to be assigned to the traits
    traitnames = c("susceptibility", "infectivity"),
    
    # What trait transformations to carry out. Options:
    # "lognormal","invprobit","invcloglog","invlogit"
    trait.transform = c("lognormal"),
    
    # Currently individual frailties are assumed lognormal, but this could be
    # made more flexible, e.g. allowing gamma-distributed epidemiological trait
    # values
    frailties = "lognormal",
    
    # Whether this is an SI or SIR epidemic. More options could be added
    Epi.type = "SI"
)

cat("Running EpiSim3() with SIR ...\n")

epitest_sir = EpiSim3(
    nind = 1000,
    group.and.index = NULL,
    ngroups = 5,
    PairwiseBeta = 0.1,
    RRgamma = 1,
    GRM = pedtest.fullhalfsib$GRM,
    GRM.from.pedigree = TRUE,
    ntraits = 3,
    tmeans = rep(0, 3),
    tVA = c(0.5, 0, 0),
    trhoG = diag(3),
    tVE = c(0.5, 0, 0),
    trhoE = diag(3),
    traitnames = c("susceptibility", "infectivity", "recoverability"),
    trait.transform = c("lognormal"),
    frailties = "lognormal",
    Epi.type = "SIR"
)


# Look at the output
#-------------------

epitest

# The object contains the same two tables loaded at the start,'pop' and 'SIRts',
# but with new data

# We now have the pedigree, trait values and 1 epidemic simulation.




# Split the epidemic into discrete time intervals ready for GLMM analysis
#------------------------------------------------------------------------

cat("Running discretise.time() ...\n")
epitest_disc = discretize.time(
    # The pop file produced above
    pop = epitest$pop,
    
    # No. of individuals
    nind = 1000,
    
    # No. of groups
    ngroups = 5,
    
    # The duration of each discrete time interval. If not specified, a new
    # interval starts at each transition event
    interval.duration = 0.05,
    
    # If TRUE, removes all time points after infection for each individual
    # (none removed for individuals that didn't become infected). Needed for
    # GLMM analysis
    truncate = TRUE,
    
    # Epidemic type, "SI" or "SIR"
    Epi.type = "SI"
)

# This 2nd version is for working with SIRE 2.0
cat("Running discretise.time() for SIRE2.0 ...\n")
epitest_disc2 = discretize.time(
    pop = epitest_sir$pop,
    nind = 1000,
    ngroups = 5,
    interval.duration = 0.05,
    truncate = FALSE, # FALSE to work with SIRE 2.0
    Epi.type = "SI"
)



# Look at the contents
#---------------------

epitest_disc

epitest_disc2


# Includes the same $pop table as previously with a few extra columns, and a new
# $discpop table with one row per individual per time interval. 'pit' is the
# response variable for GLMM analysis, with a 0 for each time interval before
# infection for each individual, and a 1 for the interval during which infection
# occurred. The mean of these for each individual is its poisson mean number of
# infection events per time unit.


# We are now ready to run a GLMM analysis.
#-----------------------------------------

# Use the R package MCMCglmm
#---------------------------


# First specify the fixed, random and residual effects priors
#------------------------------------------------------------

# The fixed effects include the intercept and the force of infection offset,
# which is column 'logoffset_deltat' in the $discpop table.

prior_offset = list(
    # Fixed effects with high variances. The offset is modified to very low variance, below
    B = list(mu = c(0, 1), V = diag(1e6, 2)),
    
    # Residual priors
    R = list(V = 1, nu = 1),
    
    # Priors for additive, nonadditive and group random effects
    G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
             G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
             G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

# Modifying the variance for the offset prior. MCMCglmm doesn't have an 'offset'
# option, so I just minimize the prior variance
diag(prior_offset$B$V)[2] = 1e-8



# Run the GLMM
#-------------

cat("Running MCMCglmm() ...\n")
model_epitest = MCMCglmm(
    # Response and fixed effects
    pit ~ logoffset_deltat,
    
    # Random effects
    random = ~ animal + id + group,
    
    # poisson model. Default log link
    family = "poisson",
    
    # The prior specified above
    prior = prior_offset,
    
    # The pedigree, including sires and dams. By default MCMCglmm recognises
    # the random effect 'animal' as associated with a pedigree
    pedigree = pedtest2.fullhalfsib$pedigree,
    
    # Dataset created above
    data = epitest_disc$discpop,
    
    # I think this is only for binomial models. It improves the MCMC
    # algorithm
    slice = TRUE,
    
    # This prevents the MCMC from searching extreme unrealistic values
    trunc = TRUE,
    
    # Total iterations including burnin
    nitt = 20000,
    
    # Length of burnin
    burnin = 5000,
    
    # Thinning, determining the number of posterior values stored
    thin = 8
)


# Look at the results
#--------------------

summary(model_epitest)

# The intercept is the estimate of log(Beta).

# 'animal' is the estimated additive genetic (i.e. heritable) variance in
# susceptibility.
