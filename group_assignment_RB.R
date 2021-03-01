# Function for random assignment of individuals into groups and as index cases
#-------------------------------------------------------------------------------

# This version is just copying Osvaldo's code for randomly assigning groups and
# randomly assigning a single index case to each group#

group.assign = function(pop, ngroups) {
    N <- nrow(pop)
    gr.size <- N / ngroups
    pop[, IndexCase_0 := sample(c(rep(0, ngroups), rep(1, N - ngroups)))]
    # Ones are randomly selected index cases.
    pop[IndexCase_0 == 0, group := 1:ngroups]
    pop[IndexCase_0 == 1, group := sample(rep(1:ngroups, gr.size - 1))]
    pop[IndexCase_0 == 0, Tinf := 0]
    
    return(pop)
}

