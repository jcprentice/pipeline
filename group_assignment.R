# Randomly assign individuals into groups and as index cases

group_assign <- function(pop, ngroups) {
    N <- nrow(pop)

    dt <- data.table(
        # shuffle id to randomise
        id = sample(pop$id),

        # ngroups doesn't have to divide N evenly
        group = factor(rep(1:ngroups, times = ngroups, length.out = N)),

        # 1 index case per group, denoted "1", the rest are "0"
        index_case = c(rep(1, ngroups), rep(0, N - ngroups))
    )
    setorder(dt, id)

    pop <- merge(pop, dt, by = "id")
    return(pop)
}
