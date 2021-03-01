library(purrr)

# Functional Programming in R

# Find all the CSV files in the current path
paths <- dir(pattern = "\\.xml$")

# and read them in as data frames
data <- vector("list", length(paths))
for (i in seq_along(paths)) {
    # use [[ whenever you have a single element
    data[[i]] <- read.csv(paths[[i]])
}

# The FP equivalent is much shorter
data <- map(paths, read.csv)
# and has convenient extensions
data <- map_dfr(paths, read.csv, id = "path")

