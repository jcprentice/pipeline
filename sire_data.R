# datatable
# ------------------------------------------------------------------------------

# this is to create an actual data.table we can work with
dt_vals <- setDT(read.table("data/datatable.tsv", header = FALSE))
names(dt_vals) <- c("id", "group", "It", "Rt", "qg", "qf", "qr")


# this simply creates a string to write to an XML file
dt_chris_xml <- structure(
    # list(paste(readLines("data/datatable.tsv"), collapse = "\n")),
    list(data_table_to_tsv_string(dt_vals)),
    id = 1, group = 2, It = 3, Rt = 4, qg = 5, qf = 6, qr = 7
)



# Ainv_nonzero
# ------------------------------------------------------------------------------

# this is an Nx3 matrix
A1 <- setDT(read.table("data/Ainv_nonzero.tsv", header = FALSE))
A2 <- sparseMatrix(
    i = A1$V1,
    j = A1$V2,
    x = A1$V3,
    index1 = FALSE,
    giveCsparse = FALSE
)


Ainv_vals <- matrix(
    scan("data/Ainv_nonzero.tsv", sep = "\t"),
    ncol = 3, byrow = TRUE
)
 
Ainv_chris_xml <- structure(
    list(paste(readLines("data/Ainv_nonzero.tsv"), collapse = "\n"))
)


