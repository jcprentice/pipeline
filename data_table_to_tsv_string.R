# Format a data table as a TSV file for including in the XML file

data_table_to_tsv_string <- function(dt) {
    paste0(paste0(
        capture.output(
            write.table(
                dt, "", sep = "\t", eol = "\n",
                quote = FALSE, row.names = FALSE, col.names = FALSE
            )
        ),
        collapse = "\n"
    ),
    "\n")
}

# test
#
# dt <- data.table(id = sample(1:100, size = 5),
#                  x = rbinom(5, 3, 1/6),
#                  y = runif(5))
# dt_tsv <- data_table_to_tsv_string(dt)
# dt
# dt_tsv
