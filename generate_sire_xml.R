
# if we want the same data as example.xml
# handy for exploring the structure
#
# sire_ex <- read_xml("example.xml")
# sire_ex_list <- as_list(sire_ex)

generate_sire_xml <- function(dt_generated, GRM, pars) {

    # turn the data table into an XML node
    # ------------------------------------------------------------------------------

    dt_generated_xml <- structure(
        list(data_table_to_tsv_string(dt_generated)),
        id = 1, group = 2, It = 3, Rt = 4, qg = 5, qf = 6, qr = 7
    )


    # GRM
    # ------------------------------------------------------------------------------

    # clip out the dams
    nondams <- c(seq(pars$nsires), seq(pars$nprogeny) + pars$nparents)
    GRM2 <- GRM[nondams, nondams]

    # turn the GRM into a sparse matrix
    sA <- summary(Matrix(GRM2, sparse = TRUE))

    # then turn that into a matrix that we can write
    mA <- matrix(c(sA$i - 1, sA$j - 1, sA$x), ncol = 3)

    # turn this into an XML node
    A_nonzero_xml <- structure(list(data_table_to_tsv_string(mA)))



    # Prediction Accuracies
    # ------------------------------------------------------------------------------

    # build up ind = "Ind0,Ind1,...,Ind99"
    sire_vals <- seq(pars$nsire) - 1
    sire_inds <- paste0(sapply(sire_vals, function(x) paste0("Ind", x)), collapse = ",")
    pa_sire_xml <- structure(list(), name = "sire", ind = sire_inds)

    # build up ind = "Ind100,Ind101,...,Ind2099"
    prog_vals <- seq(pars$nprogeny) + pars$nsires - 1
    prog_inds = paste0(sapply(prog_vals, function(x) paste0("Ind", x)), collapse = ",")
    pa_progeny_xml <- structure(list(), name = "progeny", ind = prog_inds)


    x <- list(
        SIRE = structure(
            version = "1.0",
            list(
                mcmc = structure(list(), nsample = pars$nsample, burnin = pars$burnin),
                model = structure(list(), type = "SIR", residual = "on", groupeff = "on"),
                data = structure(list(), N = pars$nsires + pars$nprogeny),
                data = structure(list(), Z = pars$ngroups),
                inference = structure(list(), tmin = 0, tmax = tmax),
                observation = structure(list(), tmin = 0, tmax = tmax), # tmax = 100.01
                prior = structure(list(), parameter = "β",    type = "Flat", val1 = 0,    val2 = 1),
                prior = structure(list(), parameter = "γ",    type = "Flat", val1 = 0,    val2 = 1),
                prior = structure(list(), parameter = "k",    type = "Flat", val1 = 1,    val2 = 5),
                # prior = structure(list(), parameter = "k",    type = "Flat", val1 = 0.99, val2 = 1.01),
                prior = structure(list(), parameter = "G",    type = "Flat", val1 = -2.3, val2 = 2.3),
                prior = structure(list(), parameter = "σ_G",  type = "Flat", val1 = 0,    val2 = 0.5),
                prior = structure(list(), parameter = "q_g",  type = "Flat", val1 = -5,   val2 = 5),
                prior = structure(list(), parameter = "q_f",  type = "Flat", val1 = -5,   val2 = 5),
                prior = structure(list(), parameter = "q_r",  type = "Flat", val1 = -5,   val2 = 5),
                prior = structure(list(), parameter = "Ω_gg", type = "Flat", val1 = 0.01, val2 = 3),
                prior = structure(list(), parameter = "Ω_ff", type = "Flat", val1 = 0.01, val2 = 3),
                prior = structure(list(), parameter = "Ω_gf", type = "Flat", val1 = -3,   val2 = 3),
                prior = structure(list(), parameter = "Ω_rr", type = "Flat", val1 = 0.01, val2 = 3),
                prior = structure(list(), parameter = "Ω_gr", type = "Flat", val1 = -3,   val2 = 3),
                prior = structure(list(), parameter = "Ω_fr", type = "Flat", val1 = -3,   val2 = 3),
                prior = structure(list(), parameter = "ε_g",  type = "Flat", val1 = -5,   val2 = 5),
                prior = structure(list(), parameter = "ε_f",  type = "Flat", val1 = -5,   val2 = 5),
                prior = structure(list(), parameter = "ε_r",  type = "Flat", val1 = -5,   val2 = 5),
                prior = structure(list(), parameter = "Ψ_gg", type = "Flat", val1 = 0.01, val2 = 3),
                prior = structure(list(), parameter = "Ψ_ff", type = "Flat", val1 = 0.01, val2 = 3),
                prior = structure(list(), parameter = "Ψ_gf", type = "Flat", val1 = -3,   val2 = 3),
                prior = structure(list(), parameter = "Ψ_rr", type = "Flat", val1 = 0.01, val2 = 3),
                prior = structure(list(), parameter = "Ψ_gr", type = "Flat", val1 = -3,   val2 = 3),
                prior = structure(list(), parameter = "Ψ_fr", type = "Flat", val1 = -3,   val2 = 3),

                datatable = dt_generated_xml, # generated from Richard
                # datatable = dt_chris_xml, # from Chris's example.xml
                # Ainv_nonzero = Ainv_chris_xml,
                A_nonzero = A_nonzero_xml,
                prediction_accuracy = pa_sire_xml,
                prediction_accuracy = pa_progeny_xml
            )
        )
    )


    xml_x <- as_xml_document(x)
    write_xml(xml_x, "foo.xml", options = c("as_xml", "format"))
}
