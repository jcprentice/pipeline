source("model_SIS.R")
source("model_SIR.R")
source("model_SEIR.R")

simulate_epidemic <- function(pars, traitdata) {
    switch(pars$epi_type,
           SIS  = return(model_SIS(pars, traitdata)),
           SIR  = return(model_SIR(pars, traitdata)),
           SEIR = return(model_SEIR(pars, traitdata))
    )
}

