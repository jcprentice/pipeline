# load the data into a DataTable
fishboost_data <- as.data.table(read.csv("../fishboost/Fishboost Masterfile V4_donors+sex.csv"))

names(fishboost_data)

# select only the columns we actually need
fb <- fishboost_data[, .(
    Individual.recoded,
    Trial,
    Donor..D....Receptor..R.,
    Sire_rec,
    Dam_rec,
    Day.of.first.symptoms,
    Day.of.death,
    Day.start.to.end..resilience.,
    Days.onset.symptoms..susceptibility.,
    Days.onset.to.death..tolerance.
)]

setnames(
    fb,
    old = c(
        "Individual.recoded",
        "Trial",
        "Donor..D....Receptor..R.",
        "Sire_rec",
        "Dam_rec",
        "Day.of.first.symptoms",
        "Day.of.death",
        "Day.start.to.end..resilience.",
        "Days.onset.symptoms..susceptibility.",
        "Days.onset.to.death..tolerance."
    ), 
    new = c("id", "trial", "dr", "sire", "dam", "It", "Rt", "res", "sus", "tol")
)

fb[, It := as.Date(It, format = "%Y-%m-%d")]
fb[, Rt := as.Date(Rt, format = "%Y-%m-%d")]

start <- as.Date(min(fb[, It], na.rm=T))
fb[, It2 := It - start]
