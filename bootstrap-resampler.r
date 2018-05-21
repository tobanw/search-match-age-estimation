n.samples <- 100
bootstrap.dir <- "data/bootstrap-samples"
resamp.prefix <- "resamp_"

message('Creating ', n.samples, ' bootstrap samples in ', file.path(bootstrap.dir, resamp.prefix))

library(DBI) # RSQLite database functions
library(data.table)

set.seed(321)

# connect to sqlite database
# table names: acs, mig2met
db <- dbConnect(RSQLite::SQLite(), 'data/acs_08-16.db')

top_msa <- '(35620, 31080, 16980, 19100, 37980, 26420, 47900, 33100, 12060, 14460, 41860, 19820, 38060, 40140, 42660, 33460, 41740, 45300, 41180, 12580)'

estim_cols <- '"MET2013", "YEAR", "SERIAL", "HHWT", "PERWT", "MARST", "MARRINYR", "DIVINYR",
               "GQ", "MIGPUMA1", "MIGPLAC1", "MIGRATE1D",
               "SEX", "AGE", "AGE_SP", "RACED", "RACED_SP", "EDUC", "EDUC_SP"'

# (YEAR,SERIAL) uniquely identify households (DATANUM is always 1 in my sample)

# base table
qry_base <- paste('select', estim_cols, 'from acs where "MET2013" in', top_msa)

base.dt <- data.table(dbGetQuery(db, qry_base))

# list of unique SERIAL to resample households
households <- unique(base.dt[, .(YEAR, MET2013, SERIAL)])

#' resample unique keys and join onto base dt to get a resampled dataset
resampler <- function(keys, dt) {
  # stratified resample within each (YEAR, MSA) combo
  new.keys <- keys[, .SD[sample(.N, size = .N, replace = TRUE)], .(YEAR, MET2013)]
  new.dt <- merge(new.keys, dt, by = c("YEAR", "MET2013", "SERIAL"))
  return(new.dt)
}

for (i in 1:n.samples) {
  message('Creating bootstrap sample ', i, '...')

  dir.create(path = file.path(bootstrap.dir, paste0(resamp.prefix, i)),
             showWarnings = FALSE)

  fwrite(resampler(households, base.dt),
         file = file.path(bootstrap.dir, paste0(resamp.prefix, i), "acs_08-16.csv"))
}

message('DONE!')
