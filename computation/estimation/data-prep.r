library(data.table)

# Build a table to match a (MIGPLAC1, MIGPUMA1) to a single MSA
#	* Problem 1: MIGPUMA1 is a larger unit than PUMA, can contain multiple
#		* Solution: choose any contained PUMA as (unique) foreign key
#	* Problem 2: PUMA doesn't always stay within MET2013 boundaries, can belong to multiple MSAs
#		* Solution: choose any associated MSA as (unique) foreign key

# load data
mig.map <- fread('data/PUMA_MIG_conversion.csv')
met.map <- fread('data/MET_PUMA_conversion.csv')

# Use tail or head to select a single row per group:
# choose a unique puma per migpuma
mig.unq <- mig.map[, .(PUMA = tail(PUMA, 1)), by = .(MIGPLAC1, MIGPUMA1, STATE)]
# choose a unique msa per puma
met.unq <- met.map[, .(MSA = tail(MSA, 1)), by = .(STATE, PUMA)]

# create mapping
mig2met.map <- merge(mig.unq, met.unq, by = c("STATE", "PUMA"))

# write to csv (to be loaded into sqlite db with `csv2sqlite.py`)
fwrite(mig2met.map[, .(MIGPLAC1, MIGPUMA1, MSA)], file = "mig2met.csv")
