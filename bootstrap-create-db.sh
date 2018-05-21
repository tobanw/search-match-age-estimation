parallel python data/csv2sqlite.py {}/acs_08-16.csv {}/acs_08-16.db acs ::: data/bootstrap-samples/resamp_*
parallel python data/csv2sqlite.py data/mig2met.csv {}/acs_08-16.db mig2met ::: data/bootstrap-samples/resamp_*
