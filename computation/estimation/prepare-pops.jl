#Population masses: total, singles, couples.
#Since these are purely from the data, they only need to be computed once and saved.

using DataTables, Query

### Array Conversion ###

"Construct and fill array for individual values in one marriage	market."
function indiv_array(val_dt::DataTable)
    # instantiate empty array: age, edu, race
	valarray = Array(Float64, (n_ages+1,2,2)) # include age 25 for inflows

    # fill counts
	rowidx = Array(Int64, 3) # initialize empty vector for reuse
    for i in 1:size(val_dt,1) # loop over rows of table
		# map values to array indices
        # rac, edu: map 0:1 to 1:2
		# age: map to 1:n_ages
		rowidx[:] = squeeze(Array(val_dt[i,1:end-1]), 1) + [1 - min_age, 1, 1]
		valarray[rowidx...] = get(val_dt[i,end]) # last dim is counts
    end

    return valarray
end # indiv_array
 
"Construct and fill array for couple masses in one marriage	market."
function marr_array(counts::DataTable)
	# assumes counts has 3+3 type cols plus masses (no MSA or other cols)
    # instantiate empty array
	massarray = Array(Float64, (n_ages+1,2,2,n_ages+1,2,2)) # include age 25 for inflows

    # fill counts
	rowidx = Array(Int64, 6) # initialize empty vector for reuse
    for i in 1:size(counts,1) # loop over rows of table
		# map values to array indices
        # rac, edu: map 0:1 to 1:2
		# age: map to 1:n_ages
		rowidx[:] = squeeze(Array(counts[i,1:end-1]), 1) + [1 - min_age, 1, 1, 1 - min_age, 1, 1]
		massarray[rowidx...] = get(counts[i,end]) # last dim is counts
    end

    return massarray
end # marr_array


# load up smoothed masses from csv
dt_men_sng = readtable("data/men-single.csv")
dt_wom_sng = readtable("data/wom-single.csv")
dt_men_tot = readtable("data/men-total.csv")
dt_wom_tot = readtable("data/wom-total.csv")
dt_marriages = readtable("data/marriages.csv")

# save arrays (by MSA) in separate dicts, to be stored
men_sng = Dict{AbstractString, Array}()
wom_sng = Dict{AbstractString, Array}()
men_tot = Dict{AbstractString, Array}()
wom_tot = Dict{AbstractString, Array}()
marriages = Dict{AbstractString, Array}()

for msa in top_msa
	# annual population arrays (NOTE: includes age 25)
	men_sng["$msa"] = indiv_array(@from i in dt_men_sng begin
                                  @where i.MSA == msa
                                  @select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.MASS}
                                  @collect DataTable
                                  end) / n_years

	wom_sng["$msa"] = indiv_array(@from i in dt_wom_sng begin
                                  @where i.MSA == msa
                                  @select {i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.MASS}
                                  @collect DataTable
                                  end) / n_years

	men_tot["$msa"] = indiv_array(@from i in dt_men_tot begin
                                  @where i.MSA == msa
                                  @select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.MASS}
                                  @collect DataTable
                                  end) / n_years

	wom_tot["$msa"] = indiv_array(@from i in dt_wom_tot begin
                                  @where i.MSA == msa
                                  @select {i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.MASS}
                                  @collect DataTable
                                  end) / n_years

	marriages["$msa"] = marr_array(@from i in dt_marriages begin
                                   @where i.MSA == msa
                                   @select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.MASS}
                                   @collect DataTable
                                   end) / n_years
end

# store in JLD format
jldopen("results/populations.jld", "w") do file  # open file for saving julia data
	# arrays: death arrival rates (includes age 25)
	write(file, "men_psi", indiv_array(readtable("results/men-psi.csv")))
	write(file, "wom_psi", indiv_array(readtable("results/wom-psi.csv")))

	# dictionaries: per MSA
    write(file, "men_sng", men_sng)
    write(file, "wom_sng", wom_sng)
    write(file, "men_tot", men_tot)
    write(file, "wom_tot", wom_tot)
    write(file, "marriages", marriages)
end
