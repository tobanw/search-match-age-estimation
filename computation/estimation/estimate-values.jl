# Use smoothed population masses and estimated arrival rates to recover marital production function by MSA

using DataTables

# TODO
#	* construct arrays by MSA
#	* helper functions: copy from MarriageMarkets, or just add and import?

# NOTE: must match max.age from `smooth-pop.r`
max_age = 65 # terminal age (inclusive)
min_age = 25 # initial age (inclusive)
n_ages = max_age - min_age + 1 # inclusive of both endpoints

"Construct and fill array for individual masses"
function indiv_array(counts::DataTable)
    # instantiate empty array: age, edu, race
	massarray = Array(Float64, (n_ages,2,2))

    # fill counts
	row = Array(Float64, size(counts,2)) # initialize empty vector for reuse
    for i in 1:size(counts,1) # loop over rows of table
		# map values to array indices
        # rac, edu: map 0:1 to 1:2
		# age: map to 1:n_ages

		#TODO: need Int for indices!
		row[:] = squeeze(Array(counts[i,:]), 1) + [1 - min_age, 1, 1, 0]
		massarray[row[1:end-1]...] = row[end] # last dim is counts
    end

    return massarray

end # masses_array
 
"Construct and fill array for couple masses"
function marr_array(counts::DataTable)
    # instantiate empty array
	massarray = Array(Float64, (n_ages,2,2,n_ages,2,2))

    # fill counts
	row = Array(Float64, size(counts,2)) # initialize empty vector for reuse
    for i in 1:size(counts,1) # loop over rows of table
		# map values to array indices
        # rac, edu: map 0:1 to 1:2
		# age: map to 1:n_ages
		row[:] = squeeze(Array(counts[i,:]), 1) + [1 - min_age, 1, 1, 1 - min_age, 1, 1,0]
		massarray[row[1:end-1]...] = row[end] # last dim is counts
    end

    return massarray

end # masses_array
