#Population masses: total, singles, couples.
#Since these are purely from the data, they only need to be computed once and saved.

using DataFrames, Query

### Array Conversion ###

"Construct and fill array for individual values in one marriage	market."
function indiv_array(val_dt::DataFrame)
    # instantiate array of zeros to fill: age, edu, race
	valarray = zeros(n_ages+1,2,2) # include age 25 for inflows

    # fill counts
	rowidx = Array{Int64}(3) # initialize empty vector for reuse
    for i in 1:size(val_dt,1) # loop over rows of table
		# map values to array indices
		# age: map to 1:n_ages (edu, race already 1:2)
		rowidx[:] = squeeze(Array(val_dt[i,1:end-1]), 1) + [1-min_age,0,0]
		valarray[rowidx...] = val_dt[i,end] # last dim is counts
    end

    return valarray
end # indiv_array
 
"Construct and fill array for couple masses in one marriage	market."
function marr_array(counts::DataFrame)
	# assumes counts has 3+3 type cols plus masses (no MSA or other cols)
    # instantiate array of zeros to fill
	massarray = zeros(n_ages+1,2,2,n_ages+1,2,2) # include age 25 for inflows

    # fill counts
	rowidx = Array{Int64}(6) # initialize empty vector for reuse
    for i in 1:size(counts,1) # loop over rows of table
		# map values to array indices
		# age: map to 1:n_ages (edu, race already 1:2)
		rowidx[:] = squeeze(Array(counts[i,1:end-1]), 1) + [1-min_age,0,0,1-min_age,0,0]
		massarray[rowidx...] = counts[i,end] # last dim is counts
    end

    return massarray
end # marr_array

"Outer operation function for population measures and ψ."
function outer_op(op::Function, men::Array, wom::Array)
	out = Array{Float64}(size(men)...,size(wom)...) # assumes 3+3 dims
	for xy in CartesianRange(size(out))
		x = xy.I[1:3]
		y = xy.I[4:6]
		out[xy] = op(men[x...], wom[y...]) # apply operator
	end
	return out
end

"Load data, aggregating out unwanted types."
function loader(a)
	if length(size(a)) == 3
		if no_edu && no_race
			out = sum(a, (2,3))
		elseif no_edu
			out = sum(a, 2)
		elseif no_race
			out = sum(a, 3)
		else
			out = a
		end
	elseif length(size(a)) == 6
		if no_edu && no_race
			out = sum(a, (2,3,5,6))
		elseif no_edu
			out = sum(a, (2,5))
		elseif no_race
			out = sum(a, (3,6))
		else
			out = a
		end
	else
		out = a
	end
	return out
end

# Set `data_dir` in main-estim.jl

# load up smoothed masses from csv
df_pop = readtable(joinpath(data_dir, "pop.csv"))
df_marriages = readtable(joinpath(data_dir, "marriages.csv"))
df_migration = readtable(joinpath(data_dir, "mar-migration.csv"))

# save arrays (by MSA) in separate dicts, to be stored
men_sng = Dict{AbstractString, Array}()
wom_sng = Dict{AbstractString, Array}()
men_tot = Dict{AbstractString, Array}()
wom_tot = Dict{AbstractString, Array}()
marriages = Dict{AbstractString, Array}()
mar_outflow = Dict{AbstractString, Array}()
sng_outer = Dict{AbstractString, Array}() # u_m(x)*u_f(y) arrays
pop_outer = Dict{AbstractString, Array}() # ℓ_m(x)*ℓ_f(y) arrays

# load up flows for rate parameter estimation
df_MF = readtable(joinpath(data_dir, "pair-MF.csv"))
df_DF = readtable(joinpath(data_dir, "ind-DF.csv"))

MF = Dict{AbstractString, Array}()
men_DF = Dict{AbstractString, Array}()
wom_DF = Dict{AbstractString, Array}()

for msa in top_msa
	# annual population arrays (NOTE: includes age 25)
	men_sng["$msa"] = loader(indiv_array(@from i in df_pop begin
											 @where i.MSA == msa && i.SEX == 1
											 @select {i.AGE, i.COLLEGE, i.MINORITY, i.SNG}
											 @collect DataFrame
										 end))

	wom_sng["$msa"] = loader(indiv_array(@from i in df_pop begin
											 @where i.MSA == msa && i.SEX == 2
											 @select {i.AGE, i.COLLEGE, i.MINORITY, i.SNG}
											 @collect DataFrame
										 end))

	men_tot["$msa"] = loader(indiv_array(@from i in df_pop begin
											 @where i.MSA == msa && i.SEX == 1
											 @select {i.AGE, i.COLLEGE, i.MINORITY, i.POP}
											 @collect DataFrame
										 end))

	wom_tot["$msa"] = loader(indiv_array(@from i in df_pop begin
											 @where i.MSA == msa && i.SEX == 2
											 @select {i.AGE, i.COLLEGE, i.MINORITY, i.POP}
											 @collect DataFrame
										 end))

	marriages["$msa"] = loader(marr_array(@from i in df_marriages begin
											  @where i.MSA == msa
											  @select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.MASS}
											  @collect DataFrame
										  end))

	mar_outflow["$msa"] = loader(marr_array(@from i in df_migration begin
												@where i.MSA == msa
												@select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.NET_OUTFLOW}
												@collect DataFrame
										  end))

	sng_outer["$msa"] = outer_op(*, men_sng["$msa"], wom_sng["$msa"])
	pop_outer["$msa"] = outer_op(*, men_tot["$msa"], wom_tot["$msa"])

	# annual flow arrays
	MF["$msa"] = loader(marr_array(@from i in df_MF begin
									   @where i.MSA == msa
									   @select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.FLOW}
									   @collect DataFrame
								   end))

	men_DF["$msa"] = loader(indiv_array(@from i in df_DF begin
											@where i.MSA == msa && i.SEX == 1
											@select {i.AGE, i.COLLEGE, i.MINORITY, i.FLOW}
											@collect DataFrame
										end))

	wom_DF["$msa"] = loader(indiv_array(@from i in df_DF begin
											@where i.MSA == msa && i.SEX == 2
											@select {i.AGE, i.COLLEGE, i.MINORITY, i.FLOW}
											@collect DataFrame
										end))
end

# load death arrival rates
ψ_m = loader(indiv_array(readtable(joinpath(data_dir, "men-psi.csv"))))
ψ_f = loader(indiv_array(readtable(joinpath(data_dir, "wom-psi.csv"))))
ψm_ψf = outer_op(+, ψ_m, ψ_f) # array of ψ_m(x) + ψ_f(y)

# store in JLD format
# Set `pop_file` in main-estim.jl
jldopen(pop_file, "w") do file  # open file for saving julia data
	# arrays: death arrival rates (includes age 25)
	write(file, "men_psi", ψ_m)
	write(file, "wom_psi", ψ_f)
	write(file, "psi_outer", ψm_ψf)

	# dictionaries: stocks per MSA
    write(file, "men_sng", men_sng)
    write(file, "wom_sng", wom_sng)
    write(file, "men_tot", men_tot)
    write(file, "wom_tot", wom_tot)
    write(file, "marriages", marriages)
    write(file, "mar_outflow", mar_outflow)
    write(file, "sng_outer", sng_outer)
    write(file, "pop_outer", pop_outer)

	# dictionaries: flows per MSA
    write(file, "MF", MF)
    write(file, "men_DF", men_DF)
    write(file, "wom_DF", wom_DF)
end
