#Population masses: total, singles, couples.
#Since these are purely from the data, they only need to be computed once and saved.

using DataFrames, Query

### Array Conversion ###

"Construct and fill array for individual values in one marriage	market."
function indiv_array(val_dt::DataFrame)
    # instantiate empty array: age, edu, race
	valarray = Array{Float64}(n_ages+1,2,2) # include age 25 for inflows

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
    # instantiate empty array
	massarray = Array{Float64}(n_ages+1,2,2,n_ages+1,2,2) # include age 25 for inflows

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

"Convolution / outer product function for population measures and ψ."
function convolution(op::Function, men::Array, wom::Array)
	convy = Array{Float64}(size(men)...,size(wom)...) # assumes 3+3 dims
	for xy in CartesianRange(size(convy))
		x = xy.I[1:3]
		y = xy.I[4:6]
		convy[xy] = op(men[x...], wom[y...]) # apply operator
	end
	return convy
end

# load up smoothed masses from csv
df_men_sng = readtable("data/men-single.csv")
df_wom_sng = readtable("data/wom-single.csv")
df_men_tot = readtable("data/men-total.csv")
df_wom_tot = readtable("data/wom-total.csv")
df_marriages = readtable("data/marriages.csv")

# save arrays (by MSA) in separate dicts, to be stored
men_sng = Dict{AbstractString, Array}()
wom_sng = Dict{AbstractString, Array}()
men_tot = Dict{AbstractString, Array}()
wom_tot = Dict{AbstractString, Array}()
marriages = Dict{AbstractString, Array}()
sng_conv = Dict{AbstractString, Array}() # u_m(x)*u_f(y) arrays
pop_conv = Dict{AbstractString, Array}() # ℓ_m(x)*ℓ_f(y) arrays

# load up flows for rate parameter estimation
df_MF = readtable("data/MF.csv")
df_men_DF = readtable("data/men-DF.csv")
df_wom_DF = readtable("data/wom-DF.csv")

MF = Dict{AbstractString, Array}()
men_DF = Dict{AbstractString, Array}()
wom_DF = Dict{AbstractString, Array}()

for msa in top_msa
	# annual population arrays (NOTE: includes age 25)
	men_sng["$msa"] = indiv_array(@from i in df_men_sng begin
                                  @where i.MSA == msa
                                  @select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.MASS}
                                  @collect DataFrame
                                  end) / n_years

	wom_sng["$msa"] = indiv_array(@from i in df_wom_sng begin
                                  @where i.MSA == msa
                                  @select {i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.MASS}
                                  @collect DataFrame
                                  end) / n_years

	men_tot["$msa"] = indiv_array(@from i in df_men_tot begin
                                  @where i.MSA == msa
                                  @select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.MASS}
                                  @collect DataFrame
                                  end) / n_years

	wom_tot["$msa"] = indiv_array(@from i in df_wom_tot begin
                                  @where i.MSA == msa
                                  @select {i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.MASS}
                                  @collect DataFrame
                                  end) / n_years

	marriages["$msa"] = marr_array(@from i in df_marriages begin
                                   @where i.MSA == msa
                                   @select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.MASS}
                                   @collect DataFrame
                                   end) / n_years

	sng_conv["$msa"] = convolution(*, men_sng["$msa"], wom_sng["$msa"])
	pop_conv["$msa"] = convolution(*, men_tot["$msa"], wom_tot["$msa"])

	# annual flow arrays
	MF["$msa"] = marr_array(@from i in df_MF begin
							@where i.MSA == msa
							@select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.FLOW}
							@collect DataFrame
							end) / n_years

	men_DF["$msa"] = indiv_array(@from i in df_men_DF begin
                                 @where i.MSA == msa
                                 @select {i.AGE, i.COLLEGE, i.MINORITY, i.FLOW}
                                 @collect DataFrame
                                 end) / n_years

	wom_DF["$msa"] = indiv_array(@from i in df_wom_DF begin
                                 @where i.MSA == msa
                                 @select {i.AGE, i.COLLEGE, i.MINORITY, i.FLOW}
                                 @collect DataFrame
                                 end) / n_years
end

# load death arrival rates
ψ_m = indiv_array(readtable("results/men-psi.csv"))
ψ_f = indiv_array(readtable("results/wom-psi.csv"))
ψm_ψf = convolution(+, ψ_m, ψ_f) # array of ψ_m(x) + ψ_f(y)

# store in JLD format
jldopen("results/populations.jld", "w") do file  # open file for saving julia data
	# arrays: death arrival rates (includes age 25)
	write(file, "men_psi", ψ_m)
	write(file, "wom_psi", ψ_f)
	write(file, "psi_conv", ψm_ψf)

	# dictionaries: stocks per MSA
    write(file, "men_sng", men_sng)
    write(file, "wom_sng", wom_sng)
    write(file, "men_tot", men_tot)
    write(file, "wom_tot", wom_tot)
    write(file, "marriages", marriages)
    write(file, "sng_conv", sng_conv)
    write(file, "pop_conv", pop_conv)

	# dictionaries: flows per MSA
    write(file, "MF", MF)
    write(file, "men_DF", men_DF)
    write(file, "wom_DF", wom_DF)
end
