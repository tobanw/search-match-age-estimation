# Use smoothed population masses and estimated arrival rates to recover marital production function by MSA
# Treat initial age as birth inflow (mass * phi), use u[2:end], n[2:end,2:end] for alpha, s, f

using DataTables, Query, Distributions, JLD

# TODO:
#	* optimization: declare params as constants

### Array Conversion ###

"Construct and fill array for individual masses in one marriage	market."
function indiv_array(val_dt::DataTable)
    # instantiate empty array: age, edu, race
	valarray = Array(Float64, (n_ages,2,2))

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
	massarray = Array(Float64, (n_ages,2,2,n_ages,2,2))

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


### PARAMS ### 

# arrival rates from `results/rate-param.csv`
λ = 6.2e-9
δ = 0.0117

ρ = 1.0 # aging rate

# death arrival rates: includes age 25
ψ_m = indiv_array(readtable("results/men-psi.csv")) # (AGE, COLLEGE, MINORITY, PSI)
ψ_f = indiv_array(readtable("results/wom-psi.csv")) # (AGE, COLLEGE, MINORITY, PSI)
 
# inverse cdf of standard normal distribution (for inverting α)
const STDNORMAL = Normal()
Φ_inv(x::Real) = quantile(STDNORMAL, x)

# NOTE: must match `max.age` from `smooth-pop.r`
max_age = 65 # terminal age (inclusive)
min_age = 25 # initial age (inclusive)
n_ages = max_age - min_age + 1 # inclusive of both endpoints

# NOTE: must match `top.msa` from `smooth-pop.r`
top_msa = (35620, 31080, 16980, 19100, 37980, 26420, 47900, 33100, 12060, 14460)


### Estimation Functions ###

"Compute α for a given MSA"
function compute_alpha(mar_init::Array, um_init::Array, uf_init::Array)
	# trim off initial ages
	m = mar_init[2:end,:,:,2:end,:,:] # mar_init[i] == mar[i-1] along the age dims
	u_m = um_init[2:end,:,:]
	u_f = uf_init[2:end,:,:]

	α = similar(m) # instantiate alpha array
	for xy in CartesianRange(m)
		x = xy.I[1:3]
		y = xy.I[4:6]
		α[xy] = (m[xy] * (ρ + δ + ψ_m[x] + ψ_f[y] ) - ρ * mar_init[xy]) /
		        (λ * u_m[x] * u_f[y] + δ * m[xy])
	end
	return clamp.(α, 1e-6, 1 - 1e-6) # enforce 0 < α < 1
end

"μ function, using inverse Mills ratio: E[z|z>q] = σ ϕ(q/σ) / (1 - Φ(q/σ))."
function μ(a::Real)
	st = Φ_inv(1-a) # pre-compute -s/σ = Φ^{-1}(1-a)
	return pdf(STDNORMAL, st) - a * st
end

# TODO: c, d


### Perform Estimation ###

# load up smoothed masses from csv
dt_men_sng = readtable("data/men-single.csv")
dt_wom_sng = readtable("data/wom-single.csv")
dt_men_tot = readtable("data/men-total.csv")
dt_wom_tot = readtable("data/wom-total.csv")
dt_marriages = readtable("data/marriages.csv")

# save arrays of each type in separate dicts, to be stored
men_sng = Dict{AbstractString, Array}()
wom_sng = Dict{AbstractString, Array}()
men_tot = Dict{AbstractString, Array}()
wom_tot = Dict{AbstractString, Array}()
marriages = Dict{AbstractString, Array}()
alpha = Dict{AbstractString, Array}()
surplus = Dict{AbstractString, Array}()
production = Dict{AbstractString, Array}()

for msa in top_msa

    # population arrays
	men_sng["$msa"] = indiv_array(@from i in dt_men_sng begin
                                  @where i.MSA == msa
                                  @select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.MASS}
                                  @collect DataTable
                                  end)

	wom_sng["$msa"] = indiv_array(@from i in dt_wom_sng begin
                                  @where i.MSA == msa
                                  @select {i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.MASS}
                                  @collect DataTable
                                  end)

	men_tot["$msa"] = indiv_array(@from i in dt_men_tot begin
                                  @where i.MSA == msa
                                  @select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.MASS}
                                  @collect DataTable
                                  end)

	wom_tot["$msa"] = indiv_array(@from i in dt_wom_tot begin
                                  @where i.MSA == msa
                                  @select {i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.MASS}
                                  @collect DataTable
                                  end)

	marriages["$msa"] = marr_array(@from i in dt_marriages begin
                                   @where i.MSA == msa
                                   @select {i.AGE_M, i.COLLEGE_M, i.MINORITY_M, i.AGE_F, i.COLLEGE_F, i.MINORITY_F, i.MASS}
                                   @collect DataTable
                                   end)

    # match probability (α)
	alpha["$msa"] = compute_alpha(marriages["$msa"], men_sng["$msa"], wom_sng["$msa"])

    # marital surplus (S)
	surplus["$msa"] = Φ_inv.(alpha["$msa"])

    # marital production (f)
    #production["$msa"] = 

end # for


# store in JLD format
jldopen("results/estimates.jld", "w") do file  # open file for saving julia data
    write(file, "men_sng", men_sng)
    write(file, "wom_sng", wom_sng)
    write(file, "men_tot", men_tot)
    write(file, "wom_tot", wom_tot)
    write(file, "marriages", marriages)
    write(file, "surplus", surplus)
    write(file, "production", production)
end # do
