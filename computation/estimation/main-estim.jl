using JLD, NLopt

### USAGE OPTIONS ###

reload_smooth = false # set true to reload the smoothed population csv files, otherwise use saved JLD
estimate_rates = false # set true to run SMM estimation, otherwise just use λ_0, δ_0 below


### PARAMS ### 

# initial guess of arrival rates for SMM routine
λ_0 = 2.023848578e-7
δ_0 = 0.0266481127

const r = 0.04 # discount rate
const ρ = 1.0 # aging rate
const β = 0.5 # wife share of surplus

# NOTE: must match `max.age` from `smooth-pop.r` and mortality data
const max_age = 65 # terminal age (inclusive)
const min_age = 25 # initial age (excluded)
const n_ages = max_age - min_age # excluding 25
const n_years = 7 # 2008-2014

# NOTE: must match `top.msa` from `smooth-pop.r` and `smooth-flows.r`
const top_msa = (35620, 31080, 16980, 19100, 37980, 26420, 47900, 33100, 12060, 14460, 41860, 19820, 38060, 40140, 42660, 33460, 41740, 45300, 41180, 12580)
const n_msa = length(top_msa)


### DATA ###

if reload_smooth
	# load from csv into dicts (and also save as JLD)
	include("prepare-pops.jl")
else # load from existing JLD file
	# load saved data as dict
	pp = load("results/populations.jld")

	# load and trim age 25 from psi
	ψ_m = pp["men_psi"][2:end,:,:] # ages 26-65
	ψ_f = pp["wom_psi"][2:end,:,:]

	# load population stock dicts (ages 25-65 inclusive)
	men_sng = pp["men_sng"]
	wom_sng = pp["wom_sng"]
	men_tot = pp["men_tot"]
	wom_tot = pp["wom_tot"]
	marriages = pp["marriages"]

	# load data	moments: marriage and divorce flow dicts (ages 25-65 inclusive)
	men_MF = pp["men_MF"]
	men_DF = pp["men_DF"]
	wom_MF = pp["wom_MF"]
	wom_DF = pp["wom_DF"]
end


### ESTIMATION ###

# estimation functions
include("estim-functions.jl")

### Arrival Rates: SMM estimation ###

"Objective function to minimize: distance between model and data moments."
function loss(λ::Real, δ::Real, ψ_m::Array, ψ_f::Array,
	          mar_all::Dict, um_all::Dict, uf_all::Dict,
	          dMF_m::Dict, dDF_m::Dict, dMF_f::Dict, dDF_f::Dict,
	          wgt_men::Dict, wgt_wom::Dict)

	sse = 0.0 # initialize
	for msa in top_msa
		# model_moments returns ages 26-65
		MF_m, DF_m, MF_f, DF_f = model_moments(λ, δ, ψ_m, ψ_f, mar_all["$msa"],
		                                       um_all["$msa"], uf_all["$msa"])

		# feed trimmed (age 26-64) moments and weights to loss function
		sse += loss_msa(MF_m[1:end-1,:,:], DF_m[1:end-1,:,:],
		                MF_f[1:end-1,:,:], DF_f[1:end-1,:,:],
	                    dMF_m["$msa"][2:end-1,:,:], dDF_m["$msa"][2:end-1,:,:],
		                dMF_f["$msa"][2:end-1,:,:], dDF_f["$msa"][2:end-1,:,:],
		                wgt_men["$msa"][2:end-1,:,:], wgt_wom["$msa"][2:end-1,:,:])
	end
	return sse
end

"Objective function to pass to NLopt: requires vectors for `x` and `grad`."
function loss_opt(x::Vector, grad::Vector)
	return loss(x[1], x[2], ψ_m, ψ_f, marriages, men_sng, wom_sng,
	            men_MF, men_DF, wom_MF, wom_DF,
	            men_tot, wom_tot)
end


if estimate_rates # run optimizer for MD estimation
	println("Starting optimizer for SMM estimation...")

	opt = Opt(:LN_NELDERMEAD, 2) # 2 parameters
	lower_bounds!(opt, [0,0]) # enforce non-negativity
	xtol_rel!(opt, 1e-10) # tolerance
	min_objective!(opt, loss_opt) # specify objective function
	min_f, min_x, ret = optimize(opt, [λ_0, δ_0]) # run!

	# estimates
	λ = min_x[1]
	δ = min_x[2]

	println("Done! Estimates:")
	println("λ: $λ")
	println("δ: $δ")
else # just use the initial guess
	λ = λ_0
	δ = δ_0
end


### Non-parametric objects ###

# precompute discount factors
d = compute_d(δ, ψ_m, ψ_f)
c = compute_c(d)
c1 = compute_c1(c)

avg_production = zeros(d) # initialize

# save arrays (by MSA) in separate dicts, to be stored
alpha = Dict{AbstractString, Array}()
surplus = Dict{AbstractString, Array}()
men_val = Dict{AbstractString, Array}()
wom_val = Dict{AbstractString, Array}()
production = Dict{AbstractString, Array}()

for msa in top_msa
    # match probability (α)
	alpha["$msa"] = compute_alpha(λ, δ, ψ_m, ψ_f, marriages["$msa"], men_sng["$msa"], wom_sng["$msa"])

    # marital surplus (S)
	surplus["$msa"] = invert_alpha(c1, alpha["$msa"])

	# precompute discounted μ array: solution of ∫ max{S(x,y,z),0}dG(z)
	dc1μ = compute_dc1μ(d, c1, alpha["$msa"])

	# average value functions
	men_val["$msa"], wom_val["$msa"] = compute_value_functions(λ, dc1μ, men_sng["$msa"], wom_sng["$msa"])

    # marital production (f)
	production["$msa"] = compute_production(δ, ψ_m, ψ_f, d, dc1μ, surplus["$msa"],
	                                        men_val["$msa"], wom_val["$msa"])

	# average across MSAs
	avg_production .+= production["$msa"]
end # for

# final production estimate, average of MSA estimates
avg_production = avg_production / n_msa

# store in JLD format
jldopen("results/estimates.jld", "w") do file  # open file for saving julia data
    write(file, "lambda", λ)
    write(file, "delta", δ)
    write(file, "alpha", alpha)
    write(file, "surplus", surplus)
    write(file, "men_val", men_val)
    write(file, "wom_val", wom_val)
    write(file, "production", production)
    write(file, "avg_production", avg_production)
end # do
