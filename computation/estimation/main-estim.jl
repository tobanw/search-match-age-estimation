using JLD, NLopt

### USAGE OPTIONS ###

reload_smooth = false # set true to reload the smoothed population csv files, otherwise use saved JLD
estimate_rates = true # set true to run SMM estimation, otherwise just use λ_0, δ_0 below


### PARAMS ### 

const r = 0.04 # discount rate
const ρ = 1.0 # aging rate
const β = 0.5 # wife share of surplus

# NOTE: must match `max.age` from `smooth-pop.r` and mortality data
const max_age = 65 # terminal age (inclusive)
const min_age = 25 # initial age (excluded)
const n_ages = max_age - min_age # excluding 25
const n_years = 7 # 2008-2014

# initial guess of arrival rate parameters for GMM routine
# parametric specification:
# ξ^k =  exp(-dot(ζ, [1,k,k^2,...]))
# results: total loss
#	2 params: 7.771e12
#	3 params: ?.???e12
#ζ_0 = [0.364, 0.159, 0.048, 0.093]  # fit for linear interpolation on [0,4,12,39]
#ζ_0 = [0.518, 0.368, 0.166, 0.021, 0.0052]  # fit for linear interpolation on [0,4,10,18,39]
#ζ_0 = [0.524, 0.373, 0.166, 0.022, 0.00153, 0.1] # interpolate on [0,4,10,18,27,39]
#ζ_0 = [0.527, 0.411, 0.266, 0.121, 0.0242, 0.0013, 0.1] # interpolate on [0,3,7,12,18,27,39]

# TODO: extend to k+1 knots, estimate, repeat
znodes = [0,3,7,12,18,27,39]
ζ_0 = [0.527, 0.411, 0.266, 0.121, 0.0242, 0.0013, 0.1]

δ_0 = 0.06 # arrival rate of love shocks

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


### Arrival Rates: GMM estimation ###


if estimate_rates # run optimizer for MD estimation
	println("Starting optimization routine for GMM estimation...")

	# compute moments in parallel: multiprocess
	#addprocs(4) # add a worker process per moment (male/female marriage/divorce flows)
	#println("  > Using $(nprocs()-1) worker processes")

	opt = Opt(:LN_SBPLX, length(ζ_0)+1) # number of parameters
	lower_bounds!(opt, zeros(length(ζ_0)+1)) # enforce non-negativity
	#upper_bounds!(opt, ones(length(ζ_0)+1)) # prevent it from roaming
	upper_bounds!(opt, [1.5, 1.3, 1.2, 0.4, 0.25, 0.2, 0.1, 0.2]) # prevent it from roaming
	xtol_rel!(opt, 1e-3) # tolerance
	min_objective!(opt, loss_opt) # specify objective function
	min_f, min_x, ret = optimize(opt, [ζ_0..., δ_0]) # run!

	# estimates
	ζ = min_x[1:end-1]
	δ = min_x[end]

	println("Done! Linear interpolation of ξ with raw α. Estimates:")
	println("ζ: $ζ") 
	println("δ: $δ")
else # just use the initial guess
	ζ = ζ_0
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
alpha_raw = Dict{AbstractString, Array}()
surplus = Dict{AbstractString, Array}()
men_val = Dict{AbstractString, Array}()
wom_val = Dict{AbstractString, Array}()
production = Dict{AbstractString, Array}()

for msa in top_msa
	# U_m, U_f per marriage market
	ξ = build_ξ(ζ, znodes) # construct ξ
	λ = inflate_λ(ξ ./ sqrt(sum(men_sng["$msa"]) * sum(wom_sng["$msa"])))

    # match probability (α)
	alpha["$msa"] = compute_alpha(λ, δ, ψ_m, ψ_f, marriages["$msa"], men_sng["$msa"], wom_sng["$msa"])
	# TESTING
	alpha_raw["$msa"] = compute_raw_alpha(λ, δ, ψ_m, ψ_f, marriages["$msa"], men_sng["$msa"], wom_sng["$msa"])

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
jldopen("results/estimates-spline.jld", "w") do file  # open file for saving julia data
    write(file, "zeta", ζ)
    write(file, "delta", δ)
    write(file, "alpha", alpha)
    write(file, "alpha_raw", alpha_raw)
    write(file, "surplus", surplus)
    write(file, "men_val", men_val)
    write(file, "wom_val", wom_val)
    write(file, "production", production)
    write(file, "avg_production", avg_production)
end # do
