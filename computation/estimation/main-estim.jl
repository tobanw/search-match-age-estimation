using JLD, NLopt

### USAGE OPTIONS ###

reload_smooth = false # set true to reload the smoothed population csv files, otherwise use saved JLD
grid_search = false
estimate_rates = true # set true to run GMM estimation, otherwise just use ζ_0, δ_0 below
compute_np_obj = false


# PARALLEL: compute moments on worker processes
#addprocs(2)#(Sys.CPU_CORES) # add a worker process per core
#println("  > Using $(nprocs()-1) worker processes")


### PARAMS ### 

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


### INITIAL GUESSES ###

# parametric specification for ξ: linear interpolation of ζ

#znodes = [-39,-16,-10,-5,0,5,10,16,39]
#ζ_0 = [0.03, 0.04, 0.08, 0.13, 0.24, 0.17, 0.1, 0.04, 0.03]

# parametric specification for ξ: logistic fit of ζ

# loss at 3.40732e9 (for commented param values)
# ζ = [c,m,s] for age-gap function
ζ_0 = [3.0, 1.9, 2.71] # [4.39, 2.20, 4.96]
# θ = decay parameter for total couple age
θ_0 = 7.0 # 22.2
# arrival rate of love shocks
δ_0 = 0.05 # 0.055


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
	ψm_ψf = pp["psi_conv"][2:end,:,:,2:end,:,:]

	# load population stock dicts (ages 25-65 inclusive)
	men_sng = pp["men_sng"]
	wom_sng = pp["wom_sng"]
	men_tot = pp["men_tot"]
	wom_tot = pp["wom_tot"]
	marriages = pp["marriages"]
	sng_conv = pp["sng_conv"]
	pop_conv = pp["pop_conv"]

	# load data	moments: marriage and divorce flow dicts (ages 25-65 inclusive)
	MF = pp["MF"]
	men_DF = pp["men_DF"]
	wom_DF = pp["wom_DF"]
end


### ESTIMATION ###

# estimation functions
include("estim-functions.jl")

if grid_search # call obj_landscaper manually and write the string to file
	ζ1grid = [1.0,1.5,2.0,2.5,3.0,3.5,4.0] # 7 (global multiplicative scale-factor for λ)
	ζ2grid = [0.0,0.5,1.0,1.5] # 4 (mean for age-gap meeting slowdown)
	ζ3grid = [2,3,4,5] # 4 (spread for age-gap meeting slowdown)
	θgrid = [27,30,33,36,39] # 5 (old-age search slowdown)
	δgrid = [0.01, 0.015, 0.02, 0.03] # 4 (divorce rate)

	function obj_landscaper(z1g::Vector, z2g::Vector, z3g::Vector, thg::Vector, dg::Vector, ψm_ψf, marriages, sng_conv, MF, men_DF, wom_DF, men_tot, wom_tot, pop_conv) #ζ1grid, ζ2grid, ζ3grid, θgrid # might need as kwargs
		# given a single δ, sweep through other param grids
		res = "" # results string
		for ζ1 in z1g, ζ2 in z2g, ζ3 in z3g, θ in thg, δ in dg
			val = loss([ζ1,ζ2,ζ3], θ, δ, ψm_ψf, marriages, sng_conv, MF, men_DF, wom_DF, men_tot, wom_tot, pop_conv)
			res *= string(ζ1, ",", ζ2, ",", ζ3, ",", θ, ",", δ, ",", val, "\n")
		end
		return res
	end

end # grid_search

### Arrival Rates: GMM estimation ###


if estimate_rates # run optimizer for MD estimation
	println("Starting optimization routine for GMM estimation...")

	opt = Opt(:GN_CRS2_LM, length(ζ_0)+2) # number of parameters
	lower_bounds!(opt, 0)
	#lower_bounds!(opt, 0.1 .* [ζ_0..., δ_0]) # band around initial guess
	upper_bounds!(opt, [4.0, 3.0, 5.0, 100, 0.1]) # band around initial guess
	ftol_rel!(opt, 1e-3) # tolerance
	min_objective!(opt, loss_opt) # specify objective function
	min_f, min_x, ret = optimize(opt, [ζ_0...,θ_0,δ_0]) # run!

	# estimates
	ζ = min_x[1:end-2]
	θ = min_x[end-1]
	δ = min_x[end]

	println("Done! Bivariate logistic fit of λ with raw α. Estimates:")
	println("ζ: $ζ") 
	println("θ: $θ") 
	println("δ: $δ")
else # just use the initial guess
	ζ = ζ_0
	θ = θ_0
	δ = δ_0
end # estimate_rates


### Non-parametric objects ###

if compute_np_obj

	# precompute discount factors
	d = compute_d(δ, ψm_ψf)
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
		ξ = build_ξ(ζ) # construct ξ
		λ = inflate_λ(ξ ./ sqrt(sum(sng_conv["$msa"])), θ)

		# match probability (α)
		alpha["$msa"] = compute_alpha(λ, δ, ψm_ψf, marriages["$msa"], sng_conv["$msa"][2:end,:,:,2:end,:,:])

		# non-truncated α
		alpha_raw["$msa"] = compute_raw_alpha(λ, δ, ψm_ψf, marriages["$msa"], sng_conv["$msa"][2:end,:,:,2:end,:,:])

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
		write(file, "zeta", ζ)
		write(file, "theta", θ)
		write(file, "delta", δ)
		write(file, "alpha", alpha)
		write(file, "alpha_raw", alpha_raw)
		write(file, "surplus", surplus)
		write(file, "men_val", men_val)
		write(file, "wom_val", wom_val)
		write(file, "production", production)
		write(file, "avg_production", avg_production)
	end # do
end #if
