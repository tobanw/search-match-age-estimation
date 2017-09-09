### USAGE OPTIONS ###

# select model
static_age = false # whether to use model with static age or with dynamic aging process
monte_carlo = false # use population supplies to compute equilibria with MarriageMarkets; static only for now

# select functional form for xi; if none selected, then constant ξ = ζ[1] is used
interpolated_xi = false
logistic_xi = false

reload_smooth = false # set true to reload the smoothed population csv files, otherwise use saved JLD
grid_search = false
estimate_rates = false # set true to run GMM estimation, otherwise just use ζ_0, δ_0 below
compute_np_obj = false

# select optimizer for arrival rate estimation
use_nlopt = false
use_bbopt = true


### SETUP ### 

# PARALLEL: multiprocess for grid_search
if grid_search || (estimate_rates && use_bbopt)
	addprocs() # add a worker process per core
	println("  > Using $(nprocs()-1) worker processes")
end

using JLD, Distributions # JLD: for saving / loading julia objects

if interpolated_xi
	using Interpolations
end


### PARAMS ### 

@everywhere begin # load on all process in case running parallel

	const r = 0.04 # discount rate
	const ρ = 1.0 # aging rate
	const β = 0.5 # wife share of surplus
	const STDNORMAL = Distributions.Normal()

	# NOTE: must match `max.age` from `smooth-pop.r` and mortality data
	const max_age = 65 # terminal age (inclusive)
	const min_age = 25 # initial age (excluded)
	const n_ages = max_age - min_age # excluding 25
	const n_years = 7 # 2008-2014

	# NOTE: must match `top.msa` from `smooth-pop.r` and `smooth-flows.r`
	const top_msa = (35620, 31080, 16980, 19100, 37980, 26420, 47900, 33100, 12060, 14460,
					 41860, 19820, 38060, 40140, 42660, 33460, 41740, 45300, 41180, 12580)
	const n_msa = length(top_msa)

end # begin block


### DATA ###

if reload_smooth # load from csv into dicts (and also save as JLD)
	println("  > Reloading populations...")
	@everywhere include("prepare-pops.jl")
else # load from existing JLD file
	@everywhere begin # load on all process in case running parallel
		# load saved data as dict
		pp = JLD.load("results/populations.jld")

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
	end # begin block
end # load data

### MODEL ###

# load estimation functions
@everywhere include("estim-functions.jl")

if monte_carlo # monte carlo simulation to verify estimation
	println("  > Simulating equilibria for Monte Carlo estimation...")
	using MarriageMarkets
	# overwrite equilibrium singles measures: men_sng, wom_sng, sng_conv and MF/DF
	#@everywhere include("monte-carlo-static.jl")
end

# static or dynamic aging model
if static_age
	println("  > Using static age model")
	@everywhere include("static-npobj.jl")
else
	println("  > Using dynamic age model")
	@everywhere include("compute-npobj.jl")
end

# functional form for arrival rates ξ
if interpolated_xi
	@everywhere build_ξ = interp_ξ
elseif logistic_xi
	@everywhere build_ξ = logistic_ξ
else
	@everywhere build_ξ = constant_ξ
end


### INITIAL GUESSES ###

"""
*** Static model estimates ***
TODO: Monte Carlo to make sure it can recover params from simulated eqm

Uniform λ:
ξ: 0.386
δ: 0.035
Loss: 859.02

BBopt:
ζ: [2.8223, 3.80494, 4.18702]
δ: 0.037614593278515844
Loss: 845.2647480122813

With avg age params:
ζ: [11.1629, 3.5432, 4.21858, 33.0407, 82.3968]
δ: 0.03740435275515733
Loss: 845.3516205677049

Interpolations: 1+3+6 zeta
ζ: [1.2, 1.733, 0.35, 0.65, 0.7 (max), 0.267, 0.3 (min), 0.633, 0.6, 0.6]
δ: 0.0375
Loss: 823.26

*** Dynamic model estimates ***

Uniform ξ:
ξ: 0.874752
δ: 0.02442846504409377
Loss: 1109.4199650520745

Non-sqrt moment weighting; Uniform ξ:
ξ: [0.79647]
δ: 0.02507793089870493
Loss: 1301.4755109302262
"""

ζ_0 = [0.87475] #
δ_0 = 0.024428 # arrival rate of love shocks


### ESTIMATION ###

### Arrival Rates: GMM estimation ###

if grid_search

	println("  > Running parameter grid search...")

	# with 8 processes, runs 54 evals per second: 3240/min, 195k/hour
	# with double interpolation, ran slower: 140k/hour (probably because interpolators are constructed inside the function)

	"""
	# logistic: ζ = [c,m,s1,s2] (age-gap) + [m,s] (avg age)
	ζ1grid = linspace(7., 16., 8) # (global multiplicative scale-factor for λ)
	ζ2grid = linspace(0., 4., 4) # (mean for age-gap meeting slowdown)
	ζ3grid = linspace(0.2, 0.4, 5) # (1/spread for negative age-gap meeting slowdown)
	ζ4grid = linspace(0.1, 0.25, 4) # (1/spread for positive age-gap meeting slowdown)
	ζ5grid = [1.,4.,10.,70.,100.] #linspace(10., 70., 5) # (mean for avg age search slowdown)
	ζ6grid = linspace(0.005, 0.03, 5) # (spread for avg age search slowdown)
	δgrid = [0.04] #linspace(0.03, 0.05, 3) # (divorce rate)
	"""

	# interpolations: ζ2-ζ4 avg age knots, ζ5-ζ7 age gap knots
	ζ1grid = linspace(0.1, 3.0, 300) # 8 (level)
	ζ2grid = linspace(1.2, 2.0, 4) # 4 age 3 knot (relative to level)
	ζ3grid = linspace(0.25, 0.55, 4) # 3 age 10 knot
	ζ4grid = linspace(0.55, 0.85, 4) # 3 age 40 knot
	#gap_knots = ([-39,-8,-3,2,7,12,39],)
	ζ5grid = linspace(0.1, 0.7, 4) # 3 (age gap: -39)
	ζ6grid = linspace(0.1, 0.6, 4) # 4 (age gap -8)
	ζ7grid = linspace(0.3, 0.8, 4) # 14 (age gap -3)
	ζ8grid = linspace(0.5, 0.9, 4) # 14 (age gap 7)
	ζ9grid = linspace(0.4, 0.7, 4) # 14 (age gap 12)
	ζ10grid = linspace(0.2, 0.8, 4) # 14 (age gap 39)
	δgrid = linspace(0.02, 0.05, 300) #[0.0375] # 1

	# list of jobs: for each ζ1
	gs_jobs = [(@spawn obj_landscaper(ζ1, #ζ2grid, ζ3grid, ζ4grid,
									  #ζ5grid, ζ6grid, ζ7grid, ζ8grid, ζ9grid, ζ10grid,
									  δgrid,
									  ψm_ψf, marriages, sng_conv, MF, men_DF, wom_DF,
									  men_tot, wom_tot, pop_conv)) for ζ1 in ζ1grid]

	result_list = [fetch(job) for job in gs_jobs] # list of strings
	result_str = join(result_list) # merge into single string

	open("results/loss-function/loss-grid.csv", "w") do f
		#write(f, "ζ1,ζ2,ζ3,ζ4,ζ5,ζ6,ζ7,ζ8,ζ9,ζ10,δ,LOSS\n")
		write(f, "ζ1,δ,LOSS\n")
		write(f, result_str)
	end # write to file
end # grid_search


### Optimization: GMM estimator ###

if estimate_rates # run optimizer for MD estimation
	println("  > Starting optimization routine for GMM estimation:")

	if use_nlopt
		println("  > Using NLopt...")
		using NLopt
		opt = Opt(:GN_CRS2_LM, length(ζ_0)+1) # number of parameters
		#opt = Opt(:GN_ESCH, length(ζ_0)+1) # doesn't converge...
		population!(opt, 256) # default is 10*(n+1)
		#lower_bounds!(opt, 0.1 * [ζ_0..., δ_0]) # band around initial guess
		lower_bounds!(opt, [0.01, 0.01])
		upper_bounds!(opt, [10., 0.1])
		#ftol_rel!(opt, 1e-12) # tolerance |Δf|/|f|
		ftol_abs!(opt, 1e-6) # tolerance |Δf|
		min_objective!(opt, loss_nlopt) # specify objective function
		min_f, min_x, ret = optimize(opt, [ζ_0...,δ_0]) # run!

	elseif use_bbopt
		println("  > Using BlackBoxOptim...")
		using BlackBoxOptim
		# NES algos can be run in parallel
		resbb = bboptimize(loss_bbopt; Method=:xnes, PopulationSize=64, Workers = workers(),
						   SearchRange = [(0.2, 1.6),# (1.5, 5.5), (0, 0.8), (0, 0.8),
										  #(6.,39.), (0, 0.1),
										  (0.01, 0.04)],
						   MinDeltaFitnessTolerance = 1e-3, MaxSteps = 200000,
						   TraceInterval = 30.0)

		min_x = best_candidate(resbb)
		min_f = best_fitness(resbb)
	
	else
		println("  > WARNING: No optimizer selected!")
	end

	# estimates
	ζ = min_x[1:end-1]
	δ = min_x[end]

	println("  > Done! Uniform λ with raw α.")
	println("ζ: $ζ") 
	println("δ: $δ")
	println("Loss: $min_f")
else # just use the initial guess
	ζ = ζ_0
	δ = δ_0
end # estimate_rates


### Non-parametric objects ###

if compute_np_obj

	println("  > Computing and saving non-parametric object estimates")

	# precompute discount factors
	d = compute_d(δ, ψm_ψf) # same function name for static and dynamic age
	if ~static_age
		c = compute_c(d)
		c1 = compute_c1(c)
	end

	avg_production = zeros(d) # initialize

	# save arrays (by MSA) in separate dicts, to be stored
	alpha = Dict{AbstractString, Array}()
	alpha_raw = Dict{AbstractString, Array}()
	surplus = Dict{AbstractString, Array}()
	men_val = Dict{AbstractString, Array}()
	wom_val = Dict{AbstractString, Array}()
	production = Dict{AbstractString, Array}()

	ξ = build_ξ(ζ) # construct ξ

	for msa in top_msa
		# U_m, U_f per marriage market
		λ = ξ / sqrt(sum(men_sng["$msa"]) * sum(wom_sng["$msa"]))

		# match probability (α)
		alpha["$msa"] = compute_alpha(λ, δ, ψm_ψf, marriages["$msa"], sng_conv["$msa"][2:end,:,:,2:end,:,:])

		# non-truncated α
		alpha_raw["$msa"] = compute_raw_alpha(λ, δ, ψm_ψf, marriages["$msa"], sng_conv["$msa"][2:end,:,:,2:end,:,:])

		if static_age
			# marital surplus (S)
			surplus["$msa"] = invert_alpha(alpha["$msa"])

			# precompute discounted μ array: solution of ∫ max{S(x,y,z),0}dG(z)
			dμ = compute_dμ(d, alpha["$msa"])

			# average value functions
			men_val["$msa"], wom_val["$msa"] = compute_value_functions(λ, dμ, men_sng["$msa"], wom_sng["$msa"])

			# marital production (f)
			production["$msa"] = compute_production(δ, ψ_m, ψ_f, dμ, surplus["$msa"],
													men_val["$msa"], wom_val["$msa"])
		else # dynamic aging model
			# marital surplus (S)
			surplus["$msa"] = invert_alpha(c1, alpha["$msa"])

			# precompute discounted μ array: solution of ∫ max{S(x,y,z),0}dG(z)
			dc1μ = compute_dc1μ(d, c1, alpha["$msa"])

			# average value functions
			men_val["$msa"], wom_val["$msa"] = compute_value_functions(λ, dc1μ, men_sng["$msa"], wom_sng["$msa"])

			# marital production (f)
			production["$msa"] = compute_production(δ, ψ_m, ψ_f, d, dc1μ, surplus["$msa"],
													men_val["$msa"], wom_val["$msa"])
		end

		# average across MSAs
		avg_production += production["$msa"]
	end # for MSA

	# final production estimate, average of MSA estimates
	avg_production = avg_production / n_msa

	# store in JLD format
	jldopen("results/estimates.jld", "w") do file  # open file for saving julia data
		write(file, "zeta", ζ)
		write(file, "delta", δ)
		# all arrays ages 26-65
		write(file, "alpha", alpha)
		write(file, "alpha_raw", alpha_raw)
		write(file, "surplus", surplus)
		write(file, "men_val", men_val)
		write(file, "wom_val", wom_val)
		write(file, "production", production)
		write(file, "avg_production", avg_production)
	end # do
end #if
