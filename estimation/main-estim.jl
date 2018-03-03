### USAGE OPTIONS ###

# NOTE: some settings are in the @everywhere block below
#	* set `no_edu` and `no_race`
#	* set `min_age` and `max_age`

# select model
static_age = false # whether to use model with static age or with dynamic aging process
monte_carlo = false # use population supplies to compute equilibria with MarriageMarkets; static only for now

# select functional forms for xi and delta; if none selected, then constant ζ[1] is used
interpolated_xi = false
logistic_xi = false

interpolated_delta = false
logistic_delta = false

# what to run?
reload_smooth = true # set true to reload the smoothed population csv files, otherwise use saved JLD
grid_search = false # NOTE: need to manually set grids and modify `obj_landscaper` per model
estimate_rates = true # set true to run GMM estimation, otherwise just use ζx_0, ζd_0 below
compute_np_obj = true

data_dir = "data/racedu24/"
pop_file = "results/populations.jld"
save_path = "results/est-dyn-const-racedu-bw24.jld"
csv_dir = "results/estimates-csv/dynamic-const/racedu24/"

# select optimizer for arrival rate estimation
use_nlopt = true
use_bbopt = false


### SETUP ### 

# PARALLEL: multiprocess for grid_search or bbopt
if grid_search || (estimate_rates && use_bbopt)
	addprocs() # add a worker process per core
	info("Using $(nprocs()-1) worker processes")
end

using JLD, Distributions, Interpolations, DataFrames, Query


### PARAMS ### 

@everywhere begin # load on all process in case running parallel

	const no_edu = false # whether to exclude college as a static type
	const no_race = false # whether to exclude race as a static type

	const r = 0.04 # discount rate
	const ρ = 1.0 # aging rate
	const β = 0.5 # wife share of surplus
	const STDNORMAL = Distributions.Normal()

	# NOTE: must match `max.age` from `smooth-pop.r` and mortality data
	const max_age = 65 # terminal age (inclusive)
	const min_age = 25 # initial age (excluded: 18 or 25 with edu)
	const n_ages = max_age - min_age # excluding min_age (inflow group)

	# NOTE: must match `top.msa` from `smooth-pop.r` and `smooth-flows.r`
	const top_msa = (35620, 31080, 16980, 19100, 37980, 26420, 47900, 33100, 12060, 14460,
	                 41860, 19820, 38060, 40140, 42660, 33460, 41740, 45300, 41180, 12580)
	const n_msa = length(top_msa)

	if no_edu && no_race
		const dim_ind = (n_ages,1,1)
		const dim_mar = (n_ages,1,1,n_ages,1,1)
	elseif no_edu
		const dim_ind = (n_ages,1,2)
		const dim_mar = (n_ages,1,2,n_ages,1,2)
	elseif no_race
		const dim_ind = (n_ages,2,1)
		const dim_mar = (n_ages,2,1,n_ages,2,1)
	else
		const dim_ind = (n_ages,2,2)
		const dim_mar = (n_ages,2,2,n_ages,2,2)
	end
end # begin block


### DATA ###

if reload_smooth # load from csv into dicts (and also save as JLD)
	info("Reloading populations...")
	include("prepare-pops.jl")
end

@everywhere begin # load on all process in case running parallel
	# load saved data as dict
	pp = JLD.load(pop_file)


	# load and trim age 25 from psi, summing out edu
	ψ_m = pp["men_psi"][2:end,:,:] # ages 26-65
	ψ_f = pp["wom_psi"][2:end,:,:]
	ψm_ψf = pp["psi_outer"][2:end,:,:,2:end,:,:]

	# load population stock dicts (ages 25-65 inclusive)
	men_sng = pp["men_sng"]
	wom_sng = pp["wom_sng"]
	men_tot = pp["men_tot"]
	wom_tot = pp["wom_tot"]
	marriages = pp["marriages"]
	mar_outflow = pp["mar_outflow"]
	sng_outer = pp["sng_outer"]
	pop_outer = pp["pop_outer"]

	# load data	moments: marriage and divorce flow dicts (ages 25-65 inclusive)
	MF = pp["MF"]
	men_DF = pp["men_DF"]
	wom_DF = pp["wom_DF"]
end # begin block


### MODEL ###

# load estimation functions
@everywhere include("estim-functions.jl")

# static or dynamic aging model
if static_age
	info("Using static age model")
	@everywhere include("static-npobj.jl")
else
	info("Using dynamic age model")
	@everywhere include("compute-npobj.jl")
end

info("Ages: $(min_age+1) - $max_age; Educ: $(!no_edu); Race: $(!no_race)")

# functional form for arrival rates ξ
if interpolated_xi
	info("Using interpolated ξ")
	@everywhere build_ξ = interp_rate
elseif logistic_xi
	info("Using logistic ξ")
	@everywhere build_ξ = logistic_rate
else
	info("Using constant ξ")
	@everywhere build_ξ = constant_rate
end

# functional form for arrival rates δ
if interpolated_delta
	info("Using interpolated δ")
	@everywhere build_δ = interp_rate
elseif logistic_delta
	info("Using logistic δ")
	@everywhere build_δ = logistic_rate
else
	info("Using constant δ")
	@everywhere build_δ = constant_rate
end

if monte_carlo # monte carlo simulation to verify estimation
	info("Simulating equilibria for Monte Carlo estimation...")
	using MarriageMarkets
	# overwrite equilibrium singles measures: men_sng, wom_sng, sng_outer and MF/DF
	include("monte-carlo-static.jl")
end


### INITIAL GUESSES ###

"""
*** Static model estimates ***

Uniform λ:
ξ: 0.386
δ: 0.035
Loss: 859.02

Interpolations: 1+3+6 zeta
ζ: [1.2, 1.733, 0.35, 0.65, 0.7 (max), 0.267, 0.3 (min), 0.633, 0.6, 0.6]
δ: 0.0375
Loss: 823.26

*** Dynamic model estimates ***

Uniform ξ,δ (vs age-only 08-15):
ξ: 0.797566 (0.816589)
δ: 0.025083746441886164 (0.0247069)
Loss: 1301.47469321623 (45247.2472546416)

Logistic ξ,δ (symmetric 5 param, but only x1,x3,d1,d5 are free) (vs age-only):
ζx_0: [1.61508, 2, 0.283235, 0, 0] ([1.10436, 2, 0.000192278, 0, 0])
ζd_0: [0.047462, 2, 0, 0, 0.0682615] ([0.0338854, 2, 0, 0, 0.0582711])
Loss: 984.5798779676335 (28901.21320680169)

Interpolated ξ:
ζ: [6.70401, 0.288831, 0.0849464, 0.0212671, 0.001, 0.0162982, 0.775318, 0.667261, 0.0122603, 0.0126533]
δ: 0.030691015387483426
Loss: 1193.6704652231151
"""

ζx_0 = [3.0]
ζd_0 = [0.02]


### ESTIMATION ###

### Arrival Rates: GMM estimation ###

if grid_search

	info("Running parameter grid search...")

	# with 8 processes, runs 54 evals per second: 3240/min, 195k/hour
	# with double interpolation, ran slower: 140k/hour (probably because interpolators are constructed inside the function)

	# logistic: ζ = [c,m,s1,s2] (age-gap) + [m,s] (avg age)
	ζx1grid = linspace(0.01, 10, 1000) # (global multiplicative scale-factor for λ)
	#ζx2grid = [2] #linspace(0., 4., 4) # (mean for age-gap meeting slowdown)
	#ζx3grid = linspace(0.15, 0.35, 8) # (1/spread for age-gap meeting slowdown)
	#ζx4grid = [0] #linspace(10., 70., 5) # (mean for avg age search slowdown)
	#ζx5grid = [0] #linspace(0, 0.08, 3) # (1/spread for avg age search slowdown)

	ζd1grid = linspace(0.001, 0.03, 450) # (global multiplicative scale-factor for δ)
	#ζd2grid = [2] #linspace(0., 4., 4) # (mean for age-gap slowdown)
	#ζd3grid = [0] #linspace(0, 0.3, 4) # (1/spread for age-gap slowdown)
	#ζd4grid = [0] #linspace(-12, 0, 8) # (mean for avg age slowdown)
	#ζd5grid = linspace(0.05, 0.09, 10) # (1/spread for avg age slowdown)

	"""
	# interpolations: ζ2-ζ4 avg age knots, ζ5-ζ7 age gap knots
	ζ1grid = linspace(1, 24, 8) # 8 (level)
	ζ2grid = linspace(0.1, 0.8, 6) # 4 age 8 knot (relative to level)
	ζ3grid = linspace(0.01, 0.25, 3) # 3 age 20 knot
	ζ4grid = linspace(0.001, 0.25, 3) # 3 age 40 knot
	gap_knots = ([-39, -15, -2, 2, 6, 19, 39])
	ζ5grid = linspace(0.001, 0.05, 3) # 3 (age gap: -39)
	ζ6grid = linspace(0.01, 0.1, 3) # 3 (age gap -15)
	ζ7grid = linspace(0.5, 1.2, 4) # 3 (age gap -2)
	ζ8grid = linspace(0.5, 1.5, 5) # 3 (age gap 6)
	ζ9grid = linspace(0.001, 0.05, 3) # 3 (age gap 19)
	ζ10grid = linspace(0.005, 0.05, 3) # 3 (age gap 39)
	δgrid = linspace(0.02, 0.16, 6) #[0.0375] # 1
	"""

	# list of jobs: for each ζ1
	gs_jobs = [(@spawn obj_landscaper(ζx1,# ζx2grid, ζx3grid, ζx4grid, ζx5grid,
	                                  ζd1grid,# ζd2grid, ζd3grid, ζd4grid, ζd5grid,
	                                  ψm_ψf, marriages, sng_outer, mar_outflow, MF, men_DF, wom_DF,
	                                  men_tot, wom_tot, pop_outer)) for ζx1 in ζx1grid]

	result_list = [fetch(job) for job in gs_jobs] # list of strings
	result_str = join(result_list) # merge into single string

	open("results/loss-function/loss-grid.csv", "w") do f
		#write(f, "ζx1,ζx2,ζx3,ζx4,ζx5,ζd1,ζd2,ζd3,ζd4,ζd5,LOSS_MF,LOSS_DF,LOSS\n")
		write(f, "ξ,δ,LOSS_MF,LOSS_DF,LOSS\n")
		write(f, result_str)
	end # write to file
end # grid_search


### Optimization: GMM estimator ###

if estimate_rates # run optimizer for MD estimation
	info("Starting optimization routine for GMM estimation")

	if use_nlopt
		info("Using NLopt...")
		using NLopt
		opt = Opt(:LN_COBYLA, length(ζx_0)+length(ζd_0)) # number of parameters
		#opt = Opt(:GN_CRS2_LM, length(ζx_0)+length(ζd_0)) # number of parameters
		#population!(opt, 256) # default is 10*(n+1)
		#lower_bounds!(opt, 0.1 * [ζx_0..., ζd_0...]) # band around initial guess
		#upper_bounds!(opt, 10 * [ζx_0..., ζd_0...]) # band around initial guess
		lower_bounds!(opt, [0.01, 0])
		upper_bounds!(opt, [500, 1.])
		#lower_bounds!(opt, [0.1, -10, 0, -50, 0, 0.05, -10, 0, -50, 0])
		#upper_bounds!(opt, [300,  15, 1,  50, 1, 5,     15, 1,  50, 1.])
		#ftol_rel!(opt, 1e-12) # tolerance |Δf|/|f|
		ftol_abs!(opt, 1e-8) # tolerance |Δf|
		min_objective!(opt, loss_nlopt) # specify objective function
		#min_objective!(opt, (x,gr)->loss_nlopt([x[1],2,x[2],0,0,x[3],2,0,0,x[4]], gr)) # specify objective function
		min_f, min_x, ret = optimize(opt, [ζx_0..., ζd_0...]) # run!

	elseif use_bbopt
		info("Using BlackBoxOptim...")
		using BlackBoxOptim
		# NES algos can be run in parallel
		resbb = bboptimize(loss_bbopt; Method=:xnes,# PopulationSize=32, Workers = workers(),
		                   SearchRange = [(0.5, 1.2),# (1.5, 5.5), (0, 0.8), (0, 0.8),
		                                  (6.,39.), (0, 0.1),
		                                  (0.01, 0.04)],
		                   MinDeltaFitnessTolerance = 1e-3, MaxSteps = 200000,
		                   TraceInterval = 30.0)

		min_x = best_candidate(resbb)
		min_f = best_fitness(resbb)
	
	else
		warn("No optimizer selected!")
	end

	# estimates
	ζx = min_x[1:length(ζx_0)]
	ζd = min_x[end-length(ζd_0)+1:end]

	info("Done!")
	println("ζx: $ζx") 
	println("ζd: $ζd")
	println("Loss: $min_f")
else # just use the initial guess
	ζx = ζx_0
	ζd = ζd_0
end # estimate_rates


### Non-parametric object estimates ###

if compute_np_obj
	info("Computing and saving non-parametric object estimates")

	# save arrays (by MSA) in separate dicts, to be stored
	alpha = Dict{AbstractString, Array}()
	alpha_raw = Dict{AbstractString, Array}()
	surplus = Dict{AbstractString, Array}()
	men_val = Dict{AbstractString, Array}()
	wom_val = Dict{AbstractString, Array}()
	production = Dict{AbstractString, Array}()
	mMF = Dict{AbstractString, Array}()
	mDF_m = Dict{AbstractString, Array}()
	mDF_f = Dict{AbstractString, Array}()

	ξ = build_ξ(ζx) # construct ξ
	δ = build_δ(ζd) # construct ξ

	# precompute discount factors
	d = compute_d(δ, ψm_ψf) # same function name for static and dynamic age
	if ~static_age
		c = compute_c(d)
		c1 = compute_c1(c)
	end

	avg_production = zeros(d) # initialize

	for msa in top_msa
		# U_m, U_f per marriage market
		λ = ξ / sqrt(sum(men_sng["$msa"]) * sum(wom_sng["$msa"]))

		# match probability (α)
		alpha_raw["$msa"] = compute_raw_alpha(λ, δ, ψm_ψf, marriages["$msa"],
		                                      sng_outer["$msa"][2:end,:,:,2:end,:,:],
		                                      mar_outflow["$msa"][2:end,:,:,2:end,:,:])

		# truncated α
		alpha["$msa"] = clamp.(alpha_raw["$msa"], 1e-6, 1 - 1e-6) # enforce 0 < α < 1

		mMF["$msa"], mDF_m["$msa"], mDF_f["$msa"] = model_moments(λ, δ, ψm_ψf, marriages["$msa"],
		                                                          sng_outer["$msa"][2:end,:,:,2:end,:,:],
		                                                          mar_outflow["$msa"][2:end,:,:,2:end,:,:])

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
	jldopen(save_path, "w") do file  # open file for saving julia data
		write(file, "z-xi", ζx)
		write(file, "z-delta", ζd)
		# all arrays ages 26-65
		write(file, "alpha", alpha)
		write(file, "alpha_raw", alpha_raw)
		write(file, "surplus", surplus)
		write(file, "men_val", men_val)
		write(file, "wom_val", wom_val)
		write(file, "production", production)
		write(file, "avg_production", avg_production)
		write(file, "mMF", mMF)
		write(file, "mDF_m", mDF_m)
		write(file, "mDF_f", mDF_f)
	end # do

	### Write results arrays to CSV ###

	function msa_martab(d::Dict)
		df = DataFrame(AGE_M = Int64[], COLLEGE_M = Int64[], MINORITY_M = Int64[],
		               AGE_F = Int64[], COLLEGE_F = Int64[], MINORITY_F = Int64[],
		               MSA = Int64[], VALUE = Float64[])

		for msa in top_msa
			for idx in CartesianRange(size(d["$msa"]))
				i = idx.I
				push!(df, @data([min_age+i[1],i[2],i[3],min_age+i[4],i[5],i[6],msa,d["$msa"][idx]]))
			end
		end
		return df
	end

	function msa_indtab(d::Dict)
		df = DataFrame(AGE = Int64[], COLLEGE = Int64[], MINORITY = Int64[],
		               MSA = Int64[], VALUE = Float64[])

		for msa in top_msa
			for idx in CartesianRange(size(d["$msa"]))
				i = idx.I
				push!(df, @data([min_age+i[1],i[2],i[3],msa,d["$msa"][idx]]))
			end
		end
		return df
	end

	# convert to CSV for plotting
	writetable(joinpath(csv_dir, "prod.csv"), msa_martab(production))
	writetable(joinpath(csv_dir, "alpha.csv"), msa_martab(alpha))
	writetable(joinpath(csv_dir, "surplus.csv"), msa_martab(surplus))
	writetable(joinpath(csv_dir, "men_val.csv"), msa_indtab(men_val))
	writetable(joinpath(csv_dir, "wom_val.csv"), msa_indtab(wom_val))
	writetable(joinpath(csv_dir, "mMF.csv"), msa_martab(mMF))
	writetable(joinpath(csv_dir, "mDF_m.csv"), msa_indtab(mDF_m))
	writetable(joinpath(csv_dir, "mDF_f.csv"), msa_indtab(mDF_f))
end # compute_np_obj
