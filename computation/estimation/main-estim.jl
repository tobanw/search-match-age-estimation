using JLD

# TODO:
#	* SMM: estimate arrival rates --> need to load up data moments: MF(x),DF(x)

### PARAMS ### 

# arrival rates from `results/rate-param.csv`
#λ = 6.2e-9
#δ = 0.0117

# TEST: sensible arrival rates
λ = 0.000003
δ = 0.025

const r = 0.04 # discount rate
const ρ = 1.0 # aging rate
const β = 0.5 # wife share of surplus

# NOTE: must match `max.age` from `smooth-pop.r` and mortality data
const max_age = 65 # terminal age (inclusive)
const min_age = 25 # initial age (excluded)
const n_ages = max_age - min_age # excluding 25
const n_years = 7 # 2008-2014

# NOTE: must match `top.msa` from `smooth-pop.r`
const top_msa = (35620, 31080, 16980, 19100, 37980, 26420, 47900, 33100, 12060, 14460)
const n_msa = length(top_msa)


### DATA ###

# uncomment to overwrite `populations.jld` (to reload the smoothed population csv files)
#include("prepare-pops.jl")

# load saved data as dict
popu = load("results/populations.jld")

# load and trim age 25 from psi
ψ_m = popu["men_psi"][2:end,:,:]
ψ_f = popu["wom_psi"][2:end,:,:]

# load population dicts
men_sng = popu["men_sng"]
wom_sng = popu["wom_sng"]
men_tot = popu["men_tot"]
wom_tot = popu["wom_tot"]
marriages = popu["marriages"]

# estimation functions
include("estimate-values.jl")


### ESTIMATION ###

### Arrival Rates ###

#TODO: SMM routine for λ,δ

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
	men_val["$msa"], wom_val["$msa"] = value_function(λ, dc1μ, men_sng["$msa"], wom_sng["$msa"])

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
    write(file, "alpha", alpha)
    write(file, "surplus", surplus)
    write(file, "men_val", men_val)
    write(file, "wom_val", wom_val)
    write(file, "production", production)
    write(file, "avg_production", avg_production)
end # do
