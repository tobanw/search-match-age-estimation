# set parameters
ζ_mc = [0.85]
δ_mc = 0.025

ξ_mc = ζ_mc[1] # uniform arrivals
h_mc(x::Vector, y::Vector) = 1 - (x[1]-y[1])^2 - 0.3*abs(x[2]-y[2]) - 0.2*abs(x[3]-y[3]) # production function

edutypes = no_school?[1]:[1,2]
θ_mc = Vector[collect(linspace(0,1,40)), edutypes, [1,2]] # types: age, college, minority

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

mm_eqm = [SearchClosed(ξ_mc, δ_mc, r, 1,
					   θ_mc, θ_mc,
					   men_tot["$msa"][2:end,:,:],
					   wom_tot["$msa"][2:end,:,:],
					   h_mc; CRS=true, verbose=true, step_size=0.3)
		  for msa in top_msa] 

# parallel version segfaults (nlsolve issue) even with no workers
#mm_jobs = [@spawn SearchClosed(ξ_mc, δ_mc, r, 1,...
#mm_eqm = [fetch(job) for job in mm_jobs]

for (msa, mm) in zip(top_msa, mm_eqm) # simulate a single marriage market (NYC)
	# overwrite equilibrium objects and flow moments
	# leave age 25, which gets trimmed off in the estimation
	men_sng["$msa"][2:end,:,:] = mm.u_m
	wom_sng["$msa"][2:end,:,:] = mm.u_f
	sng_conv["$msa"] = outer_op(*, men_sng["$msa"], wom_sng["$msa"])

	MF["$msa"][2:end,:,:,2:end,:,:] = compute_MF(mm.λ * ones(mm.α),
									   sng_conv["$msa"][2:end,:,:,2:end,:,:],
									   mm.α)
	marriages["$msa"][2:end,:,:,2:end,:,:] = MF["$msa"][2:end,:,:,2:end,:,:] ./ (ψm_ψf + mm.δ * (1 - mm.α)) # back out eqm marriage measures

	men_DF["$msa"][2:end,:,:] = compute_DF_m(mm.δ, mm.α, marriages["$msa"][2:end,:,:,2:end,:,:])
	wom_DF["$msa"][2:end,:,:] = compute_DF_f(mm.δ, mm.α, marriages["$msa"][2:end,:,:,2:end,:,:])
end

# FIXME: estimation fails to recover... simulation has 1e5 times higher flows and lower singles
# verify that flows are correct, hard to imagine such high MF with fewer singles

# NOTE: this only overwrites the MSAs in top_msa, so might leave some originals
jldopen("results/sim-pops-static.jld", "w") do file  # open file for saving julia data
	# dictionaries: stocks per MSA
    write(file, "men_sng", men_sng)
    write(file, "wom_sng", wom_sng)
    write(file, "marriages", marriages)
    write(file, "sng_conv", sng_conv)

	# dictionaries: flows per MSA
    write(file, "MF", MF)
    write(file, "men_DF", men_DF)
    write(file, "wom_DF", wom_DF)
end

# load on all workers
@everywhere begin # load on all process in case running parallel
	# load saved data as dict
	simpp = JLD.load("results/sim-pops-static.jld")

	men_sng = simpp["men_sng"]
	wom_sng = simpp["wom_sng"]
	marriages = simpp["marriages"]
	sng_conv = simpp["sng_conv"]

	# load data	moments: marriage and divorce flow dicts (ages 25-65 inclusive)
	MF = simpp["MF"]
	men_DF = simpp["men_DF"]
	wom_DF = simpp["wom_DF"]
end
