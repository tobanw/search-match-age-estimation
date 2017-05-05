# Use smoothed population masses and estimated arrival rates to recover marital production function by MSA
# Treat initial age as birth inflow (mass * phi), use u[2:end], n[2:end,2:end] for alpha, s, f


using Distributions

### Other Estimation Objects ###

"Compute array of d factors (discount factors on surplus)."
function compute_d(δ::Real, ψ_m::Array, ψ_f::Array)
	d = Array(Float64, (n_ages,2,2,n_ages,2,2)) # NOTE: n_ages excludes age 25
	for xy in CartesianRange(size(d))
		x = xy.I[1:3]
		y = xy.I[4:6]
		d[xy] = 1 ./ (r + ρ + δ + ψ_m[x...] + ψ_f[y...])
	end
	# terminal case: d^{T,T} has no ρ (aging stops)
	d_term = Array(Float64, (2,2,2,2))
	for xy in CartesianRange(size(d_term))
		x = xy.I[1:2]
		y = xy.I[3:4]
		d_term[xy] = 1 ./ (r + δ + ψ_m[end,x...] + ψ_f[end,y...]) #	ages (T,T)
	end
	# overwrite with terminal case
	d[end,:,:,end,:,:] = d_term

	return d
end

"Compute array of c factors (accumulated discount factors on future surpluses)."
function compute_c(d::Array)
	# 3 steps: (T,T) base case, then (a,T) and (T,b) boundaries, then interior
	c = 1 + ρ * d # initialize array: only T,T terminal value is correct
	# recursively fill boundary 
	for k in 1:n_ages-1
		c[end,:,:,end-k,:,:] = 1 + ρ * d[end,:,:,end-k,:,:] .* c[end,:,:,end-k+1,:,:] # age T husbands
		c[end-k,:,:,end,:,:] = 1 + ρ * d[end-k,:,:,end,:,:] .* c[end-k+1,:,:,end,:,:] # age T wives
	end
	# recursively fill layers from boundary inwards (via shrinking squares)
	for k in 1:n_ages-1
		c[1:end-k,:,:,1:end-k,:,:] = 1 + ρ * d[1:end-k,:,:,1:end-k,:,:] .* c[2:end-k+1,:,:,2:end-k+1,:,:]
	end

	return c
end

"Compute shifted array of c factors, c^{a+1,b+1}."
function compute_c1(c::Array)
	# 3 steps: (T,T) base case, then (a,T) and (T,b) boundaries, then interior
	c1 = ones(c) # initialize array: only T,T terminal value is correct, c^{T+1,T+1} = 1
	# recursively fill boundary 
	for k in 1:n_ages-1
		c1[end,:,:,end-k,:,:] = c[end,:,:,end-k+1,:,:] # age T husbands
		c1[end-k,:,:,end,:,:] = c[end-k+1,:,:,end,:,:] # age T wives
	end
	# fill interior
	c1[1:end-1,:,:,1:end-1,:,:] = c[2:end,:,:,2:end,:,:]

	return c1
end

# inverse cdf of standard normal distribution (for inverting α)
const STDNORMAL = Normal()
Φ_inv(x::Real) = quantile(STDNORMAL, x)

"μ function, using inverse Mills ratio: E[z|z>q] = σ ϕ(q/σ) / (1 - Φ(q/σ))."
function μ(a::Real)
	st = Φ_inv(1-a) # pre-compute -s/σ = Φ^{-1}(1-a)
	return pdf(STDNORMAL, st) - a * st
end

"Compute d*c1*μ(α) array"
function compute_dc1μ(d::Array, c1::Array, α::Array)
	return d .* c1 .* μ.(α)
end


### Estimation Functions ###

"Compute α for a given MSA."
function compute_raw_alpha(λ::Array, δ::Real, ψ_m::Array, ψ_f::Array,
	                       mar_init::Array, um_init::Array, uf_init::Array)
	# trim off age 25
	m = mar_init[2:end,:,:,2:end,:,:] # mar_init[i] == mar[i-1] along the age dims
	u_m = um_init[2:end,:,:]
	u_f = uf_init[2:end,:,:]

	α = similar(m) # instantiate alpha array
	for xy in CartesianRange(size(m))
		x = xy.I[1:3]
		y = xy.I[4:6]
		# initial age as marriage inflows
		α[xy] = ((m[xy] * (ρ + δ + ψ_m[x...] + ψ_f[y...] ) - ρ * mar_init[xy])
				 / (λ[xy] * (u_m[x...] * u_f[y...]) + δ * m[xy]))
	end
	return α # raw array may not lie within [0,1]
end

"Compute α for a given MSA."
function compute_alpha(λ::Array, δ::Real, ψ_m::Array, ψ_f::Array,
					   mar_init::Array, um_init::Array, uf_init::Array)

	α = compute_raw_alpha(λ, δ, ψ_m, ψ_f, mar_init, um_init, uf_init)

	return clamp.(α, 1e-8, 1 - 1e-8) # enforce 0 < α < 1
end

"Compute surplus s by inverting α."
function invert_alpha(c1::Array, α::Array)
	return -c1 .* Φ_inv.(1 .- α) # use shifted c
end

"Compute average value functions."
function compute_value_functions(λ::Array, dc1μ::Array, um_init::Array, uf_init::Array)
	# trim off age 25
	u_m = um_init[2:end,:,:]
	u_f = uf_init[2:end,:,:]

	# use average value functions v = (r + ρ + ψ_m + ψ_f)V 
	# initialize v: augmented with v^{T+1}=0 for convenience
	v_m = zeros(um_init)
	v_f = zeros(uf_init)

	# solve backwards because of continuation value
	for j in 0:n_ages-1
		k = n_ages - j # work back from terminal age
		# NOTE: CartesianIndex doesn't need to be splatted
		# female value function
		for y in CartesianRange(size(u_f[1,:,:])) # use trimmed size
			v_f[k,y] = ρ * v_f[k+1,y] + β * sum(λ[:,:,:,k,y] .* dc1μ[:,:,:,k,y] .* u_m)
		end
		# male value function
		for x in CartesianRange(size(u_m[1,:,:])) # use trimmed size
			v_m[k,x] = ρ * v_m[k+1,x] + (1-β) * sum(λ[k,x,:,:,:] .* dc1μ[k,x,:,:,:] .* u_f)
		end
	end

	return v_m[1:end-1,:,:], v_f[1:end-1,:,:] # trim off age T+1
end

"Compute production function."
function compute_production(δ::Real, ψ_m::Array, ψ_f::Array, d::Array, dc1μ::Array,
                            s::Array, v_m::Array, v_f::Array)
	f = similar(s) # instantiate production array

	# interior: all but (T,T)
	for xy in CartesianRange(size(f)) # loop over interior
		# unpack indices: [xy] == [a,x...,b,y...]
		a = xy.I[1]
		x = xy.I[2:3]
		b = xy.I[4]
		y = xy.I[5:6]

		# continuations from aging
		ρa = ρ
		ρA = ρ
		ρb = ρ
		ρB = ρ
		ρAB = ρ

		if a == n_ages-1 # adjust discount factor on continuation value
			ρa = 0
		elseif a == n_ages # shut off continuation value
			ρA = 0
		end
		if b == n_ages-1 # adjust discount factor on continuation value
			ρb = 0
		elseif b == n_ages # shut off continuation value
			ρB = 0
		end

		if a == b == n_ages # shut off continuation surplus
			ρAB = 0
		end

		# truncate indices at T to handle boundary (except for T,T case)
		f[xy] = (s[xy] + v_m[a,x...] + v_f[b,y...] - δ * dc1μ[xy]
				 - ρAB * d[min(a+1,end),x...,min(b+1,end),y...] * s[min(a+1,end),x...,min(b+1,end),y...] # ρAB shuts off s^{T+1,T+1} case
				 - ρA * v_m[min(a+1,end),x...] / (r + ρa + ψ_m[min(a+1,end),x...])
				 - ρB * v_f[min(b+1,end),y...] / (r + ρb + ψ_f[min(b+1,end),y...])) # ρB shuts off V^{T+1}
	end

	return f
end


### GMM Estimation of Arrival Rates ###

using Interpolations

# estimation function: inputs are data moments and pop masses
#	* batch call all msa
#	* model moments function: given params and masses, compute alpha and generate model moments
#	* criterion function: given model moments and data moments, compute weighted SSE
#	* pipe into optimizer, return results
#	* Standard Errors? curvature at optimum

"Compute age-specific arrival rate vector ξ from parameter vector ζ"
function build_ξ(ζ::Vector, knots::Vector)
	# exponential decay: ξ[k] = exp(-ζ[1] - ζ[2]*k - ζ[3]*k^2 - ...)
	# harmonic decay: ξ[k] = 1/(ζ[1] + ζ[2]*k + ζ[3]*k^2 + ...)
	# gaussian kernel: ξ[k] = ζ[1]*exp(-(k/ζ[2])^2)
	# logistic kernel: ζ[1]/ζ[2] * exp(-x/ζ[2])/(1+exp(-x/ζ[2]))^2
	""" parametric construction
	ξ = zeros(n_ages)
	for k in 0:n_ages-1
		ξ[k+1] = ζ[1]/ζ[2] * exp(-(k/ζ[2]))/(1+exp(-k/ζ[2]))^2
	end
	return ξ
	"""
	# linear interpolation on knots
	grid = (knots,) # age-nodes for interpolation
	itp = interpolate(grid, ζ, Gridded(Linear()))
	return itp[collect(0:n_ages-1)]
end

"Reconstitute full λ array from vector of age-gap arrival rates."
function inflate_λ(lv::Vector)
	lam = Array(Float64, (n_ages,2,2,n_ages,2,2)) # NOTE: n_ages excludes age 25
	for xy in CartesianRange(size(lam))
		# unpack indices: [xy] == [a,x...,b,y...]
		a = xy.I[1]
		b = xy.I[4]

		# simple age gap
		gapidx = abs(a-b) + 1 # age gap aligned to lv index

		lam[xy] = lv[gapidx]
	end
	return lam
end

#Functions for parallel computation of (for given MSA):
#	* male marriage flows
#	* female marriage flows
#	* male divorce flows
#	* female divorce flows

#@everywhere 
function compute_MF_m(λ::Array, um::Array, uf::Array, α::Array)
	MF_m = zeros(um)
	for x in CartesianRange(size(MF_m)) # men
		# MF(x) = um(x)*∫λ(x,y)*uf(y)*α(x,y)dy
		MF_m[x] = um[x] * sum(λ[x,:,:,:] .* uf .* α[x,:,:,:])
	end
	return MF_m
end

function compute_MF_f(λ::Array, um::Array, uf::Array, α::Array)
	MF_f = zeros(uf)
	for y in CartesianRange(size(MF_f)) # women
		# MF(y) = uf(y)*∫λ(x,y)*um(x)*α(x,y)dx
		MF_f[y] = uf[y] * sum(λ[:,:,:,y] .* um .* α[:,:,:,y])
	end
	return MF_f
end

function compute_DF_m(δ::Real, α::Array, m::Array)
	DF_m = zeros(m[:,:,:,1,1,1])
	for x in CartesianRange(size(DF_m)) # men
		# DF(x) = δ*∫(1-α(x,y))*m(x,y)dy
		DF_m[x] = δ * sum((1 .- α[x,:,:,:]) .* m[x,:,:,:])
	end
	return DF_m
end

function compute_DF_f(δ::Real, α::Array, m::Array)
	DF_f = zeros(m[1,1,1,:,:,:])
	for y in CartesianRange(size(DF_f)) # men
		# DF(y) = δ*∫(1-α(x,y))*m(x,y)dx
		DF_f[y] = δ * sum((1 .- α[:,:,:,y]) .* m[:,:,:,y])
	end
	return DF_f
end

"Compute model moments for a given MSA."
function model_moments(λ::Array, δ::Real, ψ_m::Array, ψ_f::Array,
					   mar_init::Array, um_init::Array, uf_init::Array)
	# use raw alpha to help optimizer avoid solutions with alpha = 1
	α = compute_raw_alpha(λ, δ, ψ_m, ψ_f, mar_init, um_init, uf_init)
	
	m = mar_init[2:end,:,:,2:end,:,:]
	u_m = um_init[2:end,:,:]
	u_f = uf_init[2:end,:,:]

	MF_m = compute_MF_m(λ, u_m, u_f, α)
	MF_f = compute_MF_f(λ, u_m, u_f, α)
	DF_m = compute_DF_m(δ, α, m)
	DF_f = compute_DF_f(δ, α, m)

	# spawn jobs on worker processes
	#job_MF_m = @spawn compute_MF_m(λ, u_m, u_f, α)
	#job_MF_f = @spawn compute_MF_f(λ, u_m, u_f, α)
	#job_DF_m = @spawn compute_DF_m(δ, α, m)
	#job_DF_f = @spawn compute_DF_f(δ, α, m)

	#println("Jobs spawned!")

	#MF_m = fetch(job_MF_m)
	#MF_f = fetch(job_MF_f)	
	#DF_m = fetch(job_DF_m)
	#DF_f = fetch(job_DF_f)

	#println("Jobs fetched!")
	# NOTE: max_age will be dropped in the estimation because of truncation
	return MF_m, DF_m, MF_f, DF_f
end

"Minimum distance loss function per MSA."
function loss_msa(MF_m::Array, DF_m::Array, MF_f::Array, DF_f::Array,
				  dMF_m::Array, dDF_m::Array, dMF_f::Array, dDF_f::Array,
				  wgt_men::Array, wgt_wom::Array)
	# assumes input arrays are from ages 26-64: already dropped min and max
	# weighted sum of squared errors
	wsse = ( sum(wgt_men .* (MF_m .- dMF_m).^2)
	       + sum(wgt_men .* (DF_m .- dDF_m).^2)
	       + sum(wgt_wom .* (MF_f .- dMF_f).^2)
	       + sum(wgt_wom .* (DF_f .- dDF_f).^2) )

	return wsse
end

"Objective function to minimize: distance between model and data moments."
function loss(ζ::Vector, δ::Real, knots::Vector,
	          ψ_m::Array, ψ_f::Array, mar_all::Dict{AbstractString,Array},
	          um_all::Dict{AbstractString,Array}, uf_all::Dict{AbstractString,Array},
	          dMF_m::Dict{AbstractString,Array}, dDF_m::Dict{AbstractString,Array},
	          dMF_f::Dict{AbstractString,Array}, dDF_f::Dict{AbstractString,Array},
	          wgt_men::Dict{AbstractString,Array}, wgt_wom::Dict{AbstractString,Array})

	sse = 0.0 # initialize
	for msa in top_msa

		ξ = build_ξ(ζ, knots) # construct ξ

		# reconstitute full λ array from age-gap ξ vector
		λ_gap = ξ ./ sqrt(sum(um_all["$msa"]) * sum(uf_all["$msa"])) # U_m, U_f per marriage market
		λ = inflate_λ(λ_gap)

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

	# Verbose progress indicator 
	println(@sprintf("sse: %.8e, δ: %.3g, ζ: ", sse, δ), ζ)
	return sse
end

"Objective function to pass to NLopt: requires vectors for `x` and `grad`."
function loss_opt(x::Vector, grad::Vector)
	return loss(x[1:end-1], x[end], znodes, ψ_m, ψ_f,
				marriages, men_sng, wom_sng,
	            men_MF, men_DF, wom_MF, wom_DF,
	            men_tot, wom_tot)
end

