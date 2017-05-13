# Use smoothed population masses and estimated arrival rates to recover marital production function by MSA
# Treat initial age as birth inflow (mass * phi), use u[2:end], n[2:end,2:end] for alpha, s, f


using Distributions

### Other Estimation Objects ###

"Compute array of d factors (discount factors on surplus)."
function compute_d(δ::Real, ψm_ψf::Array)
	d = 1 ./ (r + ρ + δ .+ ψm_ψf)
	# overwrite with terminal case: d^{T,T} has no ρ (aging stops)
	d[end,:,:,end,:,:] = 1 ./ (r + δ .+ ψm_ψf[end,:,:,end,:,:])
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
function compute_raw_alpha(λ::Array, δ::Real, ψm_ψf::Array, mar_init::Array, um_uf::Array)
	# trim off age 25
	m = mar_init[2:end,:,:,2:end,:,:] # mar_init[i] == mar[i-1] along the age dims

	""" PARALLEL VERSION
	α = convert(SharedArray, similar(m)) # instantiate alpha array, in shared memory
	# @parallel can't handle CartesianRange (5/5/2017), so need to manually index
	@sync @parallel for b in 1:n_ages # split task on wife age (faster because of column major)
		for a in 1:n_ages 
			for x1 in 1:2, x2 in 1:2, y1 in 1:2, y2 in 1:2
				# initial age as marriage inflows
				α[a,x1,x2,b,y1,y2] = ((m[a,x1,x2,b,y1,y2] * (ρ + δ + ψ_m[a,x1,x2] + ψ_f[b,y1,y2] )
									   - ρ * mar_init[a,x1,x2,b,y1,y2])
				                      / (λ[a,x1,x2,b,y1,y2] * (um_uf[a,x1,x2,b,y1,y2])
										 + δ * m[a,x1,x2,b,y1,y2]))
			end
		end
	end # parallel loop
	return α # raw array might not lie within [0,1]
	"""
	# fast vectorized version
	α = (m .* (ρ + δ .+ ψm_ψf) .- ρ .* mar_init[1:end-1,:,:,1:end-1,:,:]) ./ (λ .* um_uf .+ δ .* m)
	return α # raw array may not lie within [0,1]
end

"Compute α for a given MSA."
function compute_alpha(λ::Array, δ::Real, ψm_ψf::Array, mar_init::Array, um_uf::Array)
	α = compute_raw_alpha(λ, δ, ψm_ψf, mar_init, um_uf)
	# TODO: alternatively, could use smooth truncator
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

	for xy in CartesianRange(size(f))
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

#using Interpolations

logistic(x::Real; c=1.5, m=2.5, s=8.0) = c*exp(-(x-m)/s)/(1+exp(-(x-m)/s))^2

"Compute age-specific arrival rate vector ξ from parameter vector ζ"
function build_ξ(ζ::Vector)
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
	""" linear interpolation on knots
	grid = (knots,) # age-nodes for interpolation
	itp = interpolate(grid, ζ, Gridded(Linear()))
	return itp[collect(1-n_ages:n_ages-1)]
	"""
	# logistic in age-gap with free mean
    return logistic.(1-n_ages:n_ages-1; c=ζ[1], m=ζ[2], s=ζ[3])
end

"Reconstitute full λ array from vector of age-gap arrival rates."
function inflate_λ(lv::Vector, decay::Real)
	lam = Array(Float64, (n_ages,2,2,n_ages,2,2)) # NOTE: n_ages excludes age 25
	for b in 1:n_ages, a in 1:n_ages
		gapidx = (a-b) + n_ages # age gap aligned to lv index
		lam[a,:,:,b,:,:] = lv[gapidx] * logistic((a+b)/2, c=1, m=1, s=decay) # falls in age gap AND combined age of couple
	end
	return lam
end


#Functions for parallel computation of (for given MSA):
#	* male marriage flows
#	* female marriage flows
#	* male divorce flows
#	* female divorce flows

function compute_MF(λ::Array, um_uf::Array, α::Array)
	return λ .* um_uf .* α
end

#@everywhere 
function compute_DF_m(δ::Real, α::Array, m::Array)
	DF_m = zeros(m[:,:,:,1,1,1])
	for x in CartesianRange(size(DF_m)) # men
		# DF(x) = δ*∫(1-α(x,y))*m(x,y)dy
		DF_m[x] = δ * sum((1 .- α[x,:,:,:]) .* m[x,:,:,:])
	end
	return DF_m
end

#@everywhere 
function compute_DF_f(δ::Real, α::Array, m::Array)
	DF_f = zeros(m[1,1,1,:,:,:])
	for y in CartesianRange(size(DF_f)) # men
		# DF(y) = δ*∫(1-α(x,y))*m(x,y)dx
		DF_f[y] = δ * sum((1 .- α[:,:,:,y]) .* m[:,:,:,y])
	end
	return DF_f
end

"Compute model moments for a given MSA."
function model_moments(λ::Array, δ::Real, ψm_ψf::Array, mar_init::Array, um_uf::Array)
	# use raw alpha to help optimizer avoid solutions with alpha = 1
	α = compute_raw_alpha(λ, δ, ψm_ψf, mar_init, um_uf)
	
	m = mar_init[2:end,:,:,2:end,:,:]

	# non-parallel calls
	DF_m = compute_DF_m(δ, α, m)
	DF_f = compute_DF_f(δ, α, m)

	# spawn jobs on worker processes
	#job_DF_m = @spawn compute_DF_m(δ, α, m)
	#job_DF_f = @spawn compute_DF_f(δ, α, m)

	MF = compute_MF(λ, um_uf, α) # run locally while waiting (simple array product)

	# fetch results from workers
	#DF_m = fetch(job_DF_m)
	#DF_f = fetch(job_DF_f)

	# NOTE: max_age will be dropped in the estimation because of truncation
	return MF, DF_m, DF_f
end

"Minimum distance loss function per MSA."
function loss_msa(MF::Array, DF_m::Array, DF_f::Array,
				  dMF::Array, dDF_m::Array, dDF_f::Array,
				  wgt_men::Array, wgt_wom::Array, wgt_mar::Array)
	# assumes input arrays are from ages 26-64: already dropped min and max
	# weighted sum of squared percentage errors: half of integrated weights for DF
	# parallel on 2 workers
	#job_Dm = @spawn sum(0.5 * sum(wgt_wom) .* wgt_men .* ((DF_m .- dDF_m) ./ dDF_m).^2)
	#job_Df = @spawn sum(0.5 * sum(wgt_men) .* wgt_wom .* ((DF_f .- dDF_f) ./ dDF_f).^2)
	loss_Dm = sum(0.5 * sum(wgt_wom) .* wgt_men .* ((DF_m .- dDF_m) ./ dDF_m).^2)
	loss_Df = sum(0.5 * sum(wgt_men) .* wgt_wom .* ((DF_f .- dDF_f) ./ dDF_f).^2)

	loss_M = sum(wgt_mar .* ((MF .- dMF) ./ dMF).^2) # local while waiting

	return loss_M + loss_Dm + loss_Df #fetch(job_Dm) + fetch(job_Df)
end

"Objective function to minimize: distance between model and data moments."
function loss(ζ::Vector, θ::Real, δ::Real, ψm_ψf::Array,
			  mar_all::Dict{AbstractString,Array}, um_uf_all::Dict{AbstractString,Array},
	          dMF::Dict{AbstractString,Array},
			  dDF_m::Dict{AbstractString,Array}, dDF_f::Dict{AbstractString,Array},
	          wgt_men::Dict{AbstractString,Array}, wgt_wom::Dict{AbstractString,Array},
			  wgt_mar_all::Dict{AbstractString,Array})

	sse = 0.0 # initialize
	for msa in top_msa

		ξ = build_ξ(ζ) # construct ξ

		# reconstitute full λ array from age-gap ξ vector and θ
		# (\sum_x u_m(x))*(\sum_y u_f(y)) = U_m * U_f = \sum_x \sum_y u_m(x)*u_f(y)
		λ_gap = ξ ./ sqrt(sum(um_uf_all["$msa"]))
		λ = inflate_λ(λ_gap, θ)

		# model_moments returns ages 26-65 (only need age 25 for marriage inflows)
		MF, DF_m, DF_f = model_moments(λ, δ, ψm_ψf, mar_all["$msa"], um_uf_all["$msa"][2:end,:,:,2:end,:,:])

		# feed trimmed (age 26-64) moments and weights to loss function
		sse += loss_msa(MF[1:end-1,:,:,1:end-1,:,:], DF_m[1:end-1,:,:], DF_f[1:end-1,:,:],
						dMF["$msa"][2:end-1,:,:,2:end-1,:,:],
						dDF_m["$msa"][2:end-1,:,:], dDF_f["$msa"][2:end-1,:,:],
		                wgt_men["$msa"][2:end-1,:,:], wgt_wom["$msa"][2:end-1,:,:],
						wgt_mar_all["$msa"][2:end-1,:,:,2:end-1,:,:])
	end

	#gc() # julia bug doesn't garbage collect enough for SharedArray, runs out of memory

	# Verbose progress indicator 
	println(@sprintf("sse: %.5e, δ: %.3g, θ: %.3g, ζ: ", sse, δ, θ), ζ)

	return sse
end

"Objective function to pass to NLopt: requires vectors for `x` and `grad`."
function loss_opt(x::Vector, grad::Vector)
	return loss(x[1:end-2], x[end-1], x[end], ψm_ψf,
				marriages, sng_conv,
	            MF, men_DF, wom_DF,
	            men_tot, wom_tot, pop_conv)
end
