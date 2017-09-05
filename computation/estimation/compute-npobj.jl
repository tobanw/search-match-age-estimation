# Use smoothed population masses and estimated arrival rates to recover marital production function by MSA
# Treat initial age as birth inflow (mass * phi), use u[2:end], n[2:end,2:end] for alpha, s, f


### Non-parametric Objects ###

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
Φ_inv(x::Real) = Distributions.quantile(STDNORMAL, x)

"μ function, using inverse Mills ratio: E[z|z>q] = σ ϕ(q/σ) / (1 - Φ(q/σ))."
function μ(a::Real)
	st = Φ_inv(1-a) # pre-compute -s/σ = Φ^{-1}(1-a)
	return Distributions.pdf(STDNORMAL, st) - a * st
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
