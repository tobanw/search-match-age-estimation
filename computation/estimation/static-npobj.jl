# Use smoothed population masses and estimated arrival rates to recover marital production function by MSA
# Treat initial age as birth inflow (mass * phi), use u[2:end], n[2:end,2:end] for alpha, s, f


### Other Estimation Objects ###

"Compute array of d factors (discount factors on surplus)."
function compute_d(δ::Real, ψm_ψf::Array)
	return 1 ./ (r + δ .+ ψm_ψf)
end

# inverse cdf of standard normal distribution (for inverting α)
Φ_inv(x::Real) = Distributions.quantile(STDNORMAL, x)

"μ function, using inverse Mills ratio: E[z|z>q] = σ ϕ(q/σ) / (1 - Φ(q/σ))."
function μ(a::Real)
	st = Φ_inv(1-a) # pre-compute -s/σ = Φ^{-1}(1-a)
	return Distributions.pdf(STDNORMAL, st) - a * st
end

"Compute d*μ(α) array"
function compute_dμ(d::Array, α::Array)
	return d .* μ.(α)
end


### Estimation Functions ###

"Compute α for a given MSA."
function compute_raw_alpha(λ::Array, δ::Real, ψm_ψf::Array, mar_init::Array, um_uf::Array)
	# trim off age 25
	m = mar_init[2:end,:,:,2:end,:,:] # mar_init[i] == mar[i-1] along the age dims
	# fast vectorized version
	α = (m .* (δ + ψm_ψf)) ./ (λ .* um_uf + δ * m)
	return α # raw array may not lie within [0,1]
end

"Compute α for a given MSA."
function compute_alpha(λ::Array, δ::Real, ψm_ψf::Array, mar_init::Array, um_uf::Array)
	α = compute_raw_alpha(λ, δ, ψm_ψf, mar_init, um_uf)
	# TODO: alternatively, could use smooth truncator
	return clamp.(α, 1e-8, 1 - 1e-8) # enforce 0 < α < 1
end

"Compute surplus s by inverting α."
function invert_alpha(a::Array)
	return -Φ_inv.(1 - a)
end

"Compute average value functions."
function compute_value_functions(λ::Array, dμ::Array, um_init::Array, uf_init::Array)
	# trim off age 25
	u_m = um_init[2:end,:,:]
	u_f = uf_init[2:end,:,:]

	# use average value functions v = (r + ψ_m + ψ_f)*V 
	v_m = zeros(u_m)
	v_f = zeros(u_f)

	# NOTE: CartesianIndex doesn't need to be splatted
	# female value function
	for y in CartesianRange(size(u_f)) # use trimmed size
		v_f[y] = sum(λ[:,:,:,y] .* dμ[:,:,:,y] .* u_m)
	end
	# male value function
	for x in CartesianRange(size(u_m)) # use trimmed size
		v_m[x] = sum(λ[x,:,:,:] .* dμ[x,:,:,:] .* u_f)
	end

	return (1-β) * v_m, β * v_f
end

"Compute production function."
function compute_production(δ::Real, ψ_m::Array, ψ_f::Array, dμ::Array,
                            s::Array, v_m::Array, v_f::Array)
	f = similar(s) # instantiate production array
	for xy in CartesianRange(size(f))
		x = xy.I[1:3]
		y = xy.I[4:6]

		f[xy] = s[xy] + v_m[x...] + v_f[y...] - δ * dμ[xy]
	end
	return f
end
