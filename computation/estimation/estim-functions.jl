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
	# recursively fill layers from boundary inwards
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
function compute_alpha(λ::Real, δ::Real, ψ_m::Array, ψ_f::Array,
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
		α[xy] = (m[xy] * (ρ + δ + ψ_m[x...] + ψ_f[y...] ) - ρ * mar_init[xy]) /
		        (λ * u_m[x...] * u_f[y...] + δ * m[xy])
	end
	return clamp.(α, 1e-6, 1 - 1e-6) # enforce 0 < α < 1
end

"Compute surplus s by inverting α."
function invert_alpha(c1::Array, α::Array)
	return -c1 .* Φ_inv.(1 .- α) # then use shifted c
end

"Compute value functions."
function value_function(λ::Real, dc1μ::Array, um_init::Array, uf_init::Array)
	# trim off age 25
	u_m = um_init[2:end,:,:]
	u_f = uf_init[2:end,:,:]

	# use average value functions v = (r + ρ + ψ_m + ψ_f)V 
	# initialize v: augmented with v^{T+1}=0 for convenience
	v_m = zeros(um_init)
	v_f = zeros(uf_init)

	# female value function
	for y in CartesianRange(size(u_f)) # use trimmed size
		k = y.I[1] # age
		yy = y.I[2:end] # traits
		# NOTE: CartesianIndex y doesn't need to be splatted
		v_f[k,yy...] = ρ * v_f[k+1,yy...] + β * λ * sum(dc1μ[:,:,:,y] .* u_m)
	end
	# male value function
	for x in CartesianRange(size(u_m)) # use trimmed size
		k = x.I[1] # age
		xx = x.I[2:end] # traits
		# NOTE: CartesianIndex x doesn't need to be splatted
		v_m[k,xx...] = ρ * v_m[k+1,xx...] + β * λ * sum(dc1μ[x,:,:,:] .* u_f)
	end

	return v_m[1:end-1,:,:], v_f[1:end-1,:,:] # trim off age T+1
end

"Compute production function."
function compute_production(δ::Real, ψ_m::Array, ψ_f::Array, d::Array, dc1μ::Array,
                            s::Array, v_m::Array, v_f::Array)
	f = similar(s) # instantiate production array

	# initialize v1: augmented with v^{T+1}=0 for convenience
	v1_m = zeros(collect(size(v_m))+[1,0,0]...)
	v1_f = zeros(collect(size(v_f))+[1,0,0]...)

	# v1 == v otherwise
	v1_m[1:end-1,:,:] = v_m
	v1_f[1:end-1,:,:] = v_f

	# interior: all but (T,T)
	for xy in CartesianRange(size(f)) # loop over interior
		# unpack indices: [xy] == [a,x...,b,y...]
		a = xy.I[1]
		x = xy.I[2:3]
		b = xy.I[4]
		y = xy.I[5:6]

		# truncate indices at T to handle boundary (except for T,T case)
		f[xy] = (s[xy] + v_m[a,x...] + v_f[b,y...] - δ * dc1μ[xy]
				 - ρ * d[min(a+1,end),x...,min(b+1,end),y...] * s[min(a+1,end),x...,min(b+1,end),y...] # drop term for (T,T) case
				 - ρ * v1_m[a+1,x...] / (r + ρ + ψ_m[min(a+1,end),x...])
				 - ρ * v1_f[b+1,y...] / (r + ρ + ψ_f[min(b+1,end),y...])) # allow augmented v to attain T+1
	end

	# patch (T,T) terminal case: no ρ (aging stops)
	for xy in CartesianRange(size(f[end,:,:,end,:,:])) # dims are (2,2,2,2)
		x = xy.I[1:2]
		y = xy.I[3:4]
		f[end,x...,end,y...] = s[end,x...,end,y...] + v_m[end,x...] + v_f[end,y...] - δ * dc1μ[end,x...,end,y...]
	end

	return f
end


### SMM Estimation of Arrival Rates ###

# estimation function: inputs are data moments and pop masses
#	* batch call all msa
#	* model moments function: given params and masses, compute alpha and generate model moments
#	* criterion function: given model moments and data moments, compute weighted SSE
#	* pipe into optimizer, return results
#	* Standard Errors? curvature at optimum

"Compute model moments for a given MSA."
function model_moments(λ::Real, δ::Real, ψ_m::Array, ψ_f::Array,
					   mar_init::Array, um_init::Array, uf_init::Array)
	α = compute_alpha(λ, δ, ψ_m, ψ_f, mar_init, um_init, uf_init)
	
	m = mar_init[2:end,:,:,2:end,:,:]
	u_m = um_init[2:end,:,:]
	u_f = uf_init[2:end,:,:]

	MF_m = zeros(u_m)
	DF_m = zeros(u_m)
	MF_f = zeros(u_f)
	DF_f = zeros(u_f)

	for x in CartesianRange(size(MF_m)) # men
		# MF(x) = λ*u_m(x)*∫α(x,y)*u_f(y)dy
		MF_m[x] = λ * u_m[x] * sum(α[x,y] * u_f[y] for y in CartesianRange(size(u_f)))

		# DF(x) = δ*∫(1-α(x,y))*m(x,y)dy
		DF_m[x] = δ * sum((1-α[x,y]) * m[x,y] for y in CartesianRange(size(u_f)))
	end

	for y in CartesianRange(size(MF_f)) # women
		# MF(y) = λ*u_f(y)*∫α(x,y)*u_m(x)dx
		MF_f[y] = λ * u_f[y] * sum(α[x,y] * u_m[x] for x in CartesianRange(size(u_m)))

		# DF(y) = δ*∫(1-α(x,y))*m(x,y)dx
		DF_f[y] = δ * sum((1-α[x,y]) * m[x,y] for x in CartesianRange(size(u_m)))
	end

	return MF_m, DF_m, MF_f, DF_f
end
