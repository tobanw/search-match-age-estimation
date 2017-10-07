### GMM Estimation of Arrival Rates ###

# Construct meeting rate array

"Compute uniform arrival rate array from ζ[1]."
function constant_rate(ζ::Vector)
	return ζ[1] * ones(dim_mar)
end

# logistic pdf kernel, scaled to have peak of 1
logistic(x::Real; m=2.0, s=1/8) = 4 * exp(-(x-m)*s) / (1+exp(-(x-m)*s))^2

"Compute bivariate logistic arrival rate array from parameter vector ζ."
function logistic_rate(ζ::Vector)
	# vector: logistic in age-gap with free mean
	rt = Array{Float64}(dim_mar) # NOTE: n_ages excludes age 25
	for a in 1:n_ages, b in 1:n_ages
		rt[a,:,:,b,:,:] = ζ[1] * logistic(a-b; m=ζ[2], s=ζ[3]) * logistic((a+b)/2; m=ζ[4], s=ζ[5]) # varies in age gap AND combined age of couple
		# ternary operator selects s based on sign of centered age gap
		#rt[a,:,:,b,:,:] = ζ[1] * logistic(a-b; m=ζ[2], s=a-b≤ζ[2]?ζ[3]:ζ[4]) * logistic((a+b)/2; m=ζ[5], s=ζ[6]) # varies in age gap AND combined age of couple
	end
    return rt
end

"Compute interpolated arrival rate array from parameter vector ζ."
function interp_rate(ζ::Vector)
	# linear interpolation on knots: age gap x avg age
	# ζ1: fix global scale (1,1)
	age_knots = ([1, 8, 20, 40],)
	age_itp = Interpolations.interpolate(age_knots, [1, ζ[2:4]...], Interpolations.Gridded(Interpolations.Linear()))

	gap_knots = ([-39, -15, -2, 2, 6, 19, 39],)
	gap_itp = Interpolations.interpolate(gap_knots, [ζ[5:7]..., 1, ζ[8:10]...], Interpolations.Gridded(Interpolations.Linear()))

	rt = Array{Float64}(dim_mar) # NOTE: n_ages excludes age 25
	for a in 1:n_ages, b in 1:n_ages
		rt[a,:,:,b,:,:] = ζ[1] * age_itp[(a+b)/2] * gap_itp[a-b] # varies in age gap AND combined age of couple
	end
    return rt
end


# Moment functions

"Compute truncated α for a given MSA."
function compute_alpha(λ::Array, δ::Array, ψm_ψf::Array, mar_init::Array, um_uf::Array, mar_out::Array)
	α = compute_raw_alpha(λ, δ, ψm_ψf, mar_init, um_uf, mar_out)
	# TODO: alternatively, could use smooth truncator
	return clamp.(α, 1e-3, 1 - 1e-3) # enforce 0 < α < 1
end

function compute_MF(λ::Array, um_uf::Array, α::Array)
	return λ .* um_uf .* α
end

function compute_DF(δ::Array, α::Array, m::Array)
	return δ .* (1 - α) .* m
end

function compute_DF_ind(δ::Array, α::Array, m::Array)
	# DF(x) = δ*∫(1-α(x,y))*m(x,y)dy
	DF = compute_DF(δ, α, m)
	# integrate out opposite sex; get rid of extra dimensions
	DF_m = sum(DF, 4:6)[:,:,:,1,1,1]
	DF_f = sum(DF, 1:3)[1,1,1,:,:,:]
	return DF_m, DF_f
end

"Compute model moments for a given MSA."
function model_moments(λ::Array, δ::Array, ψm_ψf::Array, mar_init::Array, um_uf::Array, mar_out::Array)
	# trying raw alpha to help optimizer avoid solutions with alpha = 1
	α = compute_raw_alpha(λ, δ, ψm_ψf, mar_init, um_uf, mar_out)
	
	m = mar_init[2:end,:,:,2:end,:,:]

	DF_m, DF_f = compute_DF_ind(δ, α, m)
	MF = compute_MF(λ, um_uf, α)

	# NOTE: max_age will be dropped in the estimation because of truncation
	return MF, DF_m, DF_f
end

"Minimum distance loss function per MSA."
function loss_msa(MF::Array, DF_m::Array, DF_f::Array,
				  dMF::Array, dDF_m::Array, dDF_f::Array,
				  wgt_men::Array, wgt_wom::Array, wgt_mar::Array)
	# assumes input arrays are from ages 26-64: already dropped min and max

	# proportion weights sum to 1 for each set of moments
	#sqrt_wgt_mar = sqrt.(wgt_mar) # sqrt of mass products for "linear" weighting
	#prop_wgt_mar = sqrt_wgt_mar / sum(sqrt_wgt_mar)
	prop_wgt_mar = wgt_mar / sum(wgt_mar) # "non-linear" weighting

	prop_wgt_wom = wgt_wom / sum(wgt_wom)
	prop_wgt_men = wgt_men / sum(wgt_men)

	loss_MF = sum(prop_wgt_mar .* (MF - dMF).^2)

	# scale the integrated DF moments so that MF and DF are on the same scale
	# the math simplifies to dividing by the number of moments integrated over
	loss_DFm = 0.5 * sum(prop_wgt_men .* (DF_m - dDF_m).^2) / prod(size(wgt_wom))
	loss_DFf = 0.5 * sum(prop_wgt_wom .* (DF_f - dDF_f).^2) / prod(size(wgt_men))

	return loss_MF, loss_DFm, loss_DFf
end

"Objective function to minimize: distance between model and data moments."
function loss(ζx::Vector, ζd::Vector, ψm_ψf::Array,
			  mar_all::Dict{AbstractString,Array}, um_uf_all::Dict{AbstractString,Array},
			  mar_out_all::Dict{AbstractString,Array},
	          dMF::Dict{AbstractString,Array},
			  dDF_m::Dict{AbstractString,Array}, dDF_f::Dict{AbstractString,Array},
	          wgt_men::Dict{AbstractString,Array}, wgt_wom::Dict{AbstractString,Array},
			  wgt_mar_all::Dict{AbstractString,Array})

	ξ = build_ξ(ζx) # construct ξ
	δ = build_δ(ζd) # construct ξ

	loss_MF, loss_DF = 0, 0 # initialize
	for msa in top_msa
		# reconstitute full λ array from age-gap ξ vector
		# (\sum_x u_m(x))*(\sum_y u_f(y)) = U_m * U_f = \sum_x \sum_y u_m(x)*u_f(y)
		λ = ξ / sqrt(sum(um_uf_all["$msa"]))

		# model_moments returns ages 26-65 (only need age 25 for marriage inflows)
		MF, DF_m, DF_f = model_moments(λ, δ, ψm_ψf, mar_all["$msa"],
									   um_uf_all["$msa"][2:end,:,:,2:end,:,:],
									   mar_out_all["$msa"][2:end,:,:,2:end,:,:])

		# feed trimmed (age 26-64) moments and weights to loss function
		lsM, lsDm, lsDf  = loss_msa(MF[1:end-1,:,:,1:end-1,:,:], DF_m[1:end-1,:,:], DF_f[1:end-1,:,:],
									dMF["$msa"][2:end-1,:,:,2:end-1,:,:],
									dDF_m["$msa"][2:end-1,:,:], dDF_f["$msa"][2:end-1,:,:],
									wgt_men["$msa"][2:end-1,:,:], wgt_wom["$msa"][2:end-1,:,:],
									wgt_mar_all["$msa"][2:end-1,:,:,2:end-1,:,:])
		loss_MF += lsM
		loss_DF += lsDm + lsDf
	end

	# Verbose progress indicator 
	#println(@sprintf("sse: %.5e, δ: %.3g, ζ: ", sse, δ), ζ)

	return loss_MF, loss_DF
end

"Objective function to pass to NLopt: requires vectors for `x` and `grad`."
function loss_nlopt(x::Vector, grad::Vector)
	return sum(loss(x[1:length(ζx_0)], x[end-length(ζd_0)+1:end],
					ψm_ψf, marriages, sng_outer, mar_outflow,
					MF, men_DF, wom_DF,
					men_tot, wom_tot, pop_outer))
end

"Objective function to pass to BlackBoxOptim: requires vector for `x`."
function loss_bbopt(x::Vector)
	return sum(loss(x[1:length(ζx_0)], x[end-length(ζd_0)+1:end],
					ψm_ψf, marriages, sng_outer, mar_outflow,
					MF, men_DF, wom_DF,
					men_tot, wom_tot, pop_outer))
end

"Parameter grid search: Given ζ1, sweep through other parameters and return csv string."
function obj_landscaper(ζx1::Real, dg, #zx2g, zx3g, zx4g, zx5g, zd1g, zd2g, zd3g, zd4g, zd5g,
						ψm_ψf, marriages, sng_outer, mar_outflow, MF, men_DF, wom_DF,
						men_tot, wom_tot, pop_outer)
	res = "" # results string
	#for ζx2 in zx2g, ζx3 in zx3g, ζx4 in zx4g, ζx5 in zx5g,
	#	ζd1 in zd1g, ζd2 in zd2g, ζd3 in zd3g, ζd4 in zd4g, ζd5 in zd5g
	for δ in dg
		lMF, lDF = loss([ζx1], [δ], #ζx2, ζx3, ζx4, ζx5],
						#[ζd1, ζd2, ζd3, ζd4, ζd5],
						ψm_ψf, marriages, sng_outer, mar_outflow,
						MF, men_DF, wom_DF, men_tot, wom_tot, pop_outer)
		#res *= string(ζx1, ",", ζx2, ",", ζx3, ",", ζx4, ",", ζx5, ",", ζd1, ",", ζd2, ",", ζd3, ",", ζd4, ",", ζd5, ",", lMF, ",", lDF, ",", lMF+lDF, "\n")
		res *= string(ζx1, ",", δ, ",", lMF, ",", lDF, ",", lMF+lDF, "\n")
	end
	return res
end
