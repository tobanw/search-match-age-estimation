### GMM Estimation of Arrival Rates ###

# unscaled logistic pdf
logistic(x::Real; m=2.0, s=1./8) = exp(-(x-m)*s) / (1+exp(-(x-m)*s))^2

"Compute age-gap-specific arrival rate vector ξ from parameter vector ζ."
function build_ξ(ζ::Vector)
	# vector: logistic in age-gap with free mean
	xi = Array{Float64}(n_ages,2,2,n_ages,2,2) # NOTE: n_ages excludes age 25
	for a in 1:n_ages, b in 1:n_ages
		# ternary operator selects s based on sign of centered age gap
		xi[a,:,:,b,:,:] = ζ[1] * logistic(a-b; m=ζ[2], s=a-b≤ζ[2]?ζ[3]:ζ[4]) * logistic((a+b)/2; m=ζ[5], s=ζ[6]) # varies in age gap AND combined age of couple
	end
    return xi
end


# Moment functions

function compute_MF(λ::Array, um_uf::Array, α::Array)
	return λ .* um_uf .* α
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
function model_moments(λ::Array, δ::Real, ψm_ψf::Array, mar_init::Array, um_uf::Array)
	# trying raw alpha to help optimizer avoid solutions with alpha = 1
	α = compute_raw_alpha(λ, δ, ψm_ψf, mar_init, um_uf)
	
	m = mar_init[2:end,:,:,2:end,:,:]

	DF_m = compute_DF_m(δ, α, m)
	DF_f = compute_DF_f(δ, α, m)
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
	sqrt_wgt_mar = sqrt.(wgt_mar) # sqrt of mass products for "linear" weighting
	prop_wgt_mar = sqrt_wgt_mar ./ sum(sqrt_wgt_mar)

	prop_wgt_wom = wgt_wom ./ sum(wgt_wom)
	prop_wgt_men = wgt_men ./ sum(wgt_men)

	loss_M = sum(prop_wgt_mar .* (MF .- dMF).^2)

	# scale the integrated DF moments so that MF and DF are on the same scale
	# the math simplifies to dividing by the number of moments integrated over
	loss_Dm = sum(prop_wgt_men .* (DF_m .- dDF_m).^2) / prod(size(wgt_wom))
	loss_Df = sum(prop_wgt_wom .* (DF_f .- dDF_f).^2) / prod(size(wgt_men))

	return loss_M + 0.5 * (loss_Dm + loss_Df)
end

"Objective function to minimize: distance between model and data moments."
function loss(ζ::Vector, δ::Real, ψm_ψf::Array, #θ::Real,
			  mar_all::Dict{AbstractString,Array}, um_uf_all::Dict{AbstractString,Array},
	          dMF::Dict{AbstractString,Array},
			  dDF_m::Dict{AbstractString,Array}, dDF_f::Dict{AbstractString,Array},
	          wgt_men::Dict{AbstractString,Array}, wgt_wom::Dict{AbstractString,Array},
			  wgt_mar_all::Dict{AbstractString,Array})

	ξ = build_ξ(ζ) # construct ξ

	sse = 0.0 # initialize
	for msa in top_msa
		# reconstitute full λ array from age-gap ξ vector and θ
		# (\sum_x u_m(x))*(\sum_y u_f(y)) = U_m * U_f = \sum_x \sum_y u_m(x)*u_f(y)
		λ = ξ ./ sqrt(sum(um_uf_all["$msa"]))

		# model_moments returns ages 26-65 (only need age 25 for marriage inflows)
		MF, DF_m, DF_f = model_moments(λ, δ, ψm_ψf, mar_all["$msa"], um_uf_all["$msa"][2:end,:,:,2:end,:,:])

		# feed trimmed (age 26-64) moments and weights to loss function
		sse += loss_msa(MF[1:end-1,:,:,1:end-1,:,:], DF_m[1:end-1,:,:], DF_f[1:end-1,:,:],
						dMF["$msa"][2:end-1,:,:,2:end-1,:,:],
						dDF_m["$msa"][2:end-1,:,:], dDF_f["$msa"][2:end-1,:,:],
		                wgt_men["$msa"][2:end-1,:,:], wgt_wom["$msa"][2:end-1,:,:],
						wgt_mar_all["$msa"][2:end-1,:,:,2:end-1,:,:])
	end

	# Verbose progress indicator 
	#println(@sprintf("sse: %.5e, δ: %.3g, ζ: ", sse, δ), ζ) #θ: %.3g,|, θ

	return sse
end

"Objective function to pass to NLopt: requires vectors for `x` and `grad`."
function loss_nlopt(x::Vector, grad::Vector)
	return loss(x[1:end-1], x[end],
				ψm_ψf, marriages, sng_conv,
	            MF, men_DF, wom_DF,
	            men_tot, wom_tot, pop_conv)
end

"Objective function to pass to BlackBoxOptim: requires vector for `x`."
function loss_bbopt(x::Vector)
	return loss(x[1:end-1], x[end],
				ψm_ψf, marriages, sng_conv,
	            MF, men_DF, wom_DF,
	            men_tot, wom_tot, pop_conv)
end

"Parameter grid search: Given ζ1, sweep through other parameters and return csv string."
function obj_landscaper(ζ1::Real, z2g, z3g, z4g, z5g, z6g, dg,
						ψm_ψf, marriages, sng_conv, MF, men_DF, wom_DF,
						men_tot, wom_tot, pop_conv)
	res = "" # results string
	for	ζ2 in z2g, ζ3 in z3g, ζ4 in z4g, ζ5 in z5g, ζ6 in z6g, δ in dg, 
		val = loss([ζ1, ζ2, ζ3, ζ4, ζ5, ζ6], δ,
				   ψm_ψf, marriages, sng_conv, MF, men_DF, wom_DF, men_tot, wom_tot, pop_conv)
		res *= string(ζ1, ",", ζ2, ",", ζ3, ",", ζ4, ",", ζ5, ",", ζ6, ",", δ, ",", val, "\n")
	end
	return res
end
