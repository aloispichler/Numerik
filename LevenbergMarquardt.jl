"""
	Levenberg-Marquardt algorithm
	Find a solution – with ultimate accuracy – of the function f(⋅)= 0
	created: 2026, May
	author©: Alois Pichler
"""

# 	To solve the equation f(⋅)= 0 with constraints c(⋅)≤ 0,
# 	consider solving the function [f; √μ c] with (outer)penalty max(0.0, c(⋅))²
#	(or adjusted inner penalty) and adequate penalty weight μ.

using LinearAlgebra, ForwardDiff


#	╭────────────────────────────────────────────────────────────────
#	│	dispatch: provide Jacobian by automatic differentiation, if not provided explicitly …
function LevenbergMarquardt(fun::Function, x0::Vector{Float64}; maxEval= 1000, εAccuracy= 1e-7)
	# Derivative of fun, if not provided
	Jacobian = x -> ForwardDiff.jacobian(fun, x)
	return LevenbergMarquardt(fun, Jacobian, x0; maxEval= maxEval, εAccuracy=εAccuracy)
end

#	╭────────────────────────────────────────────────────────────────
#	│	Levenberg-Marquardt iteration
function LevenbergMarquardt(fun::Function, funD::Function, x0::Vector{Float64}; maxEval= 1000, εAccuracy= 1e-7)
	evalCount= 0
	improvementFound= true; direction= Vector{Float64}(undef, length(x0))
	xMin= copy(x0); fMin= fun(x0); nfMin= norm(fMin); λ::Float64= 0.0
	grad= Vector{Float64}(undef, length(x0))
	Df2= Array{Float64}(undef, length(x0), length(x0))
	while (improvementFound || nfMin > εAccuracy) && nfMin > 0.0 && evalCount < maxEval	# run until no improvement found
		# @show xMin
		# @show fMin
		if improvementFound
			Df= funD(xMin); Df2= Df'* Df
			grad= Df'* fMin	# gradient of ½‖f(xMin)‖², the (-)direction of Gradient Descent
		end
		λ0= 1e-10+ 1e-10* maximum(diag(Df2))
		if λ < λ0; λ= λ0; end 	# set initial regularization
		direction= (Df2 + λ* I) \ grad; xTrial= xMin - direction
		if xMin == xTrial
			@info "LevenbergMarquardt: no improvement ($(evalCount)): trying random direction …" 
			direction= randn(length(x0))* (1e-7 + λ)
			xTrial= xMin - direction
		end
		ReductionPredicted= grad'* direction - direction'* Df2* direction/ 2  # ignore the term f″(x)f(x), because f(x)≈ 0.
		fTrial= fun(xTrial); nfTrial= norm(fTrial); evalCount+= 1
		if nfTrial < nfMin # && ReductionPredicted > 0
			xMin= xTrial; improvementFound= true
			fMin= fTrial; ReductionActual= (nfMin^2 - nfTrial^2)/ 2; nfMin= nfTrial
			if ReductionPredicted > 0
				if ReductionActual > 0.75* ReductionPredicted	# gain_ratio ρ = ReductionActual / ReductionPredicted
					λ/= 3	# prediction was not bad
				elseif ReductionActual < 0.25* ReductionPredicted
					λ*= 2	# prediction was poor
				end
			end
		else # iteration fails, no improvement found. Try Gradient Descent
			λ*= 4; improvementFound= false
		end
		# @show evalCount, λ, nfMin #, xMin
	end
	nfMin > εAccuracy && @warn "Levenberg–Marquardt failed to converge: ‖f(x$(evalCount))‖= $(nfMin)"
#	@info "Levenberg-Marquardt: $evalCount steps; residual= $nfMin"
	return (normfMin= nfMin, evalCount= evalCount, xMin= xMin)
end
