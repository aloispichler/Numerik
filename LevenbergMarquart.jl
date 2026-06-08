"""
	Levenberg-Marquart algorithm
	Find a solution, with ultimate accuracy, of the function f(⋅)= 0 involving its explicit derivatife f'(⋅)
	created: 2026, May
	author©: Alois Pichler
"""

using LinearAlgebra, ForwardDiff

#	╭────────────────────────────────────────────────────────────────
#	│	dispatch, if Jacobian is not provided …
function LevenbergMarquart(fun::Function, x0::Vector{Float64}; maxEval= 1000, εAccuracy= 1e-7)
	# provide funD, the derivative (Jacobian) of fun
	Jacobian = x -> ForwardDiff.jacobian(fun, x)
	return LevenbergMarquart(fun, Jacobian, x0; maxEval= maxEval, εAccuracy=εAccuracy)
end

#	╭────────────────────────────────────────────────────────────────
#	│	Levenberg-Marquart iteration
function LevenbergMarquart(fun::Function, funD::Function, x0::Vector{Float64}; maxEval= 1000, εAccuracy= 1e-7)
	evalCount= 0
	improvementFound= true; direction= Vector{Float64}(undef, length(x0))
	grad= Vector{Float64}(undef, length(x0))
	xMin= copy(x0); fMin= fun(x0); nfMin= norm(fMin); λ::Float64= 0.0
	Df2= Array{Float64}(undef, length(x0), length(x0))
	while (improvementFound || nfMin > εAccuracy) && nfMin > 0.0 && evalCount < maxEval	# run until no improvement found
		# @show xMin
		# @show fMin
		if improvementFound
			Df= funD(xMin); Df2= Df'* Df
			grad= Df'* fMin	# gradient of ½‖f(xMin)‖², steepest ascent
		end
		λ0= 1e-10+ 1e-10* maximum(diag(Df2))
		if λ < λ0; λ= λ0; end 	# set initial regularization
		direction= (Df2 + λ* I) \ grad; xTrial= xMin - direction
		if xMin == xTrial
			@info "LevenbergMarquart: no improvement ($(evalCount)): trying random direction …" 
			direction= randn(length(x0))* (1e-7 + λ)
			xTrial= xMin - direction
		end
		tmpLinear= grad'* direction; tmpQuadratic= direction'* Df2* direction/ 2
		tmpDen= tmpLinear - tmpQuadratic
		fTrial= fun(xTrial); nfTrial= norm(fTrial); evalCount+= 1
		if nfTrial < nfMin # && tmpDen > 0
			xMin= xTrial; improvementFound= true
			fMin= fTrial; tmpNom= (nfMin^2 - nfTrial^2)/ 2; nfMin= nfTrial
			if tmpDen > 0
				if tmpNom > 0.75* tmpDen	# gain_ratio ρ = tmpNom / tmpDen
					λ/= 3
				elseif tmpNom < 0.25* tmpDen
					λ*= 2
				end
			end
		else # no improvement found. Try gradient descent
			λ*= 4; improvementFound= false
		end
		# @show evalCount, λ, nfMin #, xMin
	end
	nfMin > εAccuracy && @warn "Levenberg–Marquart failed to converge: ‖f(x$(evalCount))‖= $(nfMin)"
#	@info "Levenberg-Marquart: $evalCount steps; residual= $nfMin"
	return (xMin= xMin, normfMin= nfMin, evalCount= evalCount)
end
