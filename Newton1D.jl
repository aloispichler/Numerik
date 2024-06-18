#=
	Newton method in 1D
	Find a solution, with ultimate accuracy, of the function f(⋅)= 0 involving its explicit derivatife f'(⋅)
	created: 2023, September
	author©: Alois Pichler
=#

using LinearAlgebra
#	╭──	Newton iteration	────────────────────────────────────────────────────
#	│	solve the equation fun(⋅)=0 with one unknown
#	│	the function fun returns (f= …, fD= …)
function Newton(fun::Function, x0::Float64; maxEval= 1000, εAccuracy= 1e-7)
	improvementFound= true
	xMin= x0; fMin= fun(x0); evalCount= 1	# function call
	while fMin.f≠ 0 && (improvementFound || abs(fMin.f)> εAccuracy) && evalCount < maxEval	# run until no improvement found
		direction= fMin.f/ fMin.fD
		if isnan(direction) || !isfinite(direction)
			direction= 1+ abs(x0)
			α= sqrt(eps(Float64)); λ= -1.2	# increasing stepsize and alternate direction
		else
			α= 1.0; λ= .5					# stepsize reduction
		end
		improvementFound= started= false
		while !improvementFound && !(started && x0 == xMin)	&& evalCount < maxEval # never give up
			x0= xMin- α* direction; evalCount+= 1
			fx= fun(x0) 	# function evaluation
			if abs(fx.f)< abs(fMin.f)	# improvement found
				xMin= x0; fMin= fx; improvementFound= true
			else
				α*= λ; started= true
			end
#			@show evalCount, α, nfMin
		end
	end
	abs(fMin.f) > εAccuracy && @warn "Newton failed to converge: |f($(xMin))|= $(fMin.f)"
	# @info "Newton: $evalCount steps; residual= $fMin"
	return (x= xMin, f= abs(fMin.f), evalCount= evalCount)
end
