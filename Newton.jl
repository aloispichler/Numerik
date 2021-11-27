"""
	Newton method
	Find a solution, with ultimate accuracy, of the function f(⋅)= 0 involving its explicit derivatife f'(⋅)
	created: 2021, November
	author©: Alois Pichler
"""

using LinearAlgebra
#	│	Newton iteration
#	╰───────────────────────────────────────────────────────────────
function Newton(fun::Function, funD::Function, x0::Vector{Float64}; maxEval= 1000, εAccuracy= 1e-7)
	count= 0; improvementFound= true; direction= Vector{Float64}(undef, length(x0))
	xMin= x0; fMin= fun(x0); nfMin= norm(fMin)
	while (improvementFound || nfMin >εAccuracy) && count < maxEval	# run until no improvement found
		count+= 1				# count evaluations of derivatives
		if improvementFound		# Newton's direction
			direction= funD(xMin) \ fMin
		else					# nothing found: guess a direction
			direction= randn(length(x0))* (1e-7 + norm(direction))
			@info "Now trying a random direction…" 
		end
		α= 1.0; improvementFound= false
		while !improvementFound && (α > 0.6 || x0 ≠ xMin)	# never give up
			x0= xMin- α* direction	# handle NaN, Inf
			fx= fun(x0)			# function evaluation
			normf= norm(fx)
			if !all(isfinite.(x0)) || !isfinite(normf)
				@info "NaN/ Inf encountered."; break
			end
			if normf < nfMin	# improvement found
				xMin= x0; fMin= fx; nfMin= normf; improvementFound= true
			else				# half Newton step
				α/= 2
			end
#			@show count, α, nfMin
		end
	end
	nfMin > εAccuracy && @warn "Newton failed to converge: ‖f(xMin)‖= $(nfMin)"
	return (xMin= xMin, normfMin= nfMin, count= count)
end