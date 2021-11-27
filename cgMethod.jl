"""
	global optimization
	conjuaget gradient method (cg); https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
	return minimal function value for given function plus gradient

	created: 2021, November
	author©: Alois Pichler
"""


function cgMethod(fun::Function, grad::Function, x0::Vector{Float64}; maxEval= 1000)
	count= 0; found= false; found2= true;
	α= 1.0	# current learning rate
	g= Array{Float16}(undef, length(x0))	# Gradient
	xMin= x0; fMin= fun(xMin)		# minimum ever

	while count < maxEval
		gOld= g; g= grad(xMin); # gradient evaluation
		d= -g					# steepest descent search direction
		if found				# Polak–Ribière
			β= g'*(g- gOld)/ (gOld'* gOld); isfinite(β) || (β= 0.0)
@show			β= max(0, β)		# automatic direction reset
			d= -g+ β* d			# cg update
		end
		# line search
		if d'* g ≥ 0		# reset to steepest descent
			d= -g
		end
		if norm(d) ≤ 0.0	# minimum found
			break
		end
		if !found && !found2	# no improvement found, twice in a row
			d= x0 + α* randn(length(x0))	# guess, if d is useless
		end
		iCount= 0; found2= found; found= false; γ= α; x0= xMin; f0= fMin; ArmijoFlag= Bool;
		# backtracking line search
		while count < maxEval && iCount < 10	# not more than 4 steps
			count+= 1; iCount+= 1
			fγ= fun(x0+ γ* d)				# function evaluation
			if fγ ≤ fMin; xMin= x0+ γ* d; fMin= fγ; α= γ; found= true; end
			Armijo= fγ ≤ f0+ 0.0001* γ* d'*g; if iCount== 1; ArmijoFlag= Armijo; end
			if ArmijoFlag
				if !Armijo; break; end
				γ*= 1.3	# increase the stepwidth
			else
				if Armijo; break; end
				γ*= 0.6
			end
		end
#@show		α*= 1.3			# increase α in any case
		@show iCount
	end
	return (xMin= xMin, fMin= fMin, countEval= count)
end