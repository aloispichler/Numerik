"""
	Root finding algorithm.
	created: 2021, November
    author©: Alois Pichler
"""

function regulaFalsi(fun::Function, x0::Float64...)	# secant method 
	maxCount= 2000; count= 2; Δx= Inf
	if length(x0) > 1; x1= x0[2]; x0= x0[1]   # catch the input
	else               x1= x0[1]; x0= x1; end 

	if abs(x0-x1)< max(1e-9, 1e-4*(abs(x0)+abs(x1)))
		if x0< 0; x1= 0.95* x0+ 1e-3
		else      x1= 1.05* x0+ 1e-3; end
	end
	# evaluate the function at starting values
	f0= fun(x0); f1= fun(x1)
	@assert !isnan(f0) && !isnan(f1) "find Root: f(⋅) is not a number at x0." 

	#	solve
	while abs(x1-x0) < Δx && count < maxCount
		count+= 1; if sign(f0) ≠ sign(f1); Δx= abs(x1-x0); end	# watch improvements
		λ= f0/ (f0-f1); if f0==f1 && isodd(count); λ= -λ end	# secant method
		if     isnan(λ); λ= 0.5
		elseif λ > 1; λ= min( 2.1, max(λ,  1.1)) 	# avoid 1 and too large steps
		elseif λ < 0; λ= max(-1.1, min(λ, -0.1))	# avoid 0 and too large backward steps
		else   λ= max(0.1, min(λ, 0.9)); end		# avoid 0 and 1
		xλ= (1-λ)* x0+ λ* x1
		fλ= fun(xλ)
		if sign(fλ) == sign(f1)
			x1= xλ; f1= fλ	# x0 remains, as f0 is smaller
		else
			x0= xλ; f0= fλ; # auswechseln
		end
		if abs(f1) < abs(f0); x1, x0, f1, f0 = x0, x1, f0, f1; end
	end
	return (x= x0, f= f0, count= count)
end