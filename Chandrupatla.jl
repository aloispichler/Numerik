"""
	Root finding algorithm.
	Find a root of the function f(⋅)= 0 in one dimension without derivaties
	created: 2021, October
	author©: Alois Pichler 
"""

#	Expoit: Chandraputlas method http://dx.doi.org/10.1016/S0965-9978(96)00051-8
function Chandrupatla(fun::Function, x1::Float64...)	# solve fun(⋅)=0.
	maxCount= 2000; count= 2
	if length(x1) > 1; x2= x1[2]; x1= x1[1]   # catch the input
	else               x2= x1[1]; x1= x2; end

	if abs(x2-x1)< max(1e-9, 1e-4*(abs(x1)+abs(x2)))
		if x1< 0; x2= 0.95* x1+ 1e-3		# make sure that starting points are distinct
		else      x2= 1.05* x1+ 1e-3; end
	end
	f1= fun(x1); @assert !isnan(f1) "find Root: f(⋅) is not a number at x1."
	f2= fun(x2); @assert !isnan(f2) "find Root: f(⋅) is not a number at x2."

#	part I: find a bracketing pair with opposite signs
	while sign(f1)== sign(f2) && count < maxCount
		flag= false;
		if abs(f2) < abs(f1); x2, x1, f2, f1 = x1, x2, f1, f2; flag= true; end
		if f1 == f2
			λ= 1.3						# fall back to oscillation
		else
			λ= f1/ (f2-f1)				# note that λ > 0: secant method
			λ= min(1.3, max(0.3, λ))	# avoid too small and too large steps
			flag && (λ= max(λ, 1.3))	# force larger steps
		end
		x2= (1+λ)* x1 - λ* x2		# starting from smallest, extrapolate to new point
		f2= fun(x2); count+= 1
	end

#	part II: Chandraputlas method
	λ= 0.5; Δx= Inf
	while abs(x2-x1) < Δx && count < maxCount	# f1 and f2 have opposite signs
		Δx= abs(x2-x1)
		xλ= (1-λ)*x1 + λ*x2
		fλ= fun(xλ); count+= 1
		if sign(fλ) == sign(f1)
			x1, xλ= xλ, x1		# fλ and f2 have opposite sign, replace x1 by new xλ
			f1, fλ= fλ, f1		# xλ ∉ [x1,x2]
		else
			x1, x2, xλ= xλ, x1, x2	# fλ and f1 have opposite sign, replace x2 by new xλ
			f1, f2, fλ= fλ, f1, f2	# xλ ∉ [x1,x2]
		end			# it holds that (x2 ≤ x1 ≤ xλ) || (x2 ≥ x1 ≥ xλ)
		ξ= (x1-x2)/(xλ-x2); φ= (f1-f2)/(fλ-f2)	# it holds that ξ, φ ∈ [0,1].
		if (1-φ)^2 < 1-ξ && φ^2 < ξ 	# if inverse quadratic interpolation is applicable…
			λ= f1 / (f2-f1) * fλ / (f2-fλ) + (xλ-x1)/(x2-x1)*f1/(fλ-f1)*f2/(fλ-f2)
			!isfinite(λ) && (λ= 0.5)			# ensure that λ is a useful number
		else
			λ= f1/ (f1-f2)		# fall back regula falsi
		end
		λ= min(0.9, max(0.1, λ))	# one decimal digit per iteration
	end
	if abs(f2) < abs(f1); x1, f1 = x2, f2; end
	return (x= x1, f= f1, count= count)
end