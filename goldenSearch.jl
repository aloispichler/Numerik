"""
	golden search to locate the minimum of the function fun.
	created: 2021, May
	author©: Alois Pichler

	start the iteration at x1,
	return the optimal point and value.
"""

function goldenSearch(fun::Function, x1::Float64...)
	if length(x1) > 1; x2= x1[2]; x1= x1[1]   # catch the input
	else               x2= x1[1]; x1= x2; end
	# set the start variables
	if abs(x2-x1)< max(1e-9, 1e-4*(abs(x1)+abs(x2)))
		if x1< 0; x2= 1.05* x1- 1e-3
		else      x2= 1.05* x1+ 1e-3; end
	end
	# evaluate the functions
	f1= fun(x1); f2= fun(x2); f3= -Inf; x3= x1
	if isnan(f1) || isnan(f2)
		return (NaN, NaN)
	end
	phi= (5.0^0.5- 1)/ 2  							# golden ratio, phi= 0.618
	if f1 > f2
		x1, x2= x2, x1; f1, f2= f2, f1		# now: f1 <= f2
	end
	while f3 <= f1   # finde zuerst zwei Punkte, wo dazwischen ein Minimum ist.
		x3= (2+ phi)* x1- (1+ phi)* x2				# extrapolate towards minimum (ie. x1) by keeping golden ratio
		if isnan(x3); return x1, f1; end
		f3= fun(x3)
		if isnan(f3); return x1, f1; end
		if f3 <= f1; x2, f2, x1, f1= x1, f1, x3, f3	# maintain the relation f1 < f2
		else
			x1, x2, x3= x3, x1, x2; f1, f2, f3= f3, f1, f2
			if x1 > x3; x1, x3= x3, x1; end
			break
		end
	end
	found= false; deltaX= Inf
	while !found		# mit golden Section search verbessere die Lösung
		tmp= x3- x1		# run to machine epsilon
		if tmp< deltaX; deltaX= tmp
		else; found = true; end
		if x3- x2 > x2- x1; x= phi* x2+ (1- phi)* x3
		else;               x= phi* x2+ (1- phi)* x1; end
		f = fun(x)		# evaluate function
		if isnan(f); return x2, f2; end
		if f< f2
			if x3- x2 > x2- x1;
					x1, x2= x2, x; f1, f2= f2, f
				else
					x3, x2= x2, x; f3, f2= f2, f
			end
			else
			if x3- x2 > x2- x1
					x3= x; f3= f
				else
					x1= x; f1= f
			end
		end
	end
	return (x= x2, fx= f2) 		# return optimal point and optimal value
end