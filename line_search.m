# this function returns the point of minimum and the minimal value of a convex 
# function f on the interval [a,b] with the given tolerance
function [min, f_min] = line_search(f, a, b, TolX)

# estimate iteration limit to avoid infinite loops due to rounding error
iteration_limit = log2((b-a)/TolX);
iter = 0;

# evaluate f at 0 and 1
f_000 = f(a);
f_100 = f(b);
f_050_known = 0;

# repeat until the sufficient accuracy or the iteration limit is achieved
while (b-a > TolX) && (iter < iteration_limit)
	# evaluate f at 1/2 if it is not known
	if f_050_known == 0
		f_050 = f(.5*a+.5*b);
	endif

	# line search
	if f_000 <= f_100
		if f_000 <= f_050
			f_050_known = 0;
			b = .5*a+.5*b;
			f_100 = f_050;
		else
			f_050_known = 1;
			# evaluate f at 1/4
			f_025 = f(.75*a+.25*b);
			if f_025 <= f_050
				b = .5*a+.5*b;
				f_100 = f_050;
				f_050 = f_025;
			else
				# evaluate f at 3/4
				f_075 = f(.25*a+.75*b);
				if f_075 <= f_050
					a = .5*a+.5*b;
					f_000 = f_050;
					f_050 = f_075;
				else
					tmp = a; a = .75*a+.25*b; b = .25*tmp+.75*b;
					f_000 = f_025;
					f_100 = f_075;
				endif
			endif
		endif
	else
		if f_100 <= f_050
			f_050_known = 0;
			a = .5*a+.5*b;
			f_000 = f_050;
		else
			f_050_known = 1;
			# evaluate f at 3/4
			f_075 = f(.25*a+.75*b);
			if f_075 <= f_050
				a = .5*a+.5*b;
				f_000 = f_050;
				f_050 = f_075;
			else
				# evaluate f at 1/4
				f_025 = f(.75*a+.25*b);
				if f_025 <= f_050
					b = .5*a+.5*b;
					f_100 = f_050;
					f_050 = f_025;
				else
					tmp = a; a = .75*a+.25*b; b = .25*tmp+.75*b;
					f_000 = f_025;
					f_100 = f_075;
				endif
			endif
		endif
	endif
	iter += 1;
endwhile

# return the minimum point
[f_min, min] = min([f(a), f(.5*a+.5*b), f(b)]);
min = [a, .5*a+.5*b, b](min);

endfunction
