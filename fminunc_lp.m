# this function minimizes the function norm(A-B*x,p) with the tolerance TolFun 
# using gradient descent with the initial guess x_0
function [x_min, f_min] = fminunc_lp(A, B, p, x_0, TolFun=1e-7)

# define initial variables
x_min = x_0;
f_min = norm(A-B*x_min, p);
if f_min == 0
	return
endif

# initial iteration of gradient descent
grad = -B' * (sign(A-B*x_min) .* (abs(A-B*x_min)/norm(A-B*x_min,p)).^(p-1));
dir = -grad;
# line search in the direction dir
[x_diff, f_min] = line_search(@(a) norm(A-B*(x_min+a*dir),p), 0, 2*f_min/norm(B*dir,p), 1e-12);
# update x_min
x_min += x_diff * dir;

# other iterations of gradient decsent
iter = 1;
while (norm(grad) > TolFun && iter < 10)
	grad_old = grad;
	dir_old = dir;
	# calculate gradient at x_min and the search direction
	if norm(A-B*x_min,p) > 1e-16
		grad = -B' * (sign(A-B*x_min) .* (abs(A-B*x_min)/norm(A-B*x_min,p)).^(p-1));
	else
		grad = grad_old;
	endif
	%b = 0;														% without direction update
	%b = (grad' * grad) / (grad_old' * grad_old);				% Fletcher–Reeves direction update
	b = (grad' * (grad - grad_old)) / (grad_old' * grad_old);	% Polak–Ribière direction update
	dir = -grad + b * dir_old;
	# line search in the direction dir
	[x_diff, f_min] = line_search(@(a) norm(A-B*(x_min+a*dir),p), 0, 2*f_min/norm(B*dir,p), 1e-12);
	# update x_min
	x_min += x_diff * dir;
	iter += 1;
endwhile

# complete minimization with the obtained initial guess
opt = optimset('TolX',1e-16, 'TolFun',TolFun);
x_min = fminunc(@(x) norm(A-B*x,p), x_min, opt);

endfunction