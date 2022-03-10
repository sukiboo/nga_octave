# this function calculates the average and maximal errors of approximating elements of the set F 
# by the reduced basis RB in the p-norm with the tolerance eps
function [error_avg, error_max] = test_rb(RB, F, p, eps, step=1)
start_time = tic;

# initialize variables
[N_h,N_tr] = size(F);
[N_h,N_rb] = size(RB);
error = zeros(1,N_tr);
error_avg = zeros(1,N_rb);
error_max = zeros(1,N_rb);
coef = zeros(N_rb,N_tr);

# normalize the reduced basis
RB = RB ./ norm(RB,p,'columns');
# find 'magic' points for RB
Z = magic_points(RB);

### approximate training set with m basis elements
for m = step * [1 : floor(N_rb/step)]
	# find initial guess for approximating coefficients
	coef(1:m,:) = RB(Z(1:m),1:m) \ F(Z(1:m),:);
	# find actual approximating coefficients
	coef(1:m,:) = cell2mat(pararrayfun(nproc,...
					@(n)fminunc_lp(F(:,n), RB(:,1:m), p, coef(1:m,n), eps), ...
					1:N_tr, 'VerboseLevel', 0, 'UniformOutput', 0));
	# calculate approximation error
	error = norm(F - RB(:,1:m)*coef(1:m,:), p,'columns');
	error_avg(m) = mean(error);
	error_max(m) = max(error);
endfor

# report the runtime
printf('Testing completed after %.2f sec\n', toc(start_time));
endfunction