# this function constructs a reduced basis of cardinality M in p-norm for the given training set F 
# by using the New Greedy Algorithm (NGA)
# it returns the reduced basis RB, sequence R of the norms of operators R_n, and the cputime
function [RB, R, ctime] = reduced_basis_nga(F, p, M)
start_time_cpu = cputime;
start_time_wall = tic;

# initialize the parameters
[N_h,N_tr] = size(F);
RB = zeros(N_h,M);
ctime = zeros(M,1);

### select the first element
[val_max, ind_max] = max(norm(F,p,'columns'));
RB(:,1) = F(:,ind_max) / val_max;
ctime(1) = cputime-start_time_cpu;

### run NGA for the given number of iterations
for m = 2 : M
	# approximate each element of F by RB(m-1)
	F -= RB(:,m-1) * ((sign(RB(:,m-1)) .* abs(RB(:,m-1)).^(p-1))' * F);
	# find current approximation error
	[val_max, ind_max] = max(norm(F,p,'columns'));
	# choose RB(m) and update coefficients
	RB(:,m) = F(:,ind_max) / val_max;
	ctime(m) = cputime-start_time_cpu;
	printf('NGA completed %2d steps after %5.2f sec, e_nga = %.2e\n', m, toc(start_time_wall), val_max);
endfor

# measure the runtime
printf('NGA generated %d basis elements after %.2f / %.2f sec\n\n', M, cputime-start_time_cpu, toc(start_time_wall));

### find R(m) = max{||R_n||, 1 <= n <= m-1} on subspace V_m
##printf('Measuring norms of operators R_n...\n');
R = ones(1,M);
##
##r = ones(M);
##for n = 1 : M-1
##	# measure the inverse of the norm of operators R_n on subspace V_m
##	printf('Measuring operator %2d out of %2d\n', n+1, M);
##	a_min = zeros(n,1);
##	for m = n+1 : M
##		[a_min,r(n,m)] = fminunc(@(a) norm(RB(:,1:m)*a,p) / norm(RB(:,n+1:m)*a(n+1:m),p), [a_min;1]);
##	endfor
##endfor
### calculate the maximum of operator norms for each subspace
##R = 1 ./ min(r,[],1);
##printf('\n');

endfunction