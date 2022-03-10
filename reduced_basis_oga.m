# this function constructs a reduced basis of cardinality M in p-norm for the given training set F 
# by using the Greedy Algorithm (OGA)
# it returns the reduced basis RB and the runtime
function [RB, ctime] = reduced_basis_oga(F, p, M)
pkg load optim
start_time_cpu = cputime;
start_time_wall = tic;

# initialize the parameters
[N_h,N_tr] = size(F);
range = [1:N_tr];
RB = zeros(N_h,M);
coef = zeros(M-1,N_tr);
initial_guess = zeros(1,N_tr);
ctime = zeros(M,1);

### select the first element
[val_max, ind_max] = max(norm(F,p,'columns'));
RB(:,1) = F(:,ind_max) / val_max;
range(ind_max) = [];
ctime(1) = cputime-start_time_cpu;

### run OGA for the given number of iterations
for m = 2 : M
	# find approximation for each element in F
for n = range
	coef(1:m-1,n) = fminunc_lp(F(:,n), RB(:,1:m-1), p, zeros(m-1,1), val_max/100);
endfor
	F(:,range) -= RB(:,1:m-1) * coef(1:m-1,range);
	# find current approximation error
	[val_max,ind_max] = max(norm(F(:,range),p,'columns'));
	# choose RB(m) and update coefficients
	RB(:,m) = F(:,range(ind_max)) / val_max;
	range(ind_max) = [];
	ctime(m) = cputime-start_time_cpu;
	printf('OGA completed %2d steps after %7.2f sec, e_oga = %.2e\n', m, toc(start_time_wall), val_max);
endfor

# report the runtime
printf('OGA generated %d basis elements after %.2f / %.2f sec\n\n', M, cputime-start_time_cpu, toc(start_time_wall));
endfunction