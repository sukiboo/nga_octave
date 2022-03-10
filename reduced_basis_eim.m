# this function constructs a reduced basis of cardinality M for the given training set F 
# by using the Empirical Interpolation Method (EIM)
# it returns the reduced basis RB and the runtime
function [RB, ctime] = reduced_basis_eim(F, M)
start_time_cpu = cputime;
start_time_wall = tic;

# initialize the parameters
[N_h,N_tr] = size(F);
RB = zeros(N_h,M);
B = eye(M);
IP = zeros(1,M);
ctime = zeros(M,1);

### the first iteration of EIM
	# choose the first elements of RB and IP and update parameters
	[val_max,ind_max] = max(norm(F,'inf','columns'));
	[val_max,IP(1)] = max(abs(F(:,ind_max)));
	RB(:,1) = F(:,ind_max) / F(IP(1),ind_max);
	ctime(1) = cputime-start_time_cpu;
printf('EIM completed  1 steps after %5.2f sec, e_eim = %.2e\n', toc(start_time_wall), val_max);

### other iterations of EIM
for m = 2 : M
	# approximate each element of the training set F
	I = RB(:,1:m-1) * (B(1:m-1,1:m-1) \ F(IP(1:m-1),:));
	# choose the next elements of RB and IP and update parameters
	[val_max,ind_max] = max(norm(F-I,'inf','columns'));
	[val_max,IP(m)] = max(abs(F(:,ind_max) - I(:,ind_max)));
	RB(:,m) = (F(:,ind_max) - I(:,ind_max)) / (F(IP(m),ind_max) - I(IP(m),ind_max));
	B(m,1:m-1) = RB(IP(m),1:m-1);
	ctime(m) = cputime-start_time_cpu;
printf('EIM completed %2d steps after %5.2f sec, e_eim = %.2e\n', m, toc(start_time_wall), val_max);
endfor

# measure the runtime
printf('EIM generated %d basis elements after %.2f / %.2f sec\n\n', M, cputime-start_time_cpu, toc(start_time_wall));
endfunction