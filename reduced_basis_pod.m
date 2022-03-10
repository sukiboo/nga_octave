# this function constructs a reduced basis of cardinality M for the given training set F 
# by using the Proper Orthogonal Decomposition (POD)
# it returns the reduced basis RB and the runtime
function [RB, ctime] = reduced_basis_pod(F, M)
start_time_cpu = cputime;
start_time_wall = tic;

# initialize parameters
RB = zeros(size(F,1),M);

# perform svd on F
[U,Sigma,V] = svd(F,'econ');
# select reduced basis
RB = U(:,1:M); % = F * V(:,1:M) / Sigma(1:M,1:M)

# measure the runtime
ctime = cputime-start_time_cpu * ones(M,1);
printf('POD generated %d basis elements after %.2f / %.2f sec\n\n', M, cputime-start_time_cpu, toc(start_time_wall));
endfunction