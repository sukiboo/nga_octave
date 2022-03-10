function main(p, iterations, num_sim, eps)
# this program constructs reduced bases using different algorithms and compares thier performance
pkg load parallel
pkg load symbolic
format short
more off
%%%%%%%%%%%%%%%%%%%%%%%%%%clear all
close all
clc
start_time = tic;

#############################################################################################################
### specify parameters

# choose the method of generating training data: 0 -- random data / 1 -- parametric function
method = 0;

#############################################################################################################
### generate training set

### random data
if method == 0
	# select approximation algorithms
	# 1--OGA, 2--NGA, 3--EIM, 4--POD
	rb_alg = '12';
	# initial parameters
	step_size = 1; grid_step = 5;
%	num_sim = 10;
	N_h = 1000; d = 100; N_tr = 1000;
%	p = 100; iterations = 50; eps = 1e-04;
	F_tr = randn(N_h,d) * (2*rand(d,N_tr) - 1);
	F_tr /= max(norm(F_tr,p,'columns'));

### parametric function
elseif method == 1
	syms x y z t u v
	# select approximation algorithms
	# 1--OGA, 2--NGA, 3--EIM, 4--POD
	rb_alg = '1234';
	step_size = 3;

	# function, spatial and parameter domains
%	F = exp(x+2*t+3*u) * (exp(-pi*abs(x-t/2)) * asin(sin(2*pi*exp(t)*x)) + exp(-pi*abs(x+u/2)) * asin(sin(exp(pi-u)*x)));
%	F = exp(u-abs(x-t/2)) * asin(sin(u*x/(t+1))) + exp(t-abs(x+u/2)) * asin(sin((pi-t)*x/(u+1)));
%	X = [-2,2]; D = [0,2; 0,2]; N_h = 1e+05; N_tr = [32,32];
%	p = 3; iterations = 30; eps = 1e-09;

	F = (1-x) * cos((x+1)*3*pi*t) * exp(-(1+x)*t);
%	F = (1-x) * cos((x+1)*3*pi*t) * exp(-(1+x)*t) + ((x > (t-1)/1000) - (x > (t)/1000));
	X = [-1,1]; D = [1,pi]; N_h = 1e+05; N_tr = 500;
	p = 1; iterations = 30; eps = 1e-16;

%	F = sin(x*t) * cos(y*u);
%	X = [0,1; 0,1]; D = [pi/3,2*pi; pi/3,2*pi]; N_h = [300,300]; N_tr = [25,25];
%	p = 1; iterations = 30; eps = 1e-10;

%	F = sin(x*t) * cos(y*u) * exp(abs(x)*t + abs(y)*u);
%	X = [-pi,pi; -pi,pi]; D = [-pi,pi; -pi,pi]; N_h = [300,300]; N_tr = [25,25];
%	p = 1; iterations = 30; eps = 1e-10;

%	F = sin(x*t+y*u+z*v) * exp((x+t)*(y+u)*(z+v)) * (1-x)*(1-y)*(1-z);
%	F = sin(x*t+y*u+z*v) * exp(x*y*z-t*u*v) * (1-x)*(1-y)*(1-z);
%	X = [0,1; 0,1; 0,1]; D = [0,1; 0,1; 0,1]; N_h = [50,50,50]; N_tr = [8,8,8];
%	p = 1; iterations = 30; eps = 1e-10;

	# sample the parametric function
	F_tr = generate_training_set(function_handle(F),X,D,p,N_h,N_tr);
%	F_tr += max(abs(F_tr(:))) * sprand(prod(N_h),prod(N_tr),.00001);
%	F_tr += mean(abs(F_tr(:))) * sprand(prod(N_h),prod(N_tr),.01);
%	F_tr /= max(norm(F_tr,p,'columns'));
%	F_tr2 = F_tr;
endif

# check that the training set is not empty
if length(F_tr) == 0
	printf('\nTraining set is empty -- cannot approximate.\n\n');
	return
else
	printf('Training set cardinality: %d\n', prod(N_tr));
endif

#############################################################################################################
### construct reduced bases

if nnz(rb_alg) > 0
printf('Constructing reduced bases...\n\n');

# orthogonal greedy algorithm (OGA)
if any(rb_alg == '1')
	[RB_oga, ctime_oga] = reduced_basis_oga(F_tr, p, iterations);
endif
# natural greedy algorithm (NGA)
if any(rb_alg == '2')
	[RB_nga, R_nga, ctime_nga] = reduced_basis_nga(F_tr, p, iterations);
endif
# emperical interpolation method (EIM)
if any(rb_alg == '3')
	[RB_eim, ctime_eim] = reduced_basis_eim(F_tr, iterations);
endif
# proper orthogonal decomposition (POD)
if any(rb_alg == '4')
	[RB_pod, ctime_pod] = reduced_basis_pod(F_tr, iterations);
endif
endif

#############################################################################################################
### compare different reduced bases methods

printf('\nTesting constructed reduced bases...\n');
if method == 0
	# initialize variables
	error_rb_avg = zeros(2,iterations,num_sim);
	error_rb_max = zeros(2,iterations,num_sim);

	# test reduced bases
	for n = 1 : num_sim
		[error_rb_avg(1,:,n), error_rb_max(1,:,n)] = test_rb(RB_oga, F_tr, p, eps, step_size);
		[error_rb_avg(2,:,n), error_rb_max(2,:,n)] = test_rb(RB_nga, F_tr, p, eps, step_size);
		# generate training set and reduced basis for the next simulation
		if n < num_sim
			F_tr = randn(N_h,d) * (2*rand(d,N_tr) - 1);
			F_tr /= max(norm(F_tr,p,'columns'));
			[RB_oga, ctime_oga] = reduced_basis_oga(F_tr, p, iterations);
			[RB_nga, R_nga, ctime_nga] = reduced_basis_nga(F_tr, p, iterations);
		endif
	endfor

	# plot approximation error ratio
	ratio_avg(1,:) = min(error_rb_avg(1,:,:) ./ error_rb_avg(2,:,:), [], 3);
	ratio_avg(2,:) = mean(error_rb_avg(1,:,:),3) ./ mean(error_rb_avg(2,:,:),3);
	ratio_avg(3,:) = max(error_rb_avg(1,:,:) ./ error_rb_avg(2,:,:), [], 3);
	ratio_max(1,:) = min(error_rb_max(1,:,:) ./ error_rb_max(2,:,:), [], 3);
	ratio_max(2,:) = mean(error_rb_max(1,:,:),3) ./ mean(error_rb_max(2,:,:),3);
	ratio_max(3,:) = max(error_rb_max(1,:,:) ./ error_rb_max(2,:,:), [], 3);
	plot_ratio(ratio_avg, ratio_max, grid_step);

elseif method == 1
	# get denoised training set
	F_tr = generate_training_set(function_handle(F),X,D,p,N_h,N_tr);
	# initialize variables
	error_rb_avg = zeros(4,iterations);
	error_rb_max = zeros(4,iterations);

	# OGA
	if any(rb_alg == '1')
		[error_rb_avg(1,:), error_rb_max(1,:)] = test_rb(RB_oga, F_tr, p, eps, step_size);
	endif
	# NGA
	if any(rb_alg == '2')
		[error_rb_avg(2,:), error_rb_max(2,:)] = test_rb(RB_nga, F_tr, p, eps, step_size);
	endif
	# EIM
	if any(rb_alg == '3')
		[error_rb_avg(3,:), error_rb_max(3,:)] = test_rb(RB_eim, F_tr, p, eps, step_size);
	endif
	# POD
	if any(rb_alg == '4')
		[error_rb_avg(4,:), error_rb_max(4,:)] = test_rb(RB_pod, F_tr, p, eps, step_size);
	endif

	# plot and report average and maximal errors
	plot_appr_log([error_rb_avg;error_rb_max], step_size);
	# measure and plot algortihm quality
	qual_avg = 1 ./ (error_rb_avg .* [ctime_oga,ctime_nga,ctime_eim,ctime_pod]');
	qual_max = 1 ./ (error_rb_max .* [ctime_oga,ctime_nga,ctime_eim,ctime_pod]');
	plot_qual_log(qual_avg, step_size, 'average');
	plot_qual_log(qual_max, step_size, 'minimal');

format short e
for k = step_size * [1 : floor(iterations/step_size)]
	printf('%d & %.3e & %.3e & %.3e & %.3e \\\\\\hline\n', k, error_rb_avg(1,k), error_rb_avg(2,k), error_rb_avg(3,k), error_rb_avg(4,k));
endfor
for k = step_size * [1 : floor(iterations/step_size)]
	printf('%d & %.3e & %.3e & %.3e & %.3e \\\\\\hline\n', k, error_rb_max(1,k), error_rb_max(2,k), error_rb_max(3,k), error_rb_max(4,k));
endfor

endif

#############################################################################################################

# report the total program's runtime
printf('\nTotal runtime: %.2f sec\n', toc(start_time));

save(strcat('rnd',num2str(p)))
endfunction
