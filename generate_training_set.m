# this function samples input data from parametric family of functions
# it returns the semi-normalized training set F_tr
function F_tr = generate_training_set(F,X,D,p,N_h,N_tr)
pkg load linear-algebra

if nargin ~= 6
	# return an error if an unexpected number of arguments is provided
	printf('\nUnexpected number of input arguments -- cannot generate input data.\n\n');
	F_tr = [];
	return
else
	# check that dimensions of the parameters agree
	if (nargin(F) ~= length(N_h) + length(N_tr)) || (size(X,1) < length(N_h)) || (size(D,1) < length(N_tr))
		printf('\nDimensions of input parameters do not agree -- cannot generate input data.\n\n');
		F_tr = [];
		return
	else
		printf('\nGenerating training set by sampling the parametric function...\n');
		# sample parametric and spatial domains
		D_h = arrayfun(@(n)linspace(D(n,1),D(n,2),N_tr(n)), 1:length(N_tr), 'UniformOutput',0);
		X_h = arrayfun(@(n)linspace(X(n,1),X(n,2),N_h(n)), 1:length(N_h), 'UniformOutput',0);
		# generate the mesh
		mesh = cartprod(D_h{:},X_h{:});
		# generate training set by evaluating the parametric function on the mesh
		F_tr = reshape(F(arrayfun(@(n)mesh(:,n), 1:size(mesh,2), 'UniformOutput',0){:}), prod(N_tr),prod(N_h))';
	endif
endif

# semi-normalize the training set
F_tr /= max(norm(F_tr,p,'columns'));
endfunction