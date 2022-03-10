# this function chooses 'magic' points for the given reduced basis RB
# it returns the set MP of points
function MP = magic_points(RB)

# initialize the parameters
[N_h,N] = size(RB);
MP = zeros(1,N);
B = zeros(N);

### choose order-preserving 'magic' points
# first iteration
	[val_max,MP(1)] = max(abs(RB(:,1)));
	B(1,1) = RB(MP(1),1);
# other iterations
for m = 2 : N
	RB(:,m) -= RB(:,1:m-1) * (B(1:m-1,1:m-1) \ RB(MP(1:m-1),m));
	[val_max,MP(m)] = max(abs(RB(:,m)));
	B(m,1:m) = RB(MP(m),1:m);
endfor

endfunction