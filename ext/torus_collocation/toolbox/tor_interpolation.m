function y = tor_interpolation(x,angs)
% TOR_INTERPOLATION This function recovers the torus solution based on Fourier
% expansion and then return torus solution with (high) fidelty at angs
%
% Y = TOR_INTERPOLATION(X)
%
% x: a collection of trajectories at (j-1)*2pi/(2N+1), j=0,...,2N+1. x should
%    have the size nt-dim-(2N+2), where nt is the number of time points and dim
%    is the dimension of state x
% angs: interpolated solution at these angles
% y: dim dimensional amplitudes of state vector x

[~,~,nsegs] = size(x);
% construct mapping matrix F
N  = nsegs/2-1;
Th = 2*pi*(0:2*N)/(2*N+1);
Th = kron(1:N, Th');
F  = [ones(2*N+1,1) 2*reshape([cos(Th);sin(Th)], [2*N+1 2*N])]'/(2*N+1);
x  = permute(x, [3,1,2]);
% Fourier expansion to recover
y  = Fourier_recover(x(1:end-1,:,:), F, angs); % num_angs-nt-dim
y  = permute(y, [2,3,1]);

end