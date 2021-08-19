function y = tor_amplitude(x)
% TOR_AMPLITUDE This function recovers the torus solution based on Fourier
% expansion and then calcualte the amplitude of this torus based on L_infty
% norm.
%
% Y = TOR_AMPLITUDE(X)
%
% x: a collection of trajectories at (j-1)*2pi/(2N+1), j=0,...,2N+1. x should
%    have the size nt-dim-(2N+2), where nt is the number of time points and dim
%    is the dimension of state x
% y: dim dimensional amplitudes of state vector x

[~,dim,nsegs] = size(x);
% construct mapping matrix F
N  = nsegs/2-1;
Th = 2*pi*(0:2*N)/(2*N+1);
Th = kron(1:N, Th');
F  = [ones(2*N+1,1) 2*reshape([cos(Th);sin(Th)], [2*N+1 2*N])]'/(2*N+1);
x  = permute(x, [3,1,2]);
% Fourier expansion to recover
angs  = linspace(0,2*pi,129); 
tube  = Fourier_recover(x(1:end-1,:,:), F, angs); % 129-nt-dim

y = zeros(dim,1);
for i=1:dim
    xi = tube(:,:,i);
    y(i) = norm(xi(:),'inf');
end

end