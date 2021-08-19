function y = Fourier_recover(x, F, angs)
% FOURIER_RECOVER This function recovers the torus solution in phase angles
% angs based on the utilized Fourier series expansion.
%
% x: a collection of trajectories at (j-1)*2pi/(2N+1), j=0,...,2N. x should
%    have the size (2N+1)-nt-dim, where nt is the number of time points and dim
%    is the dimension of state x
% F: mapping matrix of size (2N+1)-(2N+1)
% angs: a set of angles where the trajectory is recovered
% y: a colloection of trajectories at angs. The arrangement of y is
%    consistent with x

nangs = numel(angs);
dim   = size(x,3);
nt    = size(x,2);
N     = (size(x,1)-1)/2;
angs  = angs(:);
angs = kron(1:N, angs);
Fr   = [ones(nangs,1) reshape([cos(angs);sin(angs)], [nangs 2*N])];

y = zeros(nangs,nt,dim);
% loop over dim
for i=1:dim
    % calculate trajectory of coefficients
    coeffs = F*x(:,:,i);
    % evalute trajectory at angs
    xr = Fr*coeffs;
    y(:,:,i) = xr;
end

end