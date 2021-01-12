function sol = dft_init_sol(data, t0, x0, p0)
%DFT_INIT_SOL   Build initial solution guess.
%
% Use sampled trajectory on temporal mesh to construct initial solution
% guess for 'dft' toolbox.
%
% SOL = DFT_INIT_SOL(DATA, T0, X0, P0)
%
% DATA - Toolbox data structure.
% T0   - Array of temporal mesh points.
% X0   - Array of state vectors at mesh points.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: dft_init_sol.m 2839 2015-03-05 17:09:01Z fschild $

N    = 2*data.dft.NMAX+2; % Number of mesh points
NMOD = data.dft.NMOD;     % Number of Fourier modes

t0 = t0(:);
T0 = t0(end)-t0(1); % Interval length
t0 = (t0-t0(1))/T0; % Rescaling - only if T0<>0!
x0 = interp1(t0, x0, (0:N-1)'/N)'; % Collection of basepoint values

x0 = fft(x0.').'; % Fourier transform
x0 = x0(:,1:NMOD+1)/N;
xh = [real(x0(:,2:end)); imag(x0(:,2:end))];
x0 = [x0(:,1) reshape(xh, [size(x0,1) 2*NMOD])]; % Array of Fourier coefficients

sol.x0 = x0;
sol.T0 = T0;
sol.p0 = p0;
sol.u0 = [x0(:); zeros(size(x0,1)*(N-2*NMOD-2),1); T0; p0]; % Pad with zeros

end
