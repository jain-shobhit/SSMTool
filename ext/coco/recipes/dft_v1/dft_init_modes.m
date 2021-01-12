function data = dft_init_modes(data, sol)
%DFT_INIT_MODES   Initialize adaptable data for an instance of 'dft'.
%
% Populate remaining fields of the toolbox data structure used by 'dft'
% function objects. Function may be called for adaptive updates.
%
% DATA = DFT_INIT_MODES(DATA, SOL)
%
% DATA - Toolbox data structure.
% SOL  - Solution structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: dft_init_modes.m 2839 2015-03-05 17:09:01Z fschild $

N    = 2*data.dft.NMAX+2; % Number of mesh points
NMOD = data.dft.NMOD;     % Number of Fourier modes
dim  = data.dim;          % State-space dimension

mdim = dim*(2*NMOD+1);
ddim = dim*(N-2*NMOD-2);
data.mdim   = mdim;       % Active dimensions
data.ddim   = ddim;       % Dummy dimensions
data.xf_idx = 1:data.mdim; % Index array for active modes
data.xd_idx = data.mdim+(1:data.ddim)'; % Index array for dummy modes

phs = repmat(1:NMOD, [dim 1]);
phs = [-phs; phs].*[sol.x0(:,3:2:end); sol.x0(:,2:2:end)];
data.phs = [zeros(dim,1); phs(:)]'; % Reference coefficients for phase condition
row = sparse(1:ddim, 1+mdim:ddim+mdim, ones(1,ddim), ddim, data.T_idx);
data.jac = [row; data.phs, sparse(1,ddim+1)]; % Jacobian of dummy and phase conditions.

D           = 2*pi*1i*diag([0:NMOD zeros(1,N-2*NMOD-2+1) -NMOD:-1]);
Ds          = kron(D, speye(dim));
W           = kron([1, zeros(1,2*NMOD);
  zeros(N-1,1), [kron(speye(NMOD), [1 1i]);
  zeros(N-2*NMOD-1,2*NMOD);
  kron(fliplr(speye(NMOD)), [1 -1i])]], speye(dim));
data.Wp     = Ds(1:dim*(NMOD+1),:)*W;        % Differentiation matrix
data.F      = data.Forig(:,NMOD+2:N/2+1);    % Fourier matrix
data.Fs     = data.Fsorig(1:dim*(NMOD+1),:); % Kronecker tranpose
data.FinvsW = real(data.Finvs*W);            % Interpolation matrix

end
