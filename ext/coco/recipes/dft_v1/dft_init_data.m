function data = dft_init_data(data, sol)
%DFT_INIT_DATA   Initialize toolbox data for an instance of 'dft'.
%
% Populate remaining fields of the toolbox data structure used by 'dft'
% function objects.
%
% DATA = DFT_INIT_DATA(DATA, SOL)
%
% DATA - Toolbox data structure.
% SOL  - Solution structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: dft_init_data.m 2839 2015-03-05 17:09:01Z fschild $

N = 2*data.dft.NMAX+2; % Number of mesh points

dim  = size(sol.x0,1);
pdim = numel(sol.p0);
data.dim   = dim;                   % State-space dimension
data.pdim  = pdim;                  % Number of problem parameters
data.T_idx = dim*(N-1)+1;           % Integer index for period
data.p_idx = dim*(N-1)+1+(1:pdim)'; % Index set for problem parameters
data.p_rep = [1 N];                 % Shape for vectorization
data.x_shp = [dim N];               % Shape for vectorization

data.Forig  = (1/N)*fft(eye(N));    % Fourier matrix
Finv        = N*ifft(eye(N));       % Inverse Fourier matrix
data.Fsorig = kron(data.Forig, speye(dim)); % Kronecker transpose
data.Finvs  = kron(Finv, speye(dim)); % Kronecker transpose

data.dxrows = repmat(reshape(1:dim*N, [dim N]), [dim 1]); % Row index array for vectorization
data.dxcols = repmat(1:dim*N, [dim 1]); % Column index array for vectorization
data.dprows = repmat(reshape(1:dim*N, [dim N]), [pdim 1]); % Row index array for vectorization
data.dpcols = repmat(1:pdim, [dim N]);  % Column index array for vectorization

data = dft_init_modes(data, sol); % Append adaptable toolbox data

end
