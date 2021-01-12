function mesh = coll_mesh(int, maps, tmi)
%COLL_MAPS   Initialise data depending on order and mesh distribution.
%
% Initialize properties that depend on the discretization order (=
% number of mesh intervals) and the discretization parameters (= the
% mesh distribution). Structure changes during continuation.
%
% MESH = COLL_MESH(INT, MAPS, TMI)
%
% MESH - Data structure.
% MAPS - Data structure.
% INT  - Data structure.
% TMI  - Temporal mesh.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_mesh.m 2839 2015-03-05 17:09:01Z fschild $

mesh.tmi = tmi; % Temporal mesh

dim  = int.dim;   % State-space dimension
NCOL = int.NCOL;  % Degree of interpolating polynomials
NTST = maps.NTST; % Number of mesh intervals

bpnum  = NCOL+1;        % Number of basepoints per interval
xcndim = dim*NCOL*NTST; % Number of collocation values

ka        = diff(tmi);
mesh.ka   = ka'; % Reference warping coefficients

wts       = repmat(int.wt, [dim NTST]); % Quadrature weights
kas       = kron(ka,ones(dim,NCOL));
mesh.wts1 = wts(1,:);
mesh.kas1 = kas(1,:);
mesh.wts2 = spdiags(wts(:), 0, xcndim, xcndim);
mesh.kas2 = spdiags(kas(:), 0, xcndim, xcndim);

t  = repmat(tmi(1:end-1)/NTST, [bpnum 1]);
tt = repmat((0.5/NTST)*(int.tm+1), [1 NTST]);
tt = t+repmat(ka, [bpnum 1]).*tt;
mesh.tbp = tt(:)/tt(end); % Temporal mesh with duplication

end
