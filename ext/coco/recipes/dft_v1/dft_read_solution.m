function [sol data] = dft_read_solution(oid, run, lab)
%COLL_READ_SOLUTION   Read 'dft' solution and toolbox data from disk.
%
% Extract data and chart structures associated with 'dft' toolbox instance
% identifier from solution file and construct solution structure.
%
% [SOL DATA] = DFT_READ_SOLUTION(OID, RUN, LAB)
%
% SOL  - Solution (struct).
% DATA - Toolbox data (struct).
% OID  - Object instance identifier (string).
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: dft_read_solution.m 2839 2015-03-05 17:09:01Z fschild $

tbid         = coco_get_id(oid, 'dft');
[data chart] = coco_read_solution(tbid, run, lab);

dim = data.dim;
N   = 2*data.dft.NMAX+2; % Number of mesh points

sol.t = chart.x(data.T_idx)*(0:N)/N;
sol.c = reshape(chart.x(data.xf_idx), [dim data.dft.NMOD*2+1]);
sol.x = reshape(real(data.FinvsW*chart.x(data.xf_idx)), [dim N]); % Inverse Fourier transform
sol.x = [sol.x sol.x(:,1)]'; % Repeat initial point
sol.p = chart.x(data.p_idx);
sol.T = chart.x(data.T_idx);

end
