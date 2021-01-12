function prob = var_coll_add(prob, segoid)
%VAR_COLL_ADD   Attach 'varcoll' instance to segment object.
%
% Append monitor function to 'coll' instance that stores the fundamental
% solution to the corresponding variational equations in the 'varcoll'
% instance toolbox data.
%
% PROB = VAR_COLL_ADD(PROB, SEGOID)
%
% PROB   - Continuation problem structure.
% SEGOID - Segment object instance identifier.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: var_coll_add.m 2839 2015-03-05 17:09:01Z fschild $

data = var_coll_init_data(prob, segoid);            % Build toolbox data
uidx = coco_get_func_data(prob, data.tbid, 'uidx'); % Extract 'coll' index array
tbid = coco_get_id(segoid, 'var');
prob = coco_add_func(prob, tbid, @var_coll_seg, data, ...
  'regular', {}, 'uidx', uidx); % Empty range = no continuation parameter

end

function data = var_coll_init_data(prob, segoid)
%VAR_COLL_INIT_DATA   Initialize toolbox data for an instance of 'varcoll'.
%
% Expects vectorized encoding of vector field as function of two arguments.

data.tbid   = coco_get_id(segoid, 'coll');
fdata       = coco_get_func_data(prob, data.tbid, 'data'); % Extract 'coll' toolbox data

dim         = fdata.dim;
data.dim    = dim;                                            % State-space dimension
data.M1_idx = fdata.xbp_idx(end-dim+(1:dim));                 % Jacobian of time-T flow
data.row    = [eye(dim), zeros(dim, fdata.xbp_idx(end)-dim)]; % Initial condition and rhs

end

function [data y] = var_coll_seg(prob, data, u)
%VAR_COLL_SEG  Compute fundamental solution of variational problem.
%
% Solve the variational equations with initial condition given by the
% identity matrix and store the result in the toolbox data structure.
% Return an empty array in the y output argument.

fdata = coco_get_func_data(prob, data.tbid, 'data'); % Extract 'coll' toolbox data

x = u(fdata.xbp_idx); % Extract basepoint values
T = u(fdata.T_idx);   % Extract interval length
p = u(fdata.p_idx);   % Extract problem parameters

xx = reshape(fdata.W*x, fdata.x_shp); % Values at collocation nodes
pp = repmat(p, fdata.p_rep);

if isempty(fdata.dfdxhan)
  dxode = coco_ezDFDX('f(x,p)v', fdata.fhan, xx, pp);
else
  dxode = fdata.dfdxhan(xx, pp);
end
dxode = sparse(fdata.dxrows, fdata.dxcols, dxode(:));
dxode = (0.5*T/fdata.coll.NTST)*dxode*fdata.W-fdata.Wp; % Variational equation

data.M = [data.row; dxode; fdata.Q]\data.row'; % Fundamental solution

y = []; % Empty range = no continuation parameter

end
