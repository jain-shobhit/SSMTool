function prob = coll_construct_seg(prob, tbid, data, sol)
%COLL_CONSTRUCT_SEG   Append an instance of 'coll' to problem.
%
% Add collocation and continuity conditions, monitor functions that
% evaluate to the problem parameters, and corresponding inactive
% continuation parameters.
%
% Identical to coll_v1.
%
% PROB = COLL_CONSTRUCT_SEG(PROB, TBID, DATA, SOL)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% SOL  - Initial solution guess.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_construct_seg.m 2839 2015-03-05 17:09:01Z fschild $

prob = coco_add_func(prob, tbid, @coll_F, @coll_DFDU, data, 'zero', ...
  'u0', sol.u);
if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  uidx = coco_get_func_data(prob, tbid, 'uidx'); % Context-dependent index array
  prob = coco_add_pars(prob, fid, uidx(data.p_idx), data.pnames);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

function [data y] = coll_F(prob, data, u)
%COLL_F   Collocation zero problem.
%
% Expects vectorized encoding of vector field as function of two arguments.
%
% Identical to coll_v1.

x = u(data.xbp_idx); % Extract basepoint values
T = u(data.T_idx);   % Extract interval length
p = u(data.p_idx);   % Extract problem parameters

xx = reshape(data.W*x, data.x_shp); % Values at collocation nodes
pp = repmat(p, data.p_rep);

ode = data.fhan(xx, pp);
ode = (0.5*T/data.coll.NTST)*ode(:)-data.Wp*x; % Collocation conditions
cnt = data.Q*x;                                % Continuity conditions

y = [ode; cnt];

end

function [data J] = coll_DFDU(prob, data, u)
%COLL_DFDU   Linearization of collocation zero problem.
%
% Expects vectorized encoding of Jacobians of vector field with respect to
% its arguments. When dfdxhan and/or dfdphan are empty, approximate
% Jacobians are obtained using numerical differentiation.
%
% Identical to coll_v1.

x = u(data.xbp_idx); % Extract basepoint values
T = u(data.T_idx);   % Extract interval length
p = u(data.p_idx);   % Extract problem parameters

xx = reshape(data.W*x, data.x_shp); % Values at collocation nodes
pp = repmat(p, data.p_rep);

if isempty(data.dfdxhan)
  dxode = coco_ezDFDX('f(x,p)v', data.fhan, xx, pp);
else
  dxode = data.dfdxhan(xx, pp);
end
dxode = sparse(data.dxrows, data.dxcols, dxode(:));
dxode = (0.5*T/data.coll.NTST)*dxode*data.W-data.Wp; % W.r.t. basepoint values

dTode = data.fhan(xx, pp);
dTode = (0.5/data.coll.NTST)*dTode(:); % W.r.t. interval length

if isempty(data.dfdphan)
  dpode = coco_ezDFDP('f(x,p)v', data.fhan, xx, pp);
else
  dpode = data.dfdphan(xx, pp);
end
dpode = sparse(data.dprows, data.dpcols, dpode(:));
dpode = (0.5*T/data.coll.NTST)*dpode; % W.r.t. problem parameters

J = [dxode dTode dpode; data.Q data.dTpcnt]; % data.dTpcnt = [0...0]

end
