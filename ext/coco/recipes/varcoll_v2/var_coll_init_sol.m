function M0 = var_coll_init_sol(prob, data)
%VAR_COLL_INIT_SOL   Build 'varcoll' initial solution guess.
%
% Use initial solution guess for 'coll' instance to construct initial
% solution guess for fundamental solution.
%
% M0 = VAR_COLL_INIT_SOL(prob, data)
%
% M0   - Initial solution guess for fundamental solution on basepoints.
% PROB - Continuation problem structure.
% DATA - 'varcoll' instance toolbox data.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: var_coll_init_sol.m 2839 2015-03-05 17:09:01Z fschild $

[fdata u0] = coco_get_func_data(prob, data.coll_id, 'data', 'u0'); % 'coll' data and initial solution guess

x = u0(fdata.xbp_idx); % Extract trajectory basepoint values
T = u0(fdata.T_idx);   % Extract interval length
p = u0(fdata.p_idx);   % Extract problem parameters

xx = reshape(fdata.W*x, fdata.x_shp); % Trajectory values at collocation nodes
pp = repmat(p, fdata.p_rep);

if isempty(fdata.dfdxhan)
  dxode = coco_ezDFDX('f(x,p)v', fdata.fhan, xx, pp);
else
  dxode = fdata.dfdxhan(xx, pp);
end
dxode = sparse(fdata.dxrows, fdata.dxcols, dxode(:));
dxode = (0.5*T/fdata.coll.NTST)*dxode*fdata.W-fdata.Wp;

M0 = [data.R; dxode; fdata.Q]\data.R'; % Fundamental solution

end
