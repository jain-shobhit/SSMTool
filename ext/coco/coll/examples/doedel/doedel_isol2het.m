function prob = doedel_isol2het(prob, segs, epts, dev0, th0, vec0, lam0)
%DOEDEL_ISOL2HET   Construct instance of doedel continuation problem from initial data.
%
% PROB = DOEDEL_ISOL2HET(PROB, SEGS, EPTS, DEV0, TH0, VEC0, LAM0)
%
% PROB - continuation problem structure.
% SEGS - array of initial solution information for two instances of the
%        'coll' toolbox.
% EPTS - array of initial solution information for two instances of the
%        ‘ep’ toolbox.
% DEV0 - array of deviations of heteroclinic trajectory end points from
%        corresponding equilibria.
% TH0  - initial solution guess for angle relative to horizontal
% VEC0 - initial solution guess for eigenvector
% LAM0 - initial solution guess for eigenvalue

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: doedel_isol2het.m 2896 2015-10-05 20:08:38Z hdankowicz $

% Construct two ‘coll’ instances.
prob = ode_isol2coll(prob, 'doedel1', @doedel, @doedel_DFDX, ...
  @doedel_DFDP, segs(1).t0, segs(1).x0, segs(1).p0);
prob = ode_isol2coll(prob, 'doedel2', @doedel, @doedel_DFDX, ...
  @doedel_DFDP, segs(2).t0, segs(2).x0, segs(2).p0);

% Construct two ‘ep’ instances and include a single-vector variational
% problem with the second instance.
prob = ode_isol2ep(prob, 'doedel3', @doedel, @doedel_DFDX, ...
  @doedel_DFDP, epts(1).x0, epts(1).p0);
prob = ode_isol2ep(prob, 'doedel4', @doedel, @doedel_DFDX, ...
  @doedel_DFDP, epts(2).x0, epts(2).p0, '-var', vec0);

% Apply gluing conditions
prob = doedel_close_het(prob, dev0, th0, lam0);

end
