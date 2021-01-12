function prob = doedel_isol2het(prob, segs, algs, eps0, th0, vec0, lam0)
%DOEDEL_ISOL2HET   Construct instance of doedel continuation problem from initial data.
%
% PROB = DOEDEL_ISOL2HET(PROB, SEGS, ALGS, EPS0, TH0, VEC0, LAM0)
%
% PROB - continuation problem structure.
% SEGS - array of initial solution information for two instances of the
%        'coll' toolbox.
% ALGS - array of initial solution information for two instances of the
%        'alg' toolbox.
% EPS0 - array of deviations of heteroclinic trajectory end points from
%        corresponding equilibria.
% TH0  - initial solution guess for angle relative to horizontal
% VEC0 - initial solution guess for eigenvector
% LAM0 - initial solution guess for eigenvalue

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: doedel_isol2het.m 2839 2015-03-05 17:09:01Z fschild $

prob = coll_isol2seg(prob, 'doedel1', @doedel, @doedel_DFDX, ...
  segs(1).t0, segs(1).x0, segs(1).p0); % Construct first instance of 'coll' toolbox
prob = coll_isol2seg(prob, 'doedel2', @doedel, @doedel_DFDX, ...
  segs(2).t0, segs(2).x0, segs(2).p0); % Construct second instance of 'coll' toolbox

prob = alg_isol2eqn(prob, 'doedel3', @doedel, @doedel_DFDX, ...
  algs(1).x0, algs(1).p0); % Construct first instance of 'alg' toolbox
prob = alg_isol2eqn(prob, 'doedel4', @doedel, @doedel_DFDX, ...
  algs(2).x0, algs(2).p0); % Construct second instance of 'alg' toolbox

prob = doedel_close_het(prob, eps0, th0, vec0, lam0); % Apply gluing conditions

end
