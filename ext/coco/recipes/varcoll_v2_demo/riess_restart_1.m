function prob = riess_restart_1(prob, run, lab)
%RIESS_RESTART_1   Append heteroclinic orbit problem constructed from stored data.
%
% Construct an instance of 'po', append the corresponding variational
% zero problem, add segments in the stable manifold of the periodic orbit
% and the unstable manifold of the equilibrium at the origin, and introduce
% appropriate eigenspace and boundary conditions.
%
% PROB = RIESS_RESTART_1(PROB, RUN, LAB)
%
% PROB - Continuation problem structure.
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: riess_restart_1.m 2839 2015-03-05 17:09:01Z fschild $

prob = povar_sol2orb(prob, '', run, lab);    % Construct orbit with variational zero problem
prob = coll_sol2seg(prob, 'col1', run, lab); % Construct segment in unstable manifold of origin
prob = coll_sol2seg(prob, 'col2', run, lab); % Construct segment in stable manifold of orbit

[data sol] = coco_read_solution('riess_save_1', run, lab); % Extract problem-specific function data
eps = sol.x(data.eps_idx); % Extract array of distances of end points to orbit and equilibrium
vec = sol.x(data.vec_idx); % Extract eigenvector
lam = sol.x(data.lam_idx); % Extract Floquet multiplier

prob = riess_close_het_1(prob, data, vec, lam, eps); % Append eigenspace and boundary conditions

end
