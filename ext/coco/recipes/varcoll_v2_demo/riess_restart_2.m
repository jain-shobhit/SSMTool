function prob = riess_restart_2(prob, run, lab)
%RIESS_RESTART_2   Append heteroclinic orbit problem constructed from stored data and add lin gap.
%
% Construct an instance of 'po', append the corresponding variational
% zero problem, add segments in the stable manifold of the periodic orbit
% and the unstable manifold of the equilibrium at the origin, introduce
% appropriate eigenspace and boundary conditions, and add lin gap
% conditions.
%
% PROB = RIESS_RESTART_2(PROB, RUN, LAB)
%
% PROB - Continuation problem structure.
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: riess_restart_2.m 2839 2015-03-05 17:09:01Z fschild $

prob = riess_restart_1(prob, run, lab);              % Reconstruct heteroclinic problem
data = coco_read_solution('riess_save_2', run, lab); % Extract problem-specific function data
prob = riess_close_het_2(prob, data);                % Add lin conditions

end
