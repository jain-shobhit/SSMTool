function prob = alg_sol2eqn(run, lab)
%ALG_SOL2EQN   Construct 'alg' instance from saved data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% PROB = ALG_SOL2EQN(RUN, LAB)
%
% PROB - Continuation problem structure.
% RUN  - run identifier (string).
% LAB  - solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_sol2eqn.m 2839 2015-03-05 17:09:01Z fschild $

[sol data] = alg_read_solution(run, lab);  % Extract solution and toolbox data from disk
prob       = alg_construct_eqn(data, sol); % Rebuild continuation problem

end
