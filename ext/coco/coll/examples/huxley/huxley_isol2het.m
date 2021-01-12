function prob = huxley_isol2het(prob, segs, dev0)
%HUXLEY_ISOL2HET   Construct instance of huxley continuation problem from initial data.
%
% PROB = HUXLEY_ISOL2HET(PROB, SEGS, dev0)
%
% PROB - Continuation problem structure.
% SEGS - Array of initial solution information for two instances of the
%        'coll' toolbox.
% dev0 - Array of deviations of heteroclinic trajectory end points from
%        corresponding equilibria.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: huxley_isol2het.m 2845 2015-05-12 16:52:46Z hdankowicz $

prob = ode_isol2coll(prob, 'huxley1', @huxley, ...
  segs(1).t0, segs(1).x0, segs(1).p0); % Construct first instance of 'coll' toolbox
prob = ode_isol2coll(prob, 'huxley2', @huxley, ...
  segs(2).t0, segs(2).x0, segs(2).p0); % Construct second instance of 'coll' toolbox

prob = huxley_close_het(prob, dev0); % Apply gluing conditions

end
