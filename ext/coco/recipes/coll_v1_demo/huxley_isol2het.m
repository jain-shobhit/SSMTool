function prob = huxley_isol2het(prob, segs, eps0)
%HUXLEY_ISOL2HET   Construct instance of huxley continuation problem from initial data.
%
% PROB = HUXLEY_ISOL2HET(PROB, SEGS, EPS0)
%
% PROB - Continuation problem structure.
% SEGS - Array of initial solution information for two instances of the
%        'coll' toolbox.
% EPS0 - Array of deviations of heteroclinic trajectory end points from
%        corresponding equilibria.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: huxley_isol2het.m 2839 2015-03-05 17:09:01Z fschild $

prob = coll_isol2seg(prob, 'huxley1', @huxley, ...
  segs(1).t0, segs(1).x0, segs(2).p0); % Construct first instance of 'coll' toolbox
prob = coll_isol2seg(prob, 'huxley2', @huxley, ...
  segs(2).t0, segs(2).x0, segs(2).p0); % Construct second instance of 'coll' toolbox

prob = huxley_close_het(prob, eps0); % Apply gluing conditions

end
