function prob = period_B(u0, n)
%PERIOD_B   Period-n orbit of Henon map: version 2.
%
% The continuation problem structure encoded below in prob consists of a
% family of 2n+2 zero functions in 2n+4 continuation variables,
% corresponding to a redundant encoding of the "end points" of the
% trajectory.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: period_B.m 2839 2015-03-05 17:09:01Z fschild $

prob = coco_prob();
prob = coco_add_func(prob, 'henon_1', @henon, [], 'zero', ...
  'u0', u0(1:6));
for i=2:n
  prob = coco_add_func(prob, sprintf('henon_%d', i), @henon, [], ...
    'zero', 'uidx', [1; 2; 2*i+1; 2*i+2], 'u0', u0(2*i+3:2*i+4));
end
prob = coco_add_glue(prob, 'glue', [3; 4], [2*n+3; 2*n+4]);

end
