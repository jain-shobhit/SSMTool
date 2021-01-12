function prob = period_A(u0, n)
%PERIOD_A   Period-n orbit of Henon map: version 1.
%
% The continuation problem structure encoded below in prob consists of a
% family of 2n zero functions in 2n+2 continuation variables. For n=1, the
% continuation problem structure encoded in prob consists of a family of 4
% zero functions in 6 continuation variables, corresponding to a redundant
% encoding of the "end points" of the trajectory.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: period_A.m 2839 2015-03-05 17:09:01Z fschild $

prob = coco_prob();
prob = coco_add_func(prob, 'henon_1', @henon, [], 'zero', ...
  'u0', u0(1:6));
for i=2:n-1
  prob = coco_add_func(prob, sprintf('henon_%d', i), @henon, [], ...
    'zero', 'uidx', [1; 2; 2*i+1; 2*i+2], 'u0', u0(2*i+3:2*i+4));
end
prob = coco_add_func(prob, sprintf('henon_%d', n), @henon, [], ...
  'zero', 'uidx', [1; 2; 2*n+1; 2*n+2; 3; 4]);

end
