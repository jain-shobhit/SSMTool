function y = hspo_bc(data, T, x0, x1, p)
%HSPO_BC   Multi-segment periodic boundary conditions.
%
% Trajectory segments terminate on event surface and connect in circular
% order through resets.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: hspo_bc.m 2839 2015-03-05 17:09:01Z fschild $

y = [];
for i=1:data.nsegs
  y = [y;
    data.hhan(x1(data.x1_idx{i}), p, data.events{i}); ... % events
    x0(data.x0_idx{mod(i,data.nsegs)+1})- ...
    data.ghan(x1(data.x1_idx{i}), p, data.resets{i})];    % resets
end

end
