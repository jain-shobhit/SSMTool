function prob = riess_start_2(prob, run)
%RIESS_START_2   Append heteroclinic orbit problem with lin gap constructed from stored data.
%
% Construct an instance of 'po', append the corresponding variational zero
% problem, add segments in the stable manifold of the periodic orbit and
% the unstable manifold of the equilibrium at the origin, introduce
% appropriate eigenspace and boundary conditions, and add lin conditions.
%
% PROB = RIESS_START_2(PROB, RUN)
%
% PROB - Continuation problem structure.
% RUN  - Run identifier (string).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: riess_start_2.m 2839 2015-03-05 17:09:01Z fschild $

bd   = coco_bd_read(run);       % Extract bifurcation data
labs = coco_bd_labs(bd, 'ALL'); % Extract all solution labels
endpoints = [];
labels = [];
for lab=labs
  sol       = coll_read_solution('col2', run, lab); % Extract segment in stable manifold of orbit
  endpoints = [endpoints; sol.x(1,:)];              % Collect end points at t=0
  labels    = [labels; lab];
end
sol     = coll_read_solution('col1', run, 1);                % Extract segment in unstable manifold of equilibrium
pt      = repmat(sol.x(end,:), [size(endpoints, 1) 1]);      % Collect end point at t=1
[m1 i1] = min(sqrt(sum((endpoints-pt).*(endpoints-pt), 2))); % Find closest point in stable manifold of orbit

prob = riess_restart_1(prob, run, labels(i1)); % Reconstruct heteroclinic problem using the corresponding solution

vgap        = endpoints(i1,:)-pt(i1,:);
data.gapvec = vgap/norm(vgap);                     % Lin gap vector
vphase      = endpoints(i1+1,:)-endpoints(i1-1,:);
data.vphase = vphase/norm(vphase);                 % Lin phase condition

prob = riess_close_het_2(prob, data); % Append lin conditions

end
