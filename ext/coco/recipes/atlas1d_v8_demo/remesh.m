function [prob stat xtr] = remesh(prob, data, chart, ub, Vb)
%REMESH   Remesh function for equipartitioned interpolation problem.
%
% Implement moving mesh conditions from Sect. 20.1.2 of Recipes for
% Continuation.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: remesh.m 2839 2015-03-05 17:09:01Z fschild $

f  = ub(data.x_idx);   % Extract current mesh values
df = Vb(data.x_idx,:); % Extract current tangent vector and tangent space

g  = 2*(f-f(1))/(f(end)-f(1))-1+data.s*data.t; % Evaluating stretching transformation
u  = 2*(g-g(1))/(g(end)-g(1))-1;  % Evaluate interpolant
t0 = interp1(u, data.t, data.th); % Equipartition
t0([1 end]) = [-1 1];             % Fix boundaries
ua = [interp1(data.t, f,  t0); ub(data.p_idx)];   % Updated mesh values 
Va = [interp1(data.t, df, t0); Vb(data.p_idx,:)]; % Updated tangent vector and tangent space

xtr = data.xtr; % Index array of invariant elements
N   = numel(t0);
if numel(data.t)~=numel(t0) % If mesh changed
  data.N     = N;
  data.x_idx = 1:N;         % Updated index array of mesh values
  data.p_idx = N+data.pdim; % Updated index for problem parameter. Bug: should be N+1:N+data.pdim
  xtr(end-data.pdim:end) = N:N+data.pdim; % Translation table for end point and problem parameters 
  data.xtr   = zeros(N+data.pdim,1);
  data.xtr([1 N:N+data.pdim]) = [1 N:N+data.pdim]; % New index array of invariant elements
end
data.t = t0; % New mesh
prob   = coco_change_func(prob, data, 'u0', ua, 'vecs', Va); % Update function data structure and solution

H  = max(abs(diff(data.t))); % Largest mesh interval
N2 = N;
if H>data.HINC % If above upper bound on adaptation window
  N2 = min(100, ceil(N*min((H/data.HINC), 1.1)));
elseif H<data.HDEC % If below lower bound on adaptation window
  N2 = max(10, ceil(N*max((H/data.HDEC), 0.75)));
end
if N~=N2 % If number of mesh points has changed
  data.th = linspace(-1, 1, N2)'; % Update reference mesh
  stat = 'repeat';  % Call the remesh function again
else
  stat = 'success'; % Otherwise finish
end

end
