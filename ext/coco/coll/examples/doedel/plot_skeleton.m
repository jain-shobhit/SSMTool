function plot_skeleton(varargin)
%PLOT_SKELETON    Plot heteroclinic connection and vector field
%
% Extract two-segment heteroclinic connection obtained after
% initial construction and graph this together with vector field.

v1 = [-3/sqrt(10); 1/sqrt(10)];
[t, x] = ode45(@(t,x) -doedel(x,[1;1]), [0 4.2], [1;-1]-1.0e-4*v1); %#ok<ASGLU>
plot(x(:,1), x(:,2), varargin{:})
sol = coll_read_solution('doedel1', 'doedel5', 2);
plot([-1 ; sol.xbp(:,1)], [1 ; sol.xbp(:,2)], varargin{:});
sol = coll_read_solution('doedel2', 'doedel5', 2);
plot([sol.xbp(:,1);x(1,1)], [sol.xbp(:,2);x(1,2)], varargin{:});
plot([1 1], [-3 3], varargin{:});

end
