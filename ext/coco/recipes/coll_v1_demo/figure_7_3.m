function figure_7_3
% Figure 7.3: In contrast to the example in Section 7.3.1, the construction
% of an initial solution for the continuation of a heteroclinic orbit in
% the Huxley dynamical system given by the vector field in Eq. (7.26)
% consists of a sequence of steps. We start with two short orbit segments
% close to each saddle, initial guesses to which can be obtained by Euler
% steps in the respective invariant eigenspace. The two segments after
% initial correction are shown in panel (a). In the subsequent step we grow
% both segments until they terminate on the hyperplane y_1=0.5, (b) and
% (c). At this point, the two end points are separated by a small gap, the
% so-called Lin gap, which we close in the next step, (d). In the last
% step, (e) and (f), we grow the two segments toward the saddle equilibria
% and obtain an approximation to a connecting orbit, which can be used as a
% starting point for subsequent continuation runs.

% Generate data
if coco_exist('huxley1', 'run') && coco_exist('huxley2', 'run') ...
    && coco_exist('huxley3', 'run') && coco_exist('huxley4', 'run') ...
    && coco_exist('huxley5', 'run')
  coco = @(prob, run, varargin) coco_bd_read(run); %#ok<NASGU>
end
run demo_huxley

coco_use_recipes_toolbox coll_v1

%Plot data: panels (a)-(f)
labs = [1, 5, 2, 4, 3, 3];
runs = {'huxley1', 'huxley1', 'huxley2', 'huxley3', 'huxley4', 'huxley5'};

for i=1:6
  figure(i)
  clf
  hold on
  grid on
  box on
  axis([-0.1 1.1 -0.25 0.25])
  
  plot_skeleton(p0, vs, vu, 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.7 0.7 0.7]);
  plot(0.5, 0, 'LineStyle', '-', 'LineWidth', sqrt(2), ...
    'Color', [0.5 0.5 0.5], 'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'white');
  
  sol = coll_read_solution('huxley1', runs(i), labs(i));
  plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
  sol = coll_read_solution('huxley2', runs(i), labs(i));
  plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
  
  hold off
end

coco_use_recipes_toolbox

end

function plot_skeleton(p0, vs, vu, varargin)
vs = sign(vs(1))*vs;
vu = sign(vu(1))*vu;

opts = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6, 'NormControl', 'on');

X = linspace(-0.3,0,300);
[t x] = ode45(@(t,x) huxley(x,p0), [0 11.5], -1.0e-4*vu, opts); %#ok<ASGLU>
Y = interp1([0;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})
[t x] = ode45(@(t,x) -huxley(x,p0), [0 11.5], -1.0e-4*vs, opts); %#ok<ASGLU>
Y = interp1([0;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})

X = linspace(0,0.5,500);
[t x] = ode45(@(t,x) huxley(x,p0), [0 13.5], 1.0e-4*vu, opts); %#ok<ASGLU>
Y = interp1([0;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})
[t x] = ode45(@(t,x) -huxley(x,p0), [0 13.5], 1.0e-4*vs, opts); %#ok<ASGLU>
Y = interp1([0;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})

X = linspace(1,1.3,300);
[t x] = ode45(@(t,x) huxley(x,p0), [0 11.5], [1;0]+1.0e-4*vu, opts); %#ok<ASGLU>
Y = interp1([1;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})
[t x] = ode45(@(t,x) -huxley(x,p0), [0 11.5], [1;0]+1.0e-4*vs, opts); %#ok<ASGLU>
Y = interp1([1;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})

X = linspace(0.5,1,500);
[t x] = ode45(@(t,x) huxley(x,p0), [0 13.5], [1;0]-1.0e-4*vu, opts); %#ok<ASGLU>
Y = interp1([1;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})
[t x] = ode45(@(t,x) -huxley(x,p0), [0 13.5], [1;0]-1.0e-4*vs, opts); %#ok<ASGLU>
Y = interp1([1;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})

end
