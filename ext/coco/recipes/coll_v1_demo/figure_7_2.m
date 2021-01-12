function figure_7_2
% Figure 7.2: Skeleton of the dynamics of the Huxley dynamical system given
% by the vector field in Eq. (7.26) and its direction field for p_1=1/2 and
% p_2=0. The system has three equilibria: a center at (1/2,0) and two
% saddles that are connected in a heteroclinic cycle. We demonstrate a
% continuation of the upper heteroclinic orbit in Sect. 7.3.2.

% Generate data
if coco_exist('huxley1', 'run') && coco_exist('huxley2', 'run') ...
    && coco_exist('huxley3', 'run') && coco_exist('huxley4', 'run') ...
    && coco_exist('huxley5', 'run')
  coco = @(prob, run, varargin) coco_bd_read(run); %#ok<NASGU>
end
run demo_huxley

% Plot data
N     = 20;
x     = linspace(-0.1,1.1,2*N);
y     = linspace(-0.25,0.25,N);
[X Y] = meshgrid(x,y);
fxy   = huxley([X(:) Y(:)]', repmat([0.5;0], [1 numel(X)]));
fx    = reshape(fxy(1,:), size(X));
fy    = reshape(fxy(2,:), size(X));

figure(1)
clf
hold on
grid on
box on
axis([-0.1 1.1 -0.25 0.25])

quiver(X, Y, fx, fy, 1, 'LineStyle', '-', 'LineWidth',   1, ...
  'Color', [0.3 0.3 0.3], 'ShowArrowHead', 'on',  'MaxHeadSize', 0.2)
plot_skeleton(p0, vs, vu, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black')
plot(0.5, 0, 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerFaceColor', 'white')

hold off

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
