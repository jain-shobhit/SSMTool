function figure_10_1
% Figure 10.1: Orbits of the Lorenz system given by the vector field in Eq.
% (10.87) starting in the unstable eigenspace of the equilibrum at 0,
% tracing the unstable manifold. The orbits in (a) approach equilibria,
% while the orbits in (b) seem to approach periodic orbits. This approach
% is shown in more detail in Fig. 10.2

s  = 10;
b  = 8/3;

% Generate data: panel (a)
r  = 10;
vu = [1-s+sqrt((1-s)^2 + 4*r*s) ; -2*r ; 0 ];
vu = vu/norm(vu);
x0 = 0.001*vu;
[t x1] = ode45(@(t,x) lorenz(x,[s;r;b]), [0 10], x0); %#ok<ASGLU>
[t x2] = ode45(@(t,x) lorenz(x,[s;r;b]), [0 10], -x0); %#ok<ASGLU>

% Plot data: panel(a)
figure(1)
clf
hold on
grid on
axis([-10 10 -10 10 0 15])
view([30 30])

plot3(x1(:,1), x1(:,2), x1(:,3), 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', 'black')
plot3(x2(:,1), x2(:,2), x2(:,3), 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', 'black')
plot3(0, 0, 0, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', 'x', 'MarkerSize', 9.5)

hold off

% Generate data: panel (b)
r  = 24;
vu = [1-s+sqrt((1-s)^2 + 4*r*s) ; -2*r ; 0 ];
vu = vu/norm(vu);
x0 = 0.001*vu;
[t x1] = ode45(@(t,x) lorenz(x,[s;r;b]), [0 20], x0); %#ok<ASGLU>
[t x2] = ode45(@(t,x) lorenz(x,[s;r;b]), [0 20], -x0); %#ok<ASGLU>

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
axis([-18 18 -25 25 0 40])
view([30 30])

plot3(x1(:,1), x1(:,2), x1(:,3), 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', 'black')
plot3(x2(:,1), x2(:,2), x2(:,3), 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', 'black')
plot3(0, 0, 0, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', 'x', 'MarkerSize', 9.5)

hold off

end
