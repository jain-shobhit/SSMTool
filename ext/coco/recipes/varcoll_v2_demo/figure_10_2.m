function figure_10_2
% Figure 10.2: A close-up of the potential approach and exit of W_0^u to
% and from a periodic orbit, as observed in Fig. 10.1. In panel (a), the
% manifold exits spiraling inward, while the exit is outward in (c). This
% correspond to a switch of approach from an orbit inside a stable manifold
% of a periodic orbit to an approach from an orbit outside. This suggest
% that, in between these parameter values, there exists an orbit
% approaching on the stable manifold -- a heteroclinic connection between
% the equilibrium at 0 and a periodic orbit of saddle type.

s  = 10;
b  = 8/3;

% Generate data: panel (a)
r  = 23;
vu = [1-s+sqrt((1-s)^2 + 4*r*s) ; -2*r ; 0 ];
vu = vu/norm(vu);
x0 = 0.001*vu;
opts = odeset('MaxStep', 0.05);
[t x1] = ode45(@(t,x) lorenz(x,[s;r;b]), [0 1.4], -x0); %#ok<ASGLU>
[t x1] = ode45(@(t,x) lorenz(x,[s;r;b]), [0 20], x1(end,:), opts); %#ok<ASGLU>

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
axis([-12 3 -14 -2 16 31])
view([24 18])

plot3(x1(:,1), x1(:,2), x1(:,3), 'LineStyle', '-', 'LineWidth', 1, ...
  'Color', 'black')

hold off

% Generate data: panel (b)
r  = 24;
vu = [1-s+sqrt((1-s)^2 + 4*r*s) ; -2*r ; 0 ];
vu = vu/norm(vu);
x0 = 0.001*vu;
[t x1] = ode45(@(t,x) lorenz(x,[s;r;b]), [0 1.4], -x0); %#ok<ASGLU>
[t x1] = ode45(@(t,x) lorenz(x,[s;r;b]), [0 20], x1(end,:), opts); %#ok<ASGLU>

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
axis([-12 3 -14 -2 16 31])
view([24 18])

plot3(x1(:,1), x1(:,2), x1(:,3), 'LineStyle', '-', 'LineWidth', 1, ...
  'Color', 'black')

hold off

% Generate data: panel (c)
r  = 25;
vu = [1-s+sqrt((1-s)^2 + 4*r*s) ; -2*r ; 0 ];
vu = vu/norm(vu);
x0 = 0.001*vu;
[t x1] = ode45(@(t,x) lorenz(x,[s;r;b]), [0 1.3], -x0); %#ok<ASGLU>
[t x1] = ode45(@(t,x) lorenz(x,[s;r;b]), [0 20], x1(end,:), opts); %#ok<ASGLU>

% Plot data: panel (c)
figure(3)
clf
hold on
grid on
axis([-12 3 -14 -2 16 31])
view([24 18])

plot3(x1(:,1), x1(:,2), x1(:,3), 'LineStyle', '-', 'LineWidth', 1, ...
  'Color', 'black')

hold off

end
