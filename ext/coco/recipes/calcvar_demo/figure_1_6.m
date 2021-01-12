function figure_1_6
% Figure 1.6: Comparison of the exact solution (gray), obtained from Eq.
% (1.11), with a numerical solution of the form in Eq. (1.35) of degree m=2
% for the case when b=b_+(a), Y=Y_+(a), and a=1.5. The two collocation
% nodes are marked with circles. For increasing order m, we observe
% spectral convergencen; see Fig. 1.7.

% Generate data
fp = @(x,a) [cosh(a*(x+acosh(a)/a))/a sinh(a*(x+acosh(a)/a))];
t0 = linspace(0,1,100)';
x0 = fp(t0,1.5);
Y0 = x0(end,1);

coco_use_recipes_toolbox coll_v3 bvp_v1

if ~coco_exist('run1', 'run')
  prob = coco_prob();
  prob = coco_set(prob, 'coll', 'NTST', 1, 'NCOL', 2, 'TOL', 1);
  prob = bvp_isol2seg(prob, '', @caty_ode, t0, x0, 'Y', Y0, ...
    @caty_bc, @caty_bc_DFDX);
  coco(prob, 'run1', [], 0, {'Y' 'bvp.seg.coll.err'});
end

% Plot data
figure(1)
clf
hold on
grid on
box on
axis([-0.02 1.02 0 4])

plot(t0, x0(:,1), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', [0.5 0.5 0.5], 'Marker', '.', 'MarkerSize', 15)

sol  = bvp_read_solution('', 'run1', 1);
th   = (coll_nodes(2)+1)/2;
y    = lag_interp(sol.x(:,1), sol.t, t0);
cpts = lag_interp(sol.x(:,1), sol.t, th);

plot(t0, y, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black')
plot(th, cpts, 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 9.5, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

coco_use_recipes_toolbox

end

function [nds wts] = coll_nodes(m)

n = (1:m-1)';
g = n.*sqrt(1./(4.*n.^2-1));
J = -diag(g,1)-diag(g,-1);

[w x] = eig(J);
nds   = diag(x);
wts   = 2*w(1,:).^2;

end

function y = lag_interp(c,t,x)

n = numel(t);
mask = [false true(1,n-1)];
y = zeros(size(x));
for i=1:n
  t0 = t(i);
  c0 = c(i);
  tk = t(mask);
  f  = @(tt) c0*ttprod(tt,tk)/prod(t0-tk);
  y = y + f(x);
  mask = circshift(mask,[0 1]);
end

end

function y = ttprod(tt,tk)

y = ones(size(tt));
for i=1:numel(tk)
  y = y.*(tt-tk(i));
end

end
