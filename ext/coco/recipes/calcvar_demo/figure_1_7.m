function figure_1_7
% Figure 1.7: Graphs of the approximation errors delta_m:=f(a,b)-p_1 of
% polynomial approximants p_1 of degree m as in Eq. (1.35), for the case
% when b=b_+(a), Y=Y_+(a), and a=1.5. In each case, the m collocation nodes
% are indicated by circles. Panel (b) corresponds to the graphs in Fig.
% 1.6. We observe spectral convergence, i.e., ||delta_m||<C M^(-m) for some
% 0<C<infinity and M>1.

% Generate data
NCOLS = (1:6);

fp = @(x,a) [cosh(a*(x+acosh(a)/a))/a sinh(a*(x+acosh(a)/a))];
t0 = linspace(0,1,100)';
x0 = fp(t0,1.5);
Y0 = x0(end,1);

coco_use_recipes_toolbox coll_v3 bvp_v1

for i=1:numel(NCOLS)
  run = sprintf('run2_%d', NCOLS(i));
  if ~coco_exist(run, 'run')
    prob = coco_prob();
    prob = coco_set(prob, 'coll', 'NTST', 1, 'NCOL', NCOLS(i), 'TOL', 10);
    prob = bvp_isol2seg(prob, '', @caty_ode, t0, x0, 'Y', Y0, ...
      @caty_bc, @caty_bc_DFDX);
    coco(prob, run, [], 0, {'Y' 'bvp.seg.coll.err'});
  end
end

% Plot data
for i=1:numel(NCOLS)
  figure(i)
  clf
  hold on
  grid on
  box on
  
  run = sprintf('run2_%d', NCOLS(i));
  sol  = bvp_read_solution('', run, 1);
  th   = (coll_nodes(NCOLS(i))+1)/2;
  y    = lag_interp(sol.x(:,1), sol.t, t0);
  cpts = lag_interp(sol.x(:,1), sol.t, th);
  
  plot(t0, 0*t0, 'LineStyle', '-', 'LineWidth', sqrt(2), ...
    'Color', [0.6 0.6 0.6])
  plot(t0, x0(:,1)-y, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black')
  x1 = fp(th,1.5);
  plot(th, x1(:,1)-cpts, 'LineStyle', 'none', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')
  
  axis('tight')
  lims = get(gca, 'YLim');
  axis([-0.01 1.01 lims(1)-0.01*(lims(2)-lims(1)) ...
    lims(2)+0.01*(lims(2)-lims(1))])
  hold off
end

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
