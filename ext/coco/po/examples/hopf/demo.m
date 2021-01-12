%% Continuation of periodic orbits from a Hopf bifurcation
%
% For the dynamical system 
%
%   x1' = x1*(p1 + p2*r - r^2) - x2
%   x2' = x2*(p1 + p2*r - r^2) + x1
%
% with r = x1^2 + x2^2, a family of periodic orbits emanates from a Hopf
% bifurcation of the equilibrium at the origin when p1 = 0. The bifurcation
% is supercritical if p2 < 0 and subcritical if p2 > 0. The numerical
% analysis first locates the Hopf bifurcation using an 'ep' toolbox
% constructor, and then continues the family of periodic orbits using the
% ode_HB2po constructor.

% The figure shows the initial branch of equilibria, a branch of Hopf
% bifurcations, and a several branches of periodic orbits based at a subset
% of the Hopf bifurcations.

%% Initial encoding

% The continuation problem structure encoded below includes two monitor
% functions that evaluate to each of the problem parameters, and the
% corresponding inactive continuation parameters 'p1' and 'p2'. Its
% dimensional deficit equals 0. A one-dimensional family of equilibria is
% obtained by releasing 'p1' and allowing it to vary during continuation.

% Construct equilibrium point zero problem
prob = coco_prob();
prob = ode_isol2ep(prob, '', @hopf, [0;0], {'p1', 'p2'}, [-1; 1]);

fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
  'ep_run');

bd1  = coco(prob, 'ep_run', [], 1, 'p1', [-1 1]);

%% Start continuation along family of Hopf bifurcations

% The encoding below includes two monitor functions that evaluate to each
% of the problem parameters, and the corresponding inactive continuation
% parameters 'p1' and 'p2'. Its dimensional deficit equals -1. A
% one-dimensional family of Hopf bifurcations is obtained by releasing 'p1'
% and 'p2' and allowing both to vary during continuation.

HBlab = coco_bd_labs(bd1, 'HB');

prob = coco_prob();
prob = ode_ep2HB(prob, '', 'ep_run', HBlab);
prob = coco_add_event(prob, 'UZ', 'p2', -2:0.5:3);

fprintf(...
  '\n Run=''%s'': Continue Hopf bifurcations from point %d in run ''%s''.\n', ...
  'ep_HB_run', HBlab, 'ep_run');

bd2  = coco(prob, 'ep_HB_run', [], 1, {'p1' 'p2'}, {[-1 1] [-2.01 3.01]});

%% Start continuation of periodic orbits from the Hopf bifurcation

% The continuation problem structure encoded below includes two monitor
% functions that evaluate to each of the problem parameters, and the
% corresponding inactive continuation parameters 'p1' and 'p2'. Its
% dimensional deficit equals 0. A one-dimensional family of periodic orbits
% is obtained by releasing 'p1' and allowing it to vary during
% continuation.

prob = coco_prob();
prob = ode_HB2po(prob, '', 'ep_run', HBlab);

fprintf(...
  '\n Run=''%s'': Continue periodic orbits from point %d in run ''%s''.\n', ...
  'po_run', HBlab, 'ep_run');

bd3  = coco(prob, 'po_run', [], 1, 'p1', [-1 1]);

%% Start continuation along family of saddle-node bifurcations

% The continuation problem structure encoded below includes two monitor
% functions that evaluate to each of the problem parameters, and the
% corresponding inactive continuation parameters 'p1' and 'p2'. Its
% dimensional deficit equals -1. A one-dimensional family of saddle-node
% bifurcations of periodic orbits is obtained by releasing 'p1' and 'p2'
% and allowing both to vary during continuation.

SNlab = coco_bd_labs(bd3, 'SN');

prob = coco_prob();
prob = ode_po2SN(prob, '', 'po_run', SNlab(1));

fprintf(...
  '\n Run=''%s'': Continue saddle-node bifurcations from point %d in run ''%s''.\n', ...
  'po_SN_run', SNlab(1), 'po_run');

bd4  = coco(prob, 'po_SN_run', [], 1, {'p1', 'p2'}, {[-1 1] [0.001 3]});

%% Start continuation from subset of Hopf bifurcations

% The continuation problem structure encoded below includes two monitor
% functions that evaluate to each of the problem parameters, and the
% corresponding inactive continuation parameters 'p1' and 'p2'. Its
% dimensional deficit equals 0. A one-dimensional family of periodic orbits
% is obtained by releasing 'p1' and allowing it to vary during
% continuation.

labs = coco_bd_labs(bd2, 'UZ');

for lab=labs
  runid = sprintf('po%d_run', lab);
  prob = coco_prob();
  prob = ode_HB2po(prob, '', 'ep_HB_run', lab);
  
  fprintf(...
  '\n Run=''%s'': Continue family of periodic orbits from point %d in run ''%s''.\n', ...
  runid, lab, 'ep_HB_run');

  coco(prob, runid, [], 1, {'p1' 'p2'}, [-1.01 1]);
end
  
%% Graphical representation of stored solutions

figure(1); clf; hold on; grid on; box on
thm = struct('special', {{'EP', 'HB'}});
coco_plot_bd(thm, 'ep_run', 'p1', 'p2', '||x||_2')
thm = struct('special', {{'EP'}});
coco_plot_bd(thm, 'ep_HB_run', 'p1', 'p2', '||x||_2')
thm = struct('special', {{'EP', 'SN'}});
for lab=labs
  coco_plot_bd(thm, sprintf('po%d_run',lab), 'p1', 'p2', '||x||_{2,MPD}')
end
thm = struct('special', {{'EP'}}, 'xlab', 'p_1', 'ylab', 'p_2');
coco_plot_bd(thm, 'po_SN_run', 'p1', 'p2', '||x||_{2,MPD}')
hold off
view(3)
