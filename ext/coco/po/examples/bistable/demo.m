%% Frequency response diagram for a hardening oscillator
%
% We study the effects of nonlinear hardening on a harmonically excited
% oscillator, exploring the existence of bistability in the frequency
% response diagram for sufficiently large excitation amplitudes. The
% analysis includes a computation of the nonlinear normal mode associated
% with the free response, in the absence of excitation, emanating from the
% fundamental linear resonance frequency in the limit of zero amplitude.

% Figure shows several families of periodic orbits under variations in
% excitation frequency for fixed values of the excitation amplitude, a
% family of saddle-node bifurcations that correspond to the onset and
% termination of bistability, and a backbone curve of nonlinear normal
% modes for the unforced dynamics.

%% Encoding

% The continuation problem structure encoded below includes three monitor
% functions that evaluate to the problem parameters, and the corresponding
% inactive continuation parameters 'T', 'A', and 'd', representing the
% period and amplitude of excitation, and the damping, respectively. After
% the application of gluing conditions constraining the period of the
% periodic orbit to equal the period of the excitation, the dimensional
% deficit equals -1, since the period corresponds to an initially inactive
% continuation parameter for a non-autonomous vector field. A
% one-dimensional family of periodic orbits is obtained by releasing
% 'po.period' and the period of exctiation and allowing these to vary
% during continuation. The corresponding Floquet multipliers are stored
% with the bifurcation data.

p0 = [2*pi; 0.015; 0.04];
[~, x0] = ode45(@(t,x) bistable(t,x,p0), [0 500*pi], [0; 1]); % Transients
[t0, x0] = ode45(@(t,x) bistable(t,x,p0), [0 2*pi], x0(end,:)'); % Approximate periodic orbit
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
funcs = {@bistable, @bistable_dx, @bistable_dp, @bistable_dt};
coll_args = [funcs, {t0, x0, {'T' 'A' 'd'}, p0}];
prob = ode_isol2po(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_glue(prob, 'glue', uidx(maps.T_idx), uidx(maps.p_idx(1)));
prob = po_mult_add(prob, 'po.orb');  % Store Floquet multipliers with bifurcation data
cont_args = {1, {'po.period' 'T'}, [2*pi/1.3 2*pi/0.7]};

fprintf('\n Run=''%s'': Continue primary family of periodic orbits.\n', ...
  'freq_resp');

bd1  = coco(prob, 'freq_resp', [], cont_args{:});

%% Continue curve of saddle-node bifurcations

% After imposition of gluing conditions, the continuation problem encoded
% below has a dimensional deficit of -2. A one-dimensional family of
% saddle-node bifurcations of periodic orbits is obtained by releasing
% 'po.period', 'T', and 'A', and allowing these to vary during
% continuation.

labs = coco_bd_labs(bd1, 'SN');

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 25);
prob = ode_SN2SN(prob, '', 'freq_resp', labs(1));
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_glue(prob, 'glue', uidx(maps.T_idx), uidx(maps.p_idx(1)));
prob = coco_set(prob, 'cont', 'NAdapt', 5);
prob = coco_add_event(prob, 'UZ', 'A', 0.01:0.005:0.1);
cont_args = {1, {'po.period' 'T' 'A'}, [2*pi/1.3 2*pi/0.7]};

fprintf(...
  '\n Run=''%s'': Continue saddle-node bifurcations from point %d in run ''%s''.\n', ...
  'saddle-node', labs(1), 'freq_resp');

bd2  = coco(prob, 'saddle-node', [], cont_args{:});

%% Sweep frequency-response curves

% For each of a sample of points on the curve of saddle-node bifurcations,
% the continuation problem encoded below has dimensional deficit -1.  A
% one-dimensional family of periodic orbits is obtained by releasing
% 'po.period' and 'T', and allowing these to vary during continuation.

labs  = coco_bd_labs(bd2, 'UZ');

for lab=labs
  prob = coco_prob();
  prob = coco_set(prob, 'ode', 'autonomous', false);
  prob = coco_set(prob, 'po', 'bifus', 'off');
  prob = ode_po2po(prob, '', 'saddle-node', lab);
  [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
  maps = data.coll_seg.maps;
  prob = coco_add_glue(prob, 'glue', uidx(maps.T_idx), uidx(maps.p_idx(1)));
  prob = coco_set(prob, 'cont', 'NAdapt', 10);
  cont_args = {1, {'po.period' 'T'}, [2*pi/1.3 2*pi/0.7]};
  
  fprintf(...
  '\n Run=''%s'': Continue family of periodic orbits from point %d in run ''%s''.\n', ...
  sprintf('lab=%d', lab), lab, 'saddle-node');

  coco(prob, sprintf('lab=%d', lab), [], cont_args{:});
end


%% Free vibrations - backbone

% In this continuation problem, we seek periodic responses in the absence
% of excitation. These only exist for zero damping, and constitute a
% one-dimensional family of orbits parameterized by the orbit period. The
% encoding includes four monitor functions that evaluate to the problem
% parameters and the first component of the initial point on the periodic
% orbit, and the corresponding inactive continuation parameters 'T', 'A',
% 'd', and 'y0'. Holding 'y0' fixed corresponds to the imposition of a
% Poincare condition, which results in an overall dimensional deficit of
% -1. A one-dimensional family of periodic orbits is obtained by releasing
% 'po.period' and 'd', and allowing these to vary during continuation.

t0 = (0:0.01:2*pi)';
x0 = 2e-2*[sin(t0) cos(t0)];
p0 = [2*pi; 0; 0];
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'po', 'bifus', 'off');
funcs = {@bistable, @bistable_dx, @bistable_dp};
coll_args = [funcs, {t0, x0, {'T' 'A' 'd'}, p0}];
prob = ode_isol2po(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_pars(prob, 'section', uidx(maps.x0_idx(1)), 'y0');
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', [0 100]);
cont_args = {1, {'po.period'  'd'}, [2*pi/1.3 2*pi/0.7]};

fprintf('\n Run=''%s'': Continue backbone curve of periodic orbits.\n', ...
  'backbone');

coco(prob, 'backbone', [], cont_args{:});

%% Graphical representation of stored solutions

figure(1); clf; hold on; grid on; box on
coco_plot_bd('saddle-node', 'T', @(T) 2*pi./T, ...
    '||po.orb.x||_{L_2[0,T]}')
thm = struct('ustab', '', 'lspec', {{'r', 'LineWidth', 2}});
for lab=labs
  coco_plot_bd(thm, sprintf('lab=%d', lab), 'T', @(T) 2*pi./T, ...
    '||po.orb.x||_{L_2[0,T]}')
end
thm = struct('ustab', '', 'lspec', {{'b--', 'LineWidth', 2}}, 'xlab', '2\pi/T');
coco_plot_bd(thm, 'backbone', 'po.period', @(T) 2*pi./T, ...
    '||po.orb.x||_{L_2[0,T]}')
axis([0.7 1.3 0 inf]); hold off
