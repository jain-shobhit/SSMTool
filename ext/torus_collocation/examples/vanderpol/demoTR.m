%% Continuation of quasiperiodic invariant tori
%
% Two-dimensional quasiperiodic invariant tori are characterized by a
% parallel flow with irrational rotation number. In terms of a suitably
% defined torus function, the dynamics on the torus may be described by an
% invariant circle that is mapped to itself by the flow, such that the flow
% is equivalent to a rigid rotation on the circle.
%
% We obtain an approximate description of a quasiperiodic invariant torus
% in terms of a Fourier representation of the invariant circle and a finite
% collection of trajectory segments based at points on the invariant
% circle. The rigid rotation imposes an all-to-all coupled system of
% boundary conditions on the collection of trajectory segments.

%% Encoding

% Construct the guess of initial solution

% pnames = [  om     c     a ]
p0 = [ 1.0; 0.11; 0.2 ];
T  = 2*pi/p0(1);
[~,x0] = ode45(@(t,x) vdp(t,x,p0), 0:100*T, [0; 0]); % Approximate periodic orbit of reduced system
[t0,x0] = ode45(@(t,x) vdp(t,x,p0), linspace(0,T,100), x0(end,:));
figure;
plot(t0,x0);


%% continuation of periodic orbit
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'h_max', 5);
prob = coco_set(prob, 'coll', 'NTST', 20);
prob = ode_isol2po(prob, '', @vdp, @vdp_DFDX, @vdp_DFDP, @vdp_DFDT,... 
  t0, x0, {'Om2','c','a'}, p0);
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, [], 'zero',...
    'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);

fprintf('\n Run=''%s'': Continue family of periodic orbits.\n', ...
  'po');

coco(prob, 'po', [], 1, {'Om2','po.period'}, [0.7 1.3]);

figure; coco_plot_bd('po')

%% continuation of TR points
bd     = coco_bd_read('po');
TRlabs = coco_bd_labs(bd, 'TR');
TRlab  = max(TRlabs);
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'h_max', 5, 'bi_direct', false);
prob = ode_po2TR(prob, '', 'po', TRlab);
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, [], 'zero',...
    'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);

fprintf('\n Run=''%s'': Continue family of TR points.\n', ...
  'po_TR');

coco(prob, 'po_TR', [], 1, {'Om2','po.period','a'}, [1.0 1.5]);

%% continuation of torus from TR point (from po2TR)
bd  = coco_bd_read('po_TR');
lab = coco_bd_labs(bd, 'EP');
lab = max(lab);
prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'PtMX', 100,...
    'h_max', 15, 'bi_direct', false);
prob = coco_set(prob, 'corr', 'ItMX', 20);
prob = coco_set(prob, 'coll', 'NTST', 20);
prob = ode_TR2tor(prob, '', 'po_TR', lab, 'neg', 1e-5);

fprintf(...
  '\n Run=''%s'': Continue torus from point %d in run ''%s''.\n', ...
  'tr1', lab, 'po_TR');
coco(prob, 'tr1', [], 1, {'varrho', 'Om2','om2','om1','c','a'},[-0.5,-0.3]);

for lab = 1:7
    plot_torus('','tr1',lab, [1 2]); pause(1);
end


% time integration
lab = 6;
plot_torus('','tr1', lab, [1 2]);hold on
sol = tor_read_solution('','tr1',lab);
p   = sol.p(1:end-3);
xbp = sol.xbp;
T   = 2*pi/sol.p(end-1);
[t,x] = ode45(@(t,x) vdp(t,x,p), 0:0.01:50*T, xbp(1,:,1));
plot3(x(:,1),mod(t,T),x(:,2),'r.');
