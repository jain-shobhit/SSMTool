%% construct initial periodic solution using forward simulation
p0 = [3.5; 1.5; 0]; % [om ro eps]
T  = 2*pi/p0(1);
[~,x0] = ode45(@(t,x) lang(x,p0), 0:100*T, [0.3; 0.4; 0]); % transient
[t0,x0] = ode45(@(t,x) lang(x,p0), linspace(0,T,100), x0(end,:)); % approximated periodic solution
figure;
plot(t0,x0);

%% continuation of periodic orbit
prob = coco_prob();
prob = ode_isol2po(prob, '', @lang, @lang_DFDX, @lang_DFDP, ...
  t0, x0, {'om','rho','eps'}, p0);

fprintf('\n Run=''%s'': Continue family of periodic orbits.\n', ...
  'po');

coco(prob, 'po', [], 1, 'rho', [0.2 2]);

%% continuation of torus from TR point with varied varrho
T_po   = 5.3; 
T_ret  = 2*pi/3.5;
varrho = T_ret/T_po;
bd    = coco_bd_read('po');
TRlab = coco_bd_labs(bd, 'TR');
prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'h_min',...
    1e-3, 'PtMX', 50, 'h_max', 10, 'bi_direct', false);
prob = ode_TR2tor(prob, '', 'po', TRlab, 50);
fprintf(...
  '\n Run=''%s'': Continue torus from point %d in run ''%s''.\n', ...
  'tr1', TRlab, 'po');
coco(prob, 'tr1', [], 1, {'varrho','rho','om1','om2','eps'},[varrho,0.44]);

%% continuation of torus from previous solution
bd   = coco_bd_read('tr1');
lab  = coco_bd_labs(bd, 'EP');
lab  = max(lab);
prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'h_min',...
    1e-3, 'PtMX', 30, 'h_max', 10, 'bi_direct', true);
prob = ode_tor2tor(prob, '', 'tr1', lab);
fprintf(...
  '\n Run=''%s'': Continue torus from point %d in run ''%s''.\n', ...
  'tr2', lab, 'tr1');
coco(prob, 'tr2', [], 1, {'eps','rho','om1','om2','varrho'});


%% visualization
figure; coco_plot_bd('tr2','rho','eps');
for lab = 1:5
    plot_torus('','tr1', lab, [1 2 3], 0.75); pause(1);
end

% time integration
lab = 5;
plot_torus('','tr1', lab, [1 2 3]);hold on
sol = tor_read_solution('','tr1',lab);
p   = sol.p(1:end-3);
xbp = sol.xbp;
[t,x] = ode45(@(t,x) lang(x,p), 0:0.01:100*T, xbp(1,:,1));
plot3(x(:,1),x(:,2),x(:,3),'r-');

