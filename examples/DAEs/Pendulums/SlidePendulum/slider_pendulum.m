function slider_pendulum()

options = odeset('RelTol',1e-6,'AbsTol',1e-10);
epf = 1.06; omega = 1.793; % slightly larger than omega_a=1.7927
T = 2*pi/omega;
psp = [epf; omega];
init_state = [0;-pi/2;0;0];
[~, x0] = ode45(@(t,x) ode_func(t,x,psp), [0 100*T], init_state,options); % Transients from initial state
[ttor, xtor] = ode45(@(t,x) ode_func(t,x,psp), [0 T], x0(end,:)',options);   % Approximate periodic orbit

figure;
plot(xtor(:,1),xtor(:,3),'r.-');
xlabel('$x_1$','Interpreter','latex','FontSize',14);
ylabel('$\dot{x}_1$','Interpreter','latex','FontSize',14);
title(['$\omega=',num2str(omega),'$'],'Interpreter','latex','FontSize',14);

%% continuation run using po-toolbox of coco
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'cont', 'NAdapt', 2, 'h_max', 100,'PtMX',300);
prob = coco_set(prob, 'coll', 'NTST', 50);
funcs = {@ode_func};
coll_args = [funcs, {ttor, xtor, {'eps' 'Omega'}, psp}];
prob = ode_isol2po(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, [], 'zero',...
    'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(2))]);
ampdata.dof  = 1;
ampdata.zdim = 4;
prob = coco_add_func(prob, 'amp1', @amplitude, ampdata, 'regular', 'x1',...
    'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);

cont_args = {1, {'Omega' 'po.period' 'x1' 'eps'}, [1,7]};

fprintf('\n Run=''%s'': Continue primary family of periodic orbits.\n', ...
  'freq_resp');

bd1  = coco(prob, 'freq_resp', [], cont_args{:});

% Plot results
figure; hold on
coco_plot_bd('freq_resp','Omega','x1');
xlabel('$\omega$','Interpreter','latex','FontSize',14);
ylabel('$||{x}_1||_\infty$','Interpreter','latex','FontSize',14);
grid on

end


function y = ode_func(t,x,p)
% x: (x1,phi,dot(x1),dot(phi))
% p: (f,omega)
m1 = 1;
m2 = 1;
len = 1;
g = 9.8;
c1 = 0.17;
c2 = 0.02;
k1 = 7.45;
k2 = 1;


x1   = x(1,:);
phi  = x(2,:);
v1   = x(3,:);
vphi = x(4,:);

f  = p(1,:);
om = p(2,:);

y(1,:) = v1;
y(2,:) = vphi;
% we solve \ddot{x}_1 and \ddot{\phi} explicitly
detm = (m1+m2)*m2*len^2/3-0.25*m2^2*len^2*(sin(phi)).^2;
fx1  = f.*cos(om.*t)-k1*x1-c1*v1+0.5*m2*len*cos(phi).*vphi.^2;
fphi = -k2*(pi/2+phi)-c2*vphi-0.5*m2*g*len*cos(phi);
y(3,:) = (m2*len^2/3*fx1+0.5*m2*len*sin(phi).*fphi)./detm;
y(4,:) = (0.5*m2*len*sin(phi).*fx1+(m1+m2)*fphi)./detm;

end

