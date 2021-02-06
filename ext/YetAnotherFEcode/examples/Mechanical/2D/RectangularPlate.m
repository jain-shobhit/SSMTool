%% A temperature-dependent plate model

clear all; close all; clc
%% parameters
% geometry
l  = 20e-3; % length of domain [m]
b = 10e-3; % breadth of domain [m]
t = 8e-4; % thickness of plate [m]
f = 1; % RHS
g = 1; % Neumann BC RHS

% material properties
E       = 70e9;  % 70e9 % 200e9 % Young's modulus [Pa]
rho     = 2700; % 2700 % 7850 % density [kg/m^3]
nu      = 0.33;    % Poisson's ratio 
kappa   = 1e8; % material damping modulus 1e8
alpha_T = 11.7e-6; % thermal expansion coefficient

%% Structural
% Material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'THERMAL_EXPANSION_COEFFICIENT',alpha_T);
% Element
myElementConstructor =  @()TriShellElement(t, myMaterial); % same element all across the domain

% Mesh
nl = 20;
nb = 10; 
[Nodes,Elements,bnodes] = RectangularMesh(l,b,nl,nb,0);     % Rectangular Mesh definition

myMesh = Mesh(Nodes);
myMesh.create_elements_table(Elements,myElementConstructor);
% Plot mesh
figure(1); PlotMesh(Nodes,Elements,0);

% Assemble linear stiffness and mass
nNodes = size(Nodes,1);
PlateAssembly = Assembly(myMesh);
% PlateAssembly.parallelized = true;
K = PlateAssembly.stiffness_matrix();
M = PlateAssembly.mass_matrix();
C = PlateAssembly.damping_matrix();

% Tensor Assembly
T2 = PlateAssembly.tensor('T2',[myMesh.nDOFs, myMesh.nDOFs, myMesh.nDOFs], [2,3]);
T3 = PlateAssembly.tensor('T3',[myMesh.nDOFs, myMesh.nDOFs, myMesh.nDOFs, myMesh.nDOFs], [2,3,4]);

myMesh.set_essential_boundary_condition([bnodes{1}, bnodes{2}],1:3,0) % simply supported on two opposite edges
% myMesh.set_essential_boundary_condition([bnodes{1}],1:6,0) % cantilever on one edge

% Eigenvalue problem
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 
K_red = PlateAssembly.constrain_matrix(K);
M_red = PlateAssembly.constrain_matrix(M);
[V0,omega2] = eigs(K_red,M_red,n_VMs,'SM');
omega = sqrt(diag(omega2));

V0 = PlateAssembly.unconstrain_vector(V0);
mod = 1;
v1 = reshape(V0(:,mod),6,[]);
PlotFieldonDeformedMesh(Nodes,Elements,v1(1:3,:).','factor',5)
title(['Frequency = ' num2str(omega(mod)/(2*pi)) ' Hz'] )

% constrained tensors
T2 = PlateAssembly.constrain_tensor(T2);
T3 = PlateAssembly.constrain_tensor(T3);


%% Static response under uniform pressure
Pressure = 1e6; % in Pascals 
F = Pressure*PlateAssembly.uniform_body_force();
u_lin = PlateAssembly.solve_system(K,F);
ULIN = reshape(u_lin,6,[]);
figure(2); PlotMesh(Nodes,Elements,0);
hold on; PlotFieldonDeformedMesh(Nodes,Elements,ULIN(1:3,:).','factor',1)

% Nonlinear response
u = static_equilibrium( PlateAssembly, u_lin, F );
U = reshape(u,6,[]);
hold on
PlotFieldonDeformedMesh(Nodes,Elements,U(1:3,:).','factor',1, 'color', 'w' )
colormap gray

%% Dynamic response using Implicit Newmark
% forcing frequency of the average of first two natural frequencies
omega_ext = mean(omega(1:2)); 
T =  2*pi/omega_ext; % time period of forcing

% load amplification factor
amplification_factor = 0.1;

% forcing function
F_ext = @(t) amplification_factor * F * sin(omega_ext * t);

% Initial condition: equilibrium
u0 = zeros(PlateAssembly.Mesh.nDOFs, 1);
v0 = zeros(PlateAssembly.Mesh.nDOFs, 1);
a0 = zeros(PlateAssembly.Mesh.nDOFs, 1); % a0 = M\(F_ext(0)-C*v0-F(u0)) 

q0 = PlateAssembly.constrain_vector(u0);
qd0 = PlateAssembly.constrain_vector(v0);
qdd0 = PlateAssembly.constrain_vector(a0);

% time step for integration
h = T/50;

% Precompute data for Assembly object
PlateAssembly.DATA.M = M;
PlateAssembly.DATA.K = K;
PlateAssembly.DATA.C = C;

% Instantiate object for linear time integration
TI_lin = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,PlateAssembly,F_ext);

% Linearized Time Integration
tmax = 10*T; 
TI_lin.Integrate(q0,qd0,qdd0,tmax,residual_lin);

% obtain full solution
TI_lin.Solution.u = PlateAssembly.unconstrain_vector(TI_lin.Solution.q);

% Animate solution on Mesh (very slow)
% AnimateFieldonDeformedMesh(myMesh.Nodes,myMesh.Elements,TI_lin.Solution.u ,'factor',1,'index',1:3,'filename','lineardisp')

% Instantiate object for nonlinear time integration
TI_NL = ImplicitNewmark('timestep',h,'alpha',0.005);

% Linear Residual evaluation function handle
residual = @(q,qd,qdd,t)residual_nonlinear(q,qd,qdd,t,PlateAssembly,F_ext);

% Nonlinear Time Integration
tmax = 10*T; 
TI_NL.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL.Solution.u = PlateAssembly.unconstrain_vector(TI_NL.Solution.q);

%% Generalized alpha scheme
% linear
TI_lin_alpha = GeneralizedAlpha('timestep',h,'rho_inf',0.7, 'linear',true);
TI_lin_alpha.Integrate(q0,qd0,qdd0,tmax,residual_lin);
TI_lin_alpha.Solution.u = PlateAssembly.unconstrain_vector(TI_lin_alpha.Solution.q);

% nonlinear
TI_NL_alpha = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
TI_NL_alpha.Integrate(q0,qd0,qdd0,tmax,residual);
TI_NL_alpha.Solution.u = PlateAssembly.unconstrain_vector(TI_NL_alpha.Solution.q);


%% Reduced solution Linear
m = 5; % use the first five VMs in reduction
V = V0(:,1:m);

PlateReducedAssembly = ReducedAssembly(myMesh,V);

PlateReducedAssembly.DATA.M = PlateReducedAssembly.mass_matrix();
PlateReducedAssembly.DATA.C = PlateReducedAssembly.damping_matrix();
PlateReducedAssembly.DATA.K =  PlateReducedAssembly.stiffness_matrix();

q0 = zeros(m,1);
qd0 = zeros(m,1);
qdd0 = zeros(m,1);

TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Modal linear Residual evaluation function handle
Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear(q,qd,qdd,t,PlateReducedAssembly,F_ext);

% time integration
TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;

%% Reduced solution Noninear
% For demonstration purposes, we simply reduce the nonlinear system using
% out-of-plane bending modes. This is expected to produce bad results when 
% in-plane stretching is involved in the response.


TI_NL_alpha_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);

% Modal nonlinear Residual evaluation function handle
Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear(q,qd,qdd,t,PlateReducedAssembly,F_ext);

% time integration
TI_NL_alpha_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
TI_NL_alpha_red.Solution.u = V * TI_NL_alpha_red.Solution.q;


%% Comparison
% Linear
dof = 100; % random degree of freedom at which time response is compared
figure;
plot(TI_lin.Solution.time, TI_lin.Solution.u(dof,:),'DisplayName', 'Full linear (Newmark)')
hold on
plot(TI_lin_alpha.Solution.time, TI_lin_alpha.Solution.u(dof,:),'DisplayName', 'Full linear (Generalized-$$\alpha$$)')
plot(TI_lin_red.Solution.time, TI_lin_alpha.Solution.u(dof,:),'DisplayName', 'Reduced linear (Newmark)')
xlabel('q'); ylabel('time'); grid on; axis tight; legend('show')

% Nonlinear
figure;
plot(TI_NL.Solution.time, TI_NL.Solution.u(dof,:),'DisplayName', 'Full nonlinear (Newmark)')
hold on
plot(TI_NL_alpha.Solution.time, TI_NL_alpha.Solution.u(dof,:),'DisplayName', 'Full nonlinear (Generalized-$$\alpha$$)')
plot(TI_NL_alpha_red.Solution.time, TI_NL_alpha_red.Solution.u(dof,:),'DisplayName', 'Reduced nonlinear (Generalized-$$\alpha$$)')
xlabel('q'); ylabel('time'); grid on; axis tight; legend('show')

