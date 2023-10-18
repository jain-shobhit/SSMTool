function [M,K] = build_linear_model(nDiscretization,w)
% we tune w to trigger 1:2 resonance

%% Finite Element Setup
% Geometry
l  = 1; % length of domain [m]
b = 2; % breadth of domain [m]
t = 1e-2; % thickness of plate [m]
% w = 1e-1; % curvature parameter (height of the midpoint relative to ends) [m]
% material properties
E       = 70e9;  % 70e9 % 200e9 % Young's modulus [Pa]
rho     = 2700; % 2700 % 7850 % density [kg/m^3]
nu      = 0.33;    % Poisson's ratio 
kappa   = 0; % material damping modulus 1e8

%% FE model
disp('Building FE model')
% Material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'DAMPING_MODULUS',kappa);
% Element
myElementConstructor =  @()TriShellElement(t, myMaterial); % same element all across the domain

% Meshing the geometry
nl = nDiscretization;
nb = 2*nDiscretization; 
[nodes,elements,bnodes] = RectangularMesh(l,b,nl,nb,w);     % Rectangular Mesh definition

% creating Mesh object
MyMesh = Mesh(nodes);
MyMesh.create_elements_table(elements,myElementConstructor);

% % Plot mesh
% figure('Name','Mesh'); PlotMesh(nodes,elements,0);

%% Assemble linear stiffness, mass and damping
disp('Assembling M,C,K matrices')
% % parallelized assembly
% cluster = parcluster('local');
% cluster.NumWorkers = 4;
% parpool(cluster, 4)
% MyAssembly = Assembly(myMesh,true); 

MyAssembly = Assembly(MyMesh);
K = MyAssembly.stiffness_matrix();
M = MyAssembly.mass_matrix();
% C = MyAssembly.damping_matrix();

%% apply boundary conditions
disp('Applying boundary conditions')
MyMesh.set_essential_boundary_condition([bnodes{3}, bnodes{4}],1:3,0) % simply supported on opposite ends
M = MyAssembly.constrain_matrix(M);
K = MyAssembly.constrain_matrix(K);
% C = MyAssembly.constrain_matrix(C);

end

