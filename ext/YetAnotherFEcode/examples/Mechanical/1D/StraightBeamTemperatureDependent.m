%% A temperature-dependent beam model

clear all; close all; clc
%% parameters
% geometry
l = 0.2; %1 % 0.2
h = 1e-3;
b = 1e-2; %1e-1 % 1e-2

w   = 5; %1/b*l is the width of the sin^2 graph of temperature pulse on the beam
c = 0.215; % location of the pulse center: factor of beam length
n_VMs = 5; % number of vibration modes
n_BMs = 5; % number of buckling modes
x_c = 0.08:0.04:0.92; % Center of the sin^2: factor of beam length

% mesh
nElements = 100;
dx = l/nElements;
outdof = 3*floor(nElements/4) + 2; % output degree of freedom for displaying results

% material properties
E       = 200e9;  % 70e9 % 200e9 % Young's modulus
rho     = 7850; % 2700 % 7850 % density
nu      = 0.3;    % nu
kappa   = 1e8; % material damping modulus 1e8
alpha_T = 11.7e-6; % thermal expansion coefficient
loadfactor = 1; % mechanical load factor 
%% Structural
% Material
myBeamMaterial  = KirchoffMaterial();
set(myBeamMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'DAMPING_MODULUS',kappa);
D               = myBeamMaterial.get_stress_strain_matrix_2D();

% Element
myElementConstructor = @()BeamElement(b, h, myBeamMaterial); % same element all across the domain

% Mesh
x = (0:dx:l).';
nNodes = size(x,1);
Nodes = [x, zeros(nNodes,1)];

Elements = [1:nNodes-1;2:nNodes].';
BeamMesh = Mesh(Nodes);
BeamMesh.create_elements_table(Elements,myElementConstructor);

% Assembly
u0 = zeros(nNodes*BeamMesh.nDOFPerNode,1);
BeamAssembly = Assembly(BeamMesh);
[K, F] = BeamAssembly.tangent_stiffness_and_force(u0);
M = BeamAssembly.mass_matrix();

u0 = randi(5,303,1);
F2 = BeamAssembly.vector('F2',u0);
T2 = BeamAssembly.tensor('T2',[BeamMesh.nDOFs, BeamMesh.nDOFs, BeamMesh.nDOFs], [2,3]);
F2check = ttv(T2,{u0,u0},[2,3]);
norm(F2check.data - F2)/norm(F2)

F3 = BeamAssembly.vector('F3',u0);
T3 = BeamAssembly.tensor('T3',[BeamMesh.nDOFs, BeamMesh.nDOFs, BeamMesh.nDOFs, BeamMesh.nDOFs], [2,3,4]);
F3check = ttv(T3,{u0,u0,u0},[2,3,4]);
norm(F3check.data - F3)/norm(F3)

S = BeamAssembly.scalar('strain_energy',u0);

%% Thermal beam
% Material
myThermalBeamMaterial = KirchoffMaterial();
set(myThermalBeamMaterial,'YOUNGS_MODULUS', E,'DENSITY',rho,'POISSONS_RATIO',nu,'THERMAL_EXPANSION_COEFFICIENT', alpha_T);

% Element
myElementConstructor = @()ThermalBeamElement(b, h, myThermalBeamMaterial); % same element all across the domain

% Mesh
ThermalBeamMesh = Mesh(Nodes);
ThermalBeamMesh.nDOFPerNode = nDOFperNode;
ThermalBeamMesh.create_elements_table(Elements,myElementConstructor);

% Set dirichlet DOFs
ThermalBeamMesh.set_essential_boundary_condition([1 ThermalBeamMesh.nNodes],[1 2 3],0) % Doubly clamped beam

ThermalBeamAssembly = Assembly(ThermalBeamMesh);
% force
% uniform transverse external force
f = ThermalBeamAssembly.uniform_body_force();






