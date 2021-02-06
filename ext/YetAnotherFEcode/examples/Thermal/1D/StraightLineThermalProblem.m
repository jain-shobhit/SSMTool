clear all; close all; clc
%% parameters
% geometry
l = 1;          % length of domain
f = 50;         % RHS
g = 1;          % Neumann BC RHS
r = 0;          % Robin factor dT + r*T = g 
TAway = 25;     % Temperature far away 
StBol = 5.670373e-8; 
% mesh
nElements = 10;
dx = l/nElements;

% material properties
rho     = 7850;              % 2700 % 7850 % density
k = 1;                       % conductivity
c = 1;                       % heat capacity
alpha = 10;                  % Convective constant
Emiss = 0;                 % Emissivety constant

%% Structural
% Material
myMaterial  = Material();
set(myMaterial,'THERMAL_CONDUCTIVITY', k, 'HEAT_CAPACITY', c, 'DENSITY', rho, 'CONVECTION_COEFFICIENT', alpha, 'EMISSIVITY', Emiss);

% Mesh
Nodes = (0:dx:l).';
nNodes = length(Nodes);
Elements = [1:nNodes-1;2:nNodes].';

% Element
myElementConstructor = @()TE1D(myMaterial); % same element all across the domain

myMesh = Mesh(Nodes);
myMesh.create_elements_table(Elements,myElementConstructor);

% Set dirichlet DOFs
myMesh.set_essential_boundary_condition(1,1,0); %[1; myMesh.nNodes]

% Set Neumann Elements
NeumannElements = nNodes; % 0 dimensional (point) element for Neumann BC
myNeumannElement = @()TE1D_Neumann(myMaterial,g,r);
myMesh.create_elements_table(NeumannElements,myNeumannElement,'isBoundary',true);

% Assembly
T0 = zeros(myMesh.nDOFs,1);
myAssembly = Assembly(myMesh);
[K, Fint] = myAssembly.tangent_stiffness_and_force(T0);

M = myAssembly.mass_matrix();

% force
% uniform transverse external force
F = f*myAssembly.uniform_body_force() + TAway*alpha*myAssembly.uniform_body_force()+ TAway^4*StBol*Emiss*myAssembly.uniform_body_force();

% Temperature result
mid = 2:myMesh.nNodes;
Boundary = setdiff(1:nNodes,mid);
T0(mid) = K(mid,mid)\(K(mid,Boundary)*T0(Boundary)+ F(mid)+Fint(mid));

plot(0:dx:l,T0)

