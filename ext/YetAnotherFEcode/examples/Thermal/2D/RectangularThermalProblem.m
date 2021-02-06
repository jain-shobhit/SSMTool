
clear all; close all; clc
%% parameters
% geometry
l = 1;          % length of domain
h = 1;          % length of height
f = 100;        % RHS
g = 100;        % Neumann BC RHS
r = 1;          % Activation of Robin cond.-> '1' if it takes place,'0' otherwise.
TAway = 25;     % Temperature far away 
StBol = 5.670373e-8; 

% mesh
nElements = 20;
dx = l/nElements;
dy = h/nElements;

% material properties
E = 200e9;              % 70e9 % 200e9 % Young's modulus
rho = 7850;             % 2700 % 7850 % density
k = 1;                  % conductivity
alpha = 10;              % heat convective coefficient
c = 1;                  % heat capacity
nu = 0.3;               % nu
kappa = 1e8;            % material damping modulus 1e8
Emiss = 0;              % Emissivety constant
alpha_T = 11.7e-6;      % thermal expansion coefficient
loadfactor = 50;        % mechanical load factor 

%% Structural
% Material
myMaterial  = Material();
set(myMaterial,'THERMAL_CONDUCTIVITY', k, 'HEAT_CAPACITY', c, 'DENSITY', rho, 'CONVECTION_COEFFICIENT', alpha, 'EMISSIVITY', Emiss);


% Mesh
[Nodes,Elements,bnodes] = RectangularMesh(l,h,nElements,nElements,0);     % Rectangular Mesh definition
nNodes = length(Nodes); 
NodesMesh = reshape(1:nNodes,nElements+1,nElements+1)';


% Element
myThermalElement = @()TE2D(myMaterial); % same element all across the domain
myMesh = Mesh(Nodes);
myMesh.create_elements_table(Elements,myThermalElement);

% Set dirichlet DOFs
InternalIndex = NodesMesh(1:end-1,2:end-1);
Boundary = setdiff(1:nNodes,InternalIndex(:));
myMesh.set_essential_boundary_condition(Boundary,1,0) %

% Set Neumann Elements
NeumannElements =  [NodesMesh(1,1:end-1).' NodesMesh(1,2:end).']; %TableElementNeumann(NodesMesh(1,:)); 
myNeumannElement = @()TE2D_Neumann(myMaterial,g,r);
myMesh.create_elements_table(NeumannElements,myNeumannElement,'isBoundary',true);

% Assembly
T0 = zeros(myMesh.nDOFs,1);
myAssembly = Assembly(myMesh);
[K, Fint] = myAssembly.tangent_stiffness_and_force(T0);

M = myAssembly.mass_matrix();

% force
% uniform transverse external force
F = f*myAssembly.uniform_body_force()+ TAway*alpha*myAssembly.uniform_body_force()+ TAway^4*StBol*Emiss*myAssembly.uniform_body_force();

% Temperature result
mid = InternalIndex(:);
T0(mid) = K(mid,mid)\(K(mid,Boundary)*T0(Boundary) + F(mid)-Fint(mid));    % This is without radiation  
% TRad = NonlinearSol(myAssembly,T0,F,Emiss);
% T0(mid) = TRad(mid);
T0Mesh = reshape(T0,nElements+1,nElements+1)'; 
T0MeshPlot = T0Mesh(end:-1:1,:);

X = 0:dx:l;
Y = 0:dy:h;
figure
surf(X,Y,T0MeshPlot)

