function [M,C,K,fnl,fext, outdof] = build_model(nElements)
%% Finite Element Setup
assert(mod(nElements,2)==0, 'The number of elements is not an even number.');
% Geometry

h=10;                                   %Height of beam [mm]
b=10;                                   %Width of beam [mm]
l=2700;                                 %Length of beam [mm]
% Mesh parameters

% Material properties
E=45000000;                             %Young's modulus [kPa]   
rho=1780*10^(-9);                       %Density [kg/mm^3]    
kappa = 1e3;

% E       = 70e9;  % 70e9 % 200e9 % Young's modulus
% rho     = 2700; % 2700 % 7850 % density
nu      = 0.3;    % nu
% kappa   = 1e7; % material damping modulus 1e8

%% FE model
disp('Building FE model')
% Material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'DAMPING_MODULUS',kappa);
% Element
myElementConstructor = @()BeamElement(b, h, myMaterial); % same element all across the domain

% Meshing the geometry
dx = l/nElements;
x = (0:dx:l).';
nNodes = size(x,1);
nodes = [x, zeros(nNodes,1)];
elements = [1:nNodes-1;2:nNodes].';

% creating Mesh object
MyMesh = Mesh(nodes);
MyMesh.create_elements_table(elements,myElementConstructor);

% Plot mesh
figure('Name','Mesh'); PlotMesh(nodes,elements,0);

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
C = MyAssembly.damping_matrix();

%% apply boundary conditions
disp('Applying boundary conditions')
MyMesh.set_essential_boundary_condition(1,[1 2 3],0)      % clamped at left end
MyMesh.set_essential_boundary_condition(nNodes, [1 2],0); % hinged at right end
M = MyAssembly.constrain_matrix(M);
K = MyAssembly.constrain_matrix(K);
C = MyAssembly.constrain_matrix(C);



%% Eigenvalue problem
disp('Solving undamped eigenvalue problem')
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 
[V0,omega2] = eigs(K,M,n_VMs,'SM');
omega = sqrt(diag(omega2));

V = MyAssembly.unconstrain_vector(V0);
mode = 1;
v1 = reshape(V(:,mode),3,[]);
PlotFieldonDeformedMesh(nodes,elements,v1(1:2,:).','factor',100);
title(['Mode ' num2str(mode) ', ' 'Elements ' num2str(nElements) ', ' 'Frequency = ' num2str(omega(mode)/(2*pi)) ' Hz'] )

%% Tensor Assembly
disp('Getting nonlinearity coefficients')

fnl = cell(1,2);
disp('Assembling Tensors')
fnl{1} = MyAssembly.tensor('T2',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3]);
fnl{2} = MyAssembly.tensor('T3',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3,4]);

% apply boundary conditions
for j = 1:length(fnl)
    fnl{j} = MyAssembly.constrain_tensor(fnl{j});
end

%% external force assembly
disp('Assembling external force vector')

outnode = [1 2]/4*nElements+1;
outdof = outnode*3-1; % transverse direction

outdofvec = sparse(outdof,ones(size(outdof)),1,MyMesh.nDOFs,1);
outdofvec = MyAssembly.constrain_vector(outdofvec);
outdof = find(outdofvec);

% fext = outdofvec;
weights = true(nElements,1); 
fext = MyAssembly.constrain_vector(MyAssembly.uniform_body_force('weights',weights));