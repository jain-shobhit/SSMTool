% EXAMPLE: beam meshed with 3D element
clear;
close all;
clc

% elementType = 'HEX20';
elementType = 'TET10';


%% PREPARE MODEL

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
% Element
switch elementType
    case 'HEX20'
        myElementConstructor = @()Hex20Element(myMaterial);
    case 'TET10'
        myElementConstructor = @()Tet10Element(myMaterial);
end

% MESH_____________________________________________________________________
l = 3;
w = .3;
t = .1;
nx = 10;
ny = 3;
nz = 2;
[nodes, elements, nset]=mesh_3Dparallelepiped(elementType,l,w,t,nx,ny,nz);
% % Alternatively, one can write an input file in ABAQUS and read it as:
% filename = 'Job-BeamHex';
% [nodes, elements, nset, elset] = mesh_ABAQUSread(filename); % HEX20 mesh

MyMesh = Mesh(nodes);
MyMesh.create_elements_table(elements,myElementConstructor);

% MESH > BOUNDARY CONDITONS
MyMesh.set_essential_boundary_condition([nset{1} nset{4}],1:3,0)
% MyMesh.BC.set_dirichlet_dofs([nset{2} nset{3}],1:3,0) % abaqus

% ASSEMBLY ________________________________________________________________
BeamAssembly = Assembly(MyMesh);
M = BeamAssembly.mass_matrix();
u0 = zeros(MyMesh.nDOFs,1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);


%% EXAMPLE 1

% Eigenvalue problem_______________________________________________________
n_VMs = 3; % first n_VMs modes with lowest frequency calculated
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
[V0,om] = eigs(Kc,Mc,n_VMs,'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
V0 = V0(:,ind);
for ii = 1:n_VMs
    V0(:,ii) = V0(:,ii)/max(sqrt(sum(V0(:,ii).^2,2)));
end
V0 = BeamAssembly.unconstrain_vector(V0);
mod = 1;

% PLOT
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes,elements,0);
v1 = reshape(V0(:,mod),3,[]).';
PlotFieldonDeformedMesh(nodes,elements,v1,'factor',1)
title(['$$\Phi_' num2str(mod) '$$ - Frequency = ' num2str(f0(mod),3) ' Hz'])

%% nonlinear tensors
% T2 = BeamAssembly.tensor('T2',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3]);
% T3 = BeamAssembly.tensor('T3',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3,4]);

%% EXAMPLE 2

% Define external force:
% % Body force
% Pressure = 1e6;
% F = Pressure*BeamAssembly.uniform_body_force();

% Nodal force
F = zeros(MyMesh.nDOFs,1);
nf = find_node(l/2,w/2,t/2,nodes); % node where to put the force

node_force_dofs = get_index(nf,MyMesh.nDOFPerNode);
F(node_force_dofs(3)) = 10e3;

u_lin = BeamAssembly.solve_system(K,F);
ULIN = reshape(u_lin,3,[]).';	% Linear response
u = static_equilibrium(BeamAssembly, u_lin/2, F);
UNL = reshape(u,3,[]).';        % Nonlinear response

fprintf(['\n <strong>Max displacements</strong>:\n  Linear:\t\t%.3i \n' ...
    '  Nonlinear:\t%.3i \n\n'],max(u_lin(:)),max(u(:)))

% PLOT
figure('units','normalized','position',[.2 .1 .6 .8])
scale = 100;
PlotMesh(nodes,elements,0);
PlotFieldonDeformedMesh(nodes,elements,UNL,'factor',scale,'color','w')
colormap jet
title(['NONLINEAR STATIC RESPONSE (scale factor: ' num2str(scale) 'x)'])

