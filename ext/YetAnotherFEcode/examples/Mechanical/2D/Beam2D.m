% EXAMPLE: beam meshed with 3D element
clear; 
close all; 
clc

whichModel = 'CUSTOM'; % or "ABAQUS"
elementType = 'QUAD4';
% elementType = 'QUAD8';

%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = .1;     % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
% Element
switch elementType
    case 'QUAD4'
        myElementConstructor = @()Quad4Element(thickness, myMaterial);
    case 'QUAD8'
        myElementConstructor = @()Quad8Element(thickness, myMaterial);
end

% MESH_____________________________________________________________________
Lx = 3;
Ly = .2;
nx = 30;
ny = 3;
switch upper( whichModel )
    case 'CUSTOM'
        [nodes, elements, nset] = mesh_2Drectangle(Lx,Ly,nx,ny,elementType);
    case 'ABAQUS'
        % Alternatively, one can write an input file in ABAQUS and read it as:
        filename = 'Job-BeamQuad';
        [nodes, elements, nset, elset] = mesh_ABAQUSread(filename);
end

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

% MESH > BOUNDARY CONDITONS
switch upper( whichModel )
    case 'CUSTOM'
        myMesh.set_essential_boundary_condition([nset{1} nset{3}],1:2,0)
    case 'ABAQUS'
        myMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)
end

% ASSEMBLY ________________________________________________________________
BeamAssembly = Assembly(myMesh);
M = BeamAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( myMesh.nDOFs, 1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);

% store matrices
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.M = M;

%% Tensor

nDOFs = BeamAssembly.Mesh.nDOFs;
u0 = randi(5,nDOFs,1);
F2 = BeamAssembly.vector('F2',u0,u0);
T2 = BeamAssembly.tensor('T2',[nDOFs, nDOFs, nDOFs], [2,3]);
F2check = ttv(T2,{u0,u0},[2,3]);
norm(F2check.data - F2)/norm(F2)

F3 = BeamAssembly.vector('F3',u0,u0,u0);
T3 = BeamAssembly.tensor('T3',[nDOFs, nDOFs, nDOFs, nDOFs], [2,3,4]);
F3check = ttv(T3,{u0,u0,u0},[2,3,4]);
norm(F3check.data - F3)/norm(F3)

[~,F] = BeamAssembly.tangent_stiffness_and_force(u0);
norm(K*u0 + F2 + F3 - F)/norm(F)

%% EXAMPLE 1                                                        

% Eigenvalue problem_______________________________________________________
n_VMs = 2; % first n_VMs modes with lowest frequency calculated 
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
[V0,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
V0 = V0(:,ind);
for ii = 1:n_VMs
    V0(:,ii) = V0(:,ii)/max(sqrt(sum(V0(:,ii).^2,2)));
end
V0 = BeamAssembly.unconstrain_vector(V0);

% PLOT
mod = 2;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes, elementPlot, 0);
v1 = reshape(V0(:,mod), 2, []).';
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', max(nodes(:,2)));
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])


%% EXAMPLE 2                                                        

% Define external force:
% % Body force
% Pressure = 1e6;
% F = Pressure*BeamAssembly.uniform_body_force();

% Nodal force
F = zeros(myMesh.nDOFs,1);
nf = find_node(Lx/2,Ly/2,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F(node_force_dofs(2)) = 10e5;

u_lin = BeamAssembly.solve_system(K, F);
ULIN = reshape(u_lin,2,[]).';	% Linear response
u = static_equilibrium(BeamAssembly, F, 'display', 'iter-detailed');
UNL = reshape(u,2,[]).';        % Nonlinear response

fprintf(['\n <strong>Max displacements</strong>:\n  Linear:\t\t%.3i \n' ... 
    '  Nonlinear:\t%.3i \n\n'],max(u_lin(:)),max(u(:)))

% PLOT
figure('units','normalized','position',[.2 .1 .6 .8])
scale = 2*Ly/max(abs(UNL(:)));
PlotMesh(nodes, elementPlot, 0);
PlotFieldonDeformedMesh(nodes,elementPlot,UNL,'factor',scale,'color','k');
colormap jet
title(['NONLINEAR STATIC RESPONSE (scale factor: ' num2str(scale) 'x)'])

