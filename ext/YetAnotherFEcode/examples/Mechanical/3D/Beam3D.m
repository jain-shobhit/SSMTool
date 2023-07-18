% EXAMPLE: beam meshed with 3D elements
clear;
close all;
clc

% ***** Elements with linear shape functions (very slow mesh convergence, 
%     good for fast code-testing)
% elementType = 'HEX8';
% elementType = 'TET4';

% ***** Elements with quadratic shape functions
% elementType = 'TET10';
elementType = 'WED15';
% elementType = 'HEX20';

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
    case 'HEX8'
        myElementConstructor = @()Hex8Element(myMaterial);
    case 'HEX20'
        myElementConstructor = @()Hex20Element(myMaterial);
    case 'TET4'
        myElementConstructor = @()Tet4Element(myMaterial);
    case 'TET10'
        myElementConstructor = @()Tet10Element(myMaterial);
    case 'WED15'
        myElementConstructor = @()Wed15Element(myMaterial);
%         quadrature = struct('lin', 5, 'tri', 12);
%         myElementConstructor = @()Wed15Element(myMaterial, quadrature);
end

% MESH_____________________________________________________________________
l = 3;
w = .3;
t = .05;
nx = 15;
ny = 3;
nz = 2;
[nodes, elements, nset]=mesh_3Dparallelepiped(elementType,l,w,t,nx,ny,nz);
% % Alternatively, one can write an input file in ABAQUS and read it as:
% filename = 'Job-BeamHex';
% [nodes, elements, nset, elset] = mesh_ABAQUSread(filename); % HEX20 mesh

myMesh = Mesh(nodes);
myMesh.create_elements_table(elements,myElementConstructor);

% MESH > BOUNDARY CONDITONS
myMesh.set_essential_boundary_condition([nset{1} nset{4}],1:3,0)
% myMesh.BC.set_dirichlet_dofs([nset{2} nset{3}],1:3,0) % abaqus

% ASSEMBLY ________________________________________________________________
BeamAssembly = Assembly(myMesh);
M = BeamAssembly.mass_matrix();
nNodes = size(nodes,1);
u0 = zeros( myMesh.nDOFs, 1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);

% store matrices
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.M = M;


%% EXAMPLE 1: vibration modes                                       

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

% PLOT
mod = 1;
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMesh(nodes,elements,0);
v1 = reshape(V0(:,mod),3,[]).';
PlotFieldonDeformedMesh(nodes,elements,v1,'factor',1);
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])
drawnow


%% EXAMPLE 2: static test                                           

% % Define external force:
% % Body force
% Pressure = 1e5;
% F = Pressure * BeamAssembly.uniform_body_force();

% Nodal force
F = zeros(myMesh.nDOFs,1);
nf = find_node(l/2,w/2,t/2,nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
F(node_force_dofs(3)) = 1e3;

% u_lin = BeamAssembly.solve_system(K, F);
% ULIN = reshape(u_lin,3,[]).';	% Linear response
% u = static_equilibrium(BeamAssembly, F, 'display', 'iter-detailed');
% UNL = reshape(u,3,[]).';        % Nonlinear response
% 
% fprintf(['\n <strong>Max displacements</strong>:\n  Linear:\t\t%.3i \n' ...
%     '  Nonlinear:\t%.3i \n\n'],max(u_lin(:)),max(u(:)))
% 
% % PLOT
% figure('units','normalized','position',[.2 .1 .6 .8])
% scale = 10;
% PlotMesh(nodes,elements,0);
% PlotFieldonDeformedMesh(nodes,elements,UNL,'factor',scale,'color','w');
% colormap jet
% title(['NONLINEAR STATIC RESPONSE (scale factor: ' num2str(scale) 'x)'])


%% EXAMPLE 3: compute modal derivatives                             

% take the first eigenmode
Phi1 = V0(:,1);
Phi1_c = BeamAssembly.constrain_vector(Phi1);

% compute dK/dq1
dK = BeamAssembly.stiffness_derivative(Phi1);
dK_c = BeamAssembly.constrain_matrix(dK);

% compute MD11
MD11_c = -Kc\(dK_c * Phi1_c);
MD11 = BeamAssembly.unconstrain_vector(MD11_c);
MD11 = -MD11/max(abs(MD11(:))); % normalize to 1

% PLOT
figure('units','normalized','position',[.2 .1 .6 .8])
v1 = reshape(MD11,3,[]).';
PlotFieldonDeformedMesh(nodes,elements,v1,'factor',.1,'component','U1');
title('MD_{11} (U1)')
drawnow


%% EXAMPLE 4: reduced stiffness tensors                             

RB = [V0(:,1) MD11]; % reduced basis with VM_1 and MD_11
m = size(RB,2);

RBeamAssembly = ReducedAssembly(myMesh, RB);

% reduced stiffness tensors
K2r = RB'*K*RB;
ndofs_per_element = myMesh.Elements(1).Object.nNodes * myMesh.Elements(1).Object.nDim;
if m > ndofs_per_element 
    % compute element tensor (size ndofs_per_element), then project
    mode = 'standard';
    disp(' Standard tensor construction:')
    
    tic
    K3r = RBeamAssembly.tensor('T2',[m m m], [2 3], mode);
    fprintf(' K3r... %.2f s (%.1f elem/s)\n',toc,myMesh.nElements/toc)
    
    tic
    K4r = RBeamAssembly.tensor('T3',[m m m m], [2 3 4], mode);
    fprintf(' K4r... %.2f s (%.1f elem/s)\n',toc,myMesh.nElements/toc)
else
    % use Element-Level Projection (directly computes the projected tensor, of size m)
    mode = 'ELP';
    disp(' Element-wise projection:')
    
    tic
    K3r = RBeamAssembly.tensor('T2',[m m m], [2 3], mode);
    fprintf(' K3r... %.2f s (%.1f elem/s)\n',toc,myMesh.nElements/toc)

    tic
    K4r = RBeamAssembly.tensor('T3',[m m m m], [2 3 4], mode);
    fprintf(' K4r... %.2f s (%.1f elem/s)\n',toc,myMesh.nElements/toc)
end

% compute tensors for the tangent stiffness matrix
K3rt = K3r + permute(K3r, [1 3 2]); 
K4rt = K4r + permute(K4r, [1 3 2 4]) + permute(K4r, [1 4 2 3]);




