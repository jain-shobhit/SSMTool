function [M,C,K,fnl,fext, outdof] = build_model()
%% Read Mesh
startLIN = tic;
mshfile = 'Wing.msh';
disp(['Reading mesh from ' mshfile]) 
[nodes, elements, constrainedNodes, nbcElements] = read_msh(mshfile);

%% parameters
% geometry and material properties (mm units)
% E = 70e3; nu = 0.33; t = 0.8; rho = 2700e-12; kappa = 1;

E       = 70e9;  % 70e9 % 200e9 % Young's modulus [Pa]
rho     = 2700; % 2700 % 7850 % density [kg/m^3]
nu      = 0.33;    % Poisson's ratio 
kappa   = 1e5; % material damping modulus 1e8
t = 1.5e-3; % thickness of shell element
%% FE model
disp('Building FE model')
% Material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'DAMPING_MODULUS',kappa);
% Element
myElementConstructor =  @()TriShellElement(t, myMaterial); % same element all across the domain
% Mesh object
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
% C = MyAssembly.damping_matrix();
%% apply boundary conditions
disp('Applying boundary conditions')
% myMesh.set_essential_boundary_condition([bnodes{1}, bnodes{2}],1:3,0) % simply supported on two opposite edges
MyMesh.set_essential_boundary_condition(constrainedNodes,1:6,0) % cantilever on one edge

M = MyAssembly.constrain_matrix(M);
K = MyAssembly.constrain_matrix(K);
% C = MyAssembly.constrain_matrix(C);

%% Eigenvalue problem
disp('Solving undamped eigenvalue problem')
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 
[V0,omega2] = eigs(K,M,n_VMs,'SM');
omega = sqrt(diag(omega2));

V = MyAssembly.unconstrain_vector(V0);
mod = 1;
v1 = reshape(V(:,mod),6,[]);
PlotFieldonDeformedMesh(nodes,elements,v1(1:3,:).','factor',50);
title(['Mode ' num2str(mod) ', Frequency = ' num2str(omega(mod)/(2*pi)) ' Hz'] )

%% Damping matrix
disp('Using Rayleigh damping')
W =   omega(1:2);
a = [W(1) 1/W(1);W(2) 1/W(2)]\[0.004;0.004];
C = a(2) * M + a(1) * K;

%% external force assembly
disp('Assembling external force vector')

% nForceElements = numel(nbcElements);
% weights = sparse(nbcElements,ones(nForceElements,1),true,MyMesh.nElements,1);
% fext = MyAssembly.uniform_body_force(weights);
% fext = MyAssembly.constrain_vector(fext);


outnode = [337,349];
outdof = outnode*6-4; % y direction

% forcing tip nodes
fext = sparse(outdof,ones(size(outdof)),1,MyMesh.nDOFs,1);
fext = MyAssembly.constrain_vector(fext);
outdof = find(fext);

computationTimeLIN = toc(startLIN);
%% Tensor Assembly
disp('Getting nonlinearity coefficients')
filename = ['tensors_' num2str(MyMesh.nElements) '.mat'];
try     
    load(filename,'fnl')
    disp('Loaded tensors from storage')
    load(filename, 'computationTimeTensors')
    disp(['Total time spent on model assembly = ' datestr(datenum(0,0,0,0,0,computationTimeTensors + computationTimeLIN),'HH:MM:SS')])
catch
    fnl = cell(1,2);
    disp('Assembling Tensors')
    startTensors = tic;
    fnl{1} = MyAssembly.tensor('T2',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3]);
    fnl{2} = MyAssembly.tensor('T3',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3,4]);       
    computationTimeTensors = toc(startTensors);
    disp('Saving Tensors')
    save(filename,'fnl','computationTimeTensors','-v7.3')
    disp(['Total time spent on model assembly = ' datestr(datenum(0,0,0,0,0,computationTimeTensors + computationTimeLIN),'HH:MM:SS')])
end

% apply boundary conditions
for j = 1:length(fnl)
    fnl{j} = MyAssembly.constrain_tensor(fnl{j});
end
