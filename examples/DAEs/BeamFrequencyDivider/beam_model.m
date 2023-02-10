function [M,C,K,fnl,fext, outdof, nodes, elements] = beam_model(l,b,h,BC,nElements,isViscoelastic,xOut)
%% Finite Element Setup
% Geometry
% Material properties
startLIN = tic;

pars = struct();
pars.nElements = nElements;
pars.isViscoelastic = isViscoelastic;
pars.l = l;
pars.h = h;
pars.b = b; 

% Material properties
pars.E       = 70e6;         % 70e9 Young's modulus
pars.rho     = 2700*10^(-9); % 2700 density
pars.nu      = 0.3;          % nu
pars.kappa   = 1e4;          % material damping modulus 1e7

%% FE model
disp('Building FE model')
% Material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',pars.E,'DENSITY',pars.rho,'POISSONS_RATIO',pars.nu,'DAMPING_MODULUS',pars.kappa);
% Element
if isViscoelastic
    myElementConstructor = @()BeamElementViscoelastic(pars.b, pars.h, myMaterial); % same element all across the domain
else
    myElementConstructor = @()BeamElement(pars.b, pars.h, myMaterial); % same element all across the domain
end

% Meshing the geometry
dx = pars.l/nElements;
x = (0:dx:pars.l).';
nNodes = size(x,1);
nodes = [x, zeros(nNodes,1)];
elements = [1:nNodes-1;2:nNodes].';

% creating Mesh object
MyMesh = Mesh(nodes);
MyMesh.create_elements_table(elements,myElementConstructor);

% Plot mesh
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
C = MyAssembly.damping_matrix();

%% apply boundary conditions
disp('Applying boundary conditions')
switch BC
    case 'CC' % clamped-clamped
        constrainedNodes = [1 nNodes];
        constrainedDOFs = 1:3;
        MyMesh.set_essential_boundary_condition(constrainedNodes,constrainedDOFs,0)
    case 'C'  % cantilevered
        constrainedNodes = 1;
        constrainedDOFs = 1:3;
        MyMesh.set_essential_boundary_condition(constrainedNodes,constrainedDOFs,0) 
    case 'PP' % pinned-pinned (axially restrained)
        constrainedNodes = [1 nNodes];
        constrainedDOFs = 1:2;
        MyMesh.set_essential_boundary_condition(constrainedNodes,constrainedDOFs,0) 
    case 'PPF' % pinned-pinned (axially unrestrained)
        constrainedNodes = 1;
        constrainedDOFs = [1 2];
        MyMesh.set_essential_boundary_condition(constrainedNodes,constrainedDOFs,0) 
        constrainedNodes = nNodes;
        constrainedDOFs = 2;
        MyMesh.set_essential_boundary_condition(constrainedNodes,constrainedDOFs,0) 
end
        
M = MyAssembly.constrain_matrix(M);
K = MyAssembly.constrain_matrix(K);
C = MyAssembly.constrain_matrix(C);

computationTimeLIN = toc(startLIN);

%% Tensor Assembly
%% Tensor Assembly
disp('Getting nonlinearity coefficients')
filename = ['tensors_' num2str(MyMesh.nElements) '.mat'];
computeTensors = false;
if isfile(filename)    
    load(filename,'parameters');
    if isequal(parameters,pars)
        load(filename,'fnl')        
        disp('Loaded tensors from storage')
        load(filename, 'computationTimeTensors')
    else
        computeTensors = true;        
    end
else
    computeTensors = true;
end

if computeTensors
    fnl = cell(1,2);
    disp('Assembling Tensors')
    parameters = pars;
    startTensors = tic;
    if isViscoelastic
        fnl{1} = MyAssembly.tensor('T2',[MyMesh.nDOFs, 2*MyMesh.nDOFs, 2*MyMesh.nDOFs], [2,3], MyMesh.nDOFs);
        fnl{2} = MyAssembly.tensor('T3',[MyMesh.nDOFs, 2*MyMesh.nDOFs, 2*MyMesh.nDOFs, 2*MyMesh.nDOFs], [2,3,4], MyMesh.nDOFs);       
    else
        fnl{1} = MyAssembly.tensor('T2',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3]);
        fnl{2} = MyAssembly.tensor('T3',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3,4]);
    end
    computationTimeTensors = toc(startTensors);
    disp('Saving Tensors')
    save(filename,'fnl','computationTimeTensors','parameters','-v7.3')
end

disp(['Total time spent on model assembly = ' datestr(datenum(0,0,0,0,0,computationTimeTensors + computationTimeLIN),'HH:MM:SS')])

% apply boundary conditions
for j = 1:length(fnl)
    fnl{j} = MyAssembly.constrain_tensor(fnl{j});
end


%% external force assembly
disp('Assembling external force vector')

[~,outnode] = min(abs(x- xOut));
outdof = outnode*3-2:outnode*3; % transverse direction

outdofvec = sparse(outdof,ones(size(outdof)),1,MyMesh.nDOFs,1);
outdofvec = MyAssembly.constrain_vector(outdofvec);
outdof = find(outdofvec);

% fext = outdofvec;
weights = true(nElements,1); 
fext = MyAssembly.constrain_vector(MyAssembly.uniform_body_force('weights',weights));
end





