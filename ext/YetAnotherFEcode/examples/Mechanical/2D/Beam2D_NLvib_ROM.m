% EXAMPLE: beam meshed with 2D element
%
% Contents (in brief):
% - NLvib
% - ROM
% - tensors
%
% Compute the nonlinear Frequency Response of a beam using NLvib AND a
% Reduced Order Model (ROM). The ROM is constructed projecting the
% equations of motion with a basis V = [VM MD] (Vibration Modes + Modal
% Derivatives). The internal forces are evaluated using (reduced order)
% stiffness tensors.
% --> usage example of the "custom" mode for the FE_system using NLvib

clear; 
close all; 
clc
format short g


%% PREPARE MODEL                                                    

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2800;     % density [kg/m^3]
nu      = 0.30;     % Poisson's ratio 
thickness = .2;     % [m] beam's out-of-plane thickness

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = false;	% set "true" for plane_stress
% Element
myElementConstructor = @()Quad8Element(thickness, myMaterial);

% MESH_____________________________________________________________________
Lx = 2;
Ly = .050;
h = 2;
nx = 20*h;
ny = 1*h;
[nodes, elements, nset] = mesh_2Drectangle(Lx, Ly, nx, ny, 'QUAD8');

% nominal mesh
MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);
MeshNominal.set_essential_boundary_condition([nset{1} nset{3}],1:2,0)

% new mesh
    % arch shape
    xi = 1;                                     % arch amplitude
    v = Ly * sin(pi/Lx * nodes(:,1));           % y-displacement 
    nodes_shifted = nodes + [v*0 v]*xi;         % new nodes
    arc_shape = zeros(numel(nodes),1);         
    arc_shape(2:2:end) = v;                     % vectorized defect-field
    U = arc_shape;                             	% defect basis
MeshDefected = Mesh(nodes_shifted);
MeshDefected.create_elements_table(elements, myElementConstructor);
MeshDefected.set_essential_boundary_condition([nset{1} nset{3}],1:2,0)
fprintf(' Arc shape: %.1f * beam thickness \n\n', xi)
nNodes = size(nodes,1);

% ASSEMBLY ________________________________________________________________
u0 = zeros( MeshNominal.nDOFs, 1);
myAssembly = Assembly(MeshDefected);
M = myAssembly.mass_matrix();
[K,~] = myAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    myAssembly.DATA.K = K;
    myAssembly.DATA.M = M;

% External force __________________________________________________________
node_dofs = MeshNominal.get_DOF_from_location([Lx/2, Ly/2]);
forced_dof = node_dofs(2);
Fext = zeros(nNodes*2, 1);
Fext( forced_dof ) = 1;
Fextc = myAssembly.constrain_vector( Fext );

% Let us also define the index of the forced dof in the constrained vector:
forced_dof_c = myAssembly.free2constrained_index( forced_dof );    


%% Eigenmodes                                                       

% Eigenvalue problem_______________________________________________________
n_VMs = 3; % first n_VMs modes with lowest frequency calculated

% Vibration Modes (VM)
Kc = myAssembly.constrain_matrix(K);
Mc = myAssembly.constrain_matrix(M);
[VM,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
VM = VM(:,ind);
for ii = 1:n_VMs
    VM(:,ii) = VM(:,ii)/max(sqrt(sum(VM(:,ii).^2,2)));
end
VM = myAssembly.unconstrain_vector(VM);

% PLOT
mod = 1;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure
PlotMesh(nodes_shifted, elementPlot, 0);
v1 = reshape(VM(:,mod), 2, []).';
S = 2*max(nodes_shifted(:,2));
PlotFieldonDeformedMesh(nodes_shifted, elementPlot, v1, 'factor', S);
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])
axis on; grid on; box on

% Damping _________________________________________________________________
alfa = 3.1;
beta = 6.3*1e-6;
D = alfa*M + beta*K;
myAssembly.DATA.D = D;
Dc = myAssembly.constrain_matrix(D);


%% Modal Derivatives                                                

tic
[MD, MDnames] = modal_derivatives(myAssembly, VM); % nominal
toc

% PLOT an MD
nwho = 1;
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure
v1 = reshape(MD(:, nwho), 2, []).';
S = Ly/2;
PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S, 'component', 'U');
title(['\theta_{' num2str(MDnames(nwho,1)) num2str(MDnames(nwho,2)) '}'])
axis on; grid on; box on


%% ROM                                                              

V = [VM MD]; 	% reduced order basis
m = size(V,2);

% mass-normalize
for ii = 1 : size(V, 2)
    V(:,ii) = V(:,ii) / sqrt(V(:,ii)'*M*V(:,ii));
end

Mr = V'*M*V; 	% reduced mass matrix

% Let us defined the ROM ReducedAssembly
ROM_Assembly = ReducedAssembly(MeshDefected, V);

% compute reduced stiffness tensors (using in-built functions)
mode = 'ELP'; %element-level projection mode
Q2 = V'*K*V;
Q3 = tensor( ROM_Assembly.tensor('T2',[m m m], [2 3], mode) );
Q4 = tensor( ROM_Assembly.tensor('T3',[m m m m], [2 3 4], mode) );

% compute tensors for the tangent stiffness matrix (see tensors_KF_NLvib)
Q3t = Q3 + permute(Q3, [1 3 2]); 
Q4t = Q4 + permute(Q4, [1 3 2 4]) + permute(Q4, [1 4 2 3]);

% check eigenfrequencies
f0_ROM = sort(sqrt(eig(Mr\Q2))/2/pi);
id = 1 : n_VMs;
f0_ROM = f0_ROM(id);
disp(table(f0, f0_ROM))


%% FRF - Harmonic Balance (with NLvib)                              

% PREPARE MODEL ___________________________________________________________
% ROMd: reduced matrices
Kr = Q2;                % reduced linear stiffness matrix
Mr = V' * M * V;        % reduced mass matrix
Dr = V' * D  * V;       % reduced damping matrix
Fr = V'*Fext;           % reduced external force vector

% Let us defined the ROMd Assembly (although the ReducedAssembly class is
% not strictly necessary, we want to define a different object - remember
% that the synthax obj1=obj2 does NOT copy an object)
ROM_Assembly = Assembly(MeshDefected, V); 
ROM_Assembly.DATA.K = Kr;
ROM_Assembly.DATA.M = Mr;
ROM_Assembly.DATA.D = Dr;
% the function to compute the nonlinear forces and the jacobian in NLvib:
fnl_CUSTOM = @(q) tensors_KF_NLvib(Q3, Q4, Q3t, Q4t, q);
ROM_Assembly.DATA.fnl_CUSTOM = fnl_CUSTOM;
% create MechanicalSystem object for NLvib
ROM_System = FE_system(ROM_Assembly, Fr, 'custom');


% ANALYSIS PARAMETERS _____________________________________________________
imod = 1;               % eigenfreq to study
omi = 2*pi*f0(imod); 	% linear eigenfrequency
n = size(V, 2);         % number of DOFs of the reduced system
H = 7;                  % harmonic order
N = 3*H+1;              % number of time samples per period
Om_s = omi * 0.80;   	% start frequency
Om_e = omi * 1.1;    	% end frequency
ds =  5;                % Path continuation step size
exc_lev = 4000;

% COMPUTE FRs _____________________________________________________________
fprintf('\n\n FRF from %.2f to %.2f rad/s \n\n', Om_s, Om_e)
r2 = cell(length(exc_lev),1);
for iex = 1 : length(exc_lev)
    % Set excitation level
    ROM_System.Fex1 = V' * Fext * exc_lev(iex);
    
    % Initial guess (solution of underlying linear system)
    Q1 = (-Om_s^2*Mr + 1i*Om_s*Dr + Kr) \ ROM_System.Fex1;
    y0 = zeros( (2*H+1)*n , 1);
    y0( n + (1:2*n) ) = [real(Q1);-imag(Q1)];
    
    % stuff for scaling
    qscl = max(abs((-omi^2*Mr + 1i*omi*Dr + Kr) \ ROM_System.Fex1));
    dscale = [y0*0+qscl; omi];
    Sopt = struct('Dscale', dscale, 'dynamicDscale', 1);
    
    % Solve and continue w.r.t. Om	
    [X, Solinfo, Sol] = solve_and_continue(y0, ...
        @(X) HB_residual(X, ROM_System, H, N, 'FRF'), ...
        Om_s, Om_e, ds, Sopt);
    
    % Interpret solver output
    r2{iex} = nlvib_decode(X, Solinfo, Sol, 'FRF', 'HB', n, H);
    
    results.FRF.HB{iex} = r2{iex};
end


%% PLOT FRs                                                         

r2 = results.FRF.HB;

figure
h = 1;
for iex = 1 : length(exc_lev)
    % 1st harmonic amplitude of the forced dof (use force_dof_c!)
    A = r2{iex}.Qre(:, :, h);
    B = r2{iex}.Qim(:, :, h);
    
    c = V * (A + 1i*B);  % project back reduced solution to full space
    c = c(forced_dof, :);
    
    W = r2{iex}.omega;
    plot(W, abs(c) / Ly, '.-', 'linewidth', 1); hold on
end
grid on
axis tight
xlabel('\omega [rad/s]')
ylabel('|Q_1| / L_y [-]')
title('FRF with HB (beam mid-span)')


% LINEAR RESPONSE
% compute the linear FRF for comparison
nw = 501;
w_linear = linspace(Om_s, Om_e, nw);
for iex = 1 : length(exc_lev)
    fr_linear = zeros(nNodes*2, nw);
    for ii = 1:nw
        w = w_linear(ii);
        frl = (-w^2*Mr + 1i*w*Dr + Kr) \ Fr * exc_lev(iex);
        fr_linear(:, ii) = V * frl;
    end
    plot(w_linear, abs(fr_linear(forced_dof, :))/Ly, 'k--')
end
drawnow




%% auxiliary function (MDs and NLvib-tensors)                       

function [MD, names] = modal_derivatives(myAssembly, Phi)
    
    n = size(Phi,1);
    n_VMs = size(Phi,2);

    K0 = myAssembly.DATA.K;
    K0 = myAssembly.constrain_matrix( K0 );

    MD = zeros(n, n_VMs*(n_VMs+1)/2);
    names = zeros(n_VMs*(n_VMs+1)/2, 2);
    kk = 1;
    for jj = 1 : n_VMs

        Phi_j = Phi(:, jj);
        
        dK_deta_j = myAssembly.stiffness_derivative(Phi_j);
        dK_deta_j = myAssembly.constrain_matrix( dK_deta_j );

        for ii = 1 : n_VMs
            if ii < jj
                continue
            end

            Phi_i = myAssembly.constrain_vector( Phi(:, ii) );
            dPhi_i_deta_j = -K0\(dK_deta_j * Phi_i); 

            th =  dPhi_i_deta_j / max(abs(dPhi_i_deta_j));
            MD(:,kk) = myAssembly.unconstrain_vector( th );
            names(kk, :) = [ii jj];
            kk = kk + 1;
        end
    end
end

function [Kt,fi] = tensors_KF_NLvib(Q3,Q4,Q3t,Q4t,q)
% NB: in NLvib, the nonlinear function must contain ONLY the nonlinear
% terms (i.e. WITHOUT the linear ones)
fi = ttsv(Q3,q,-1)  + ttsv(Q4,q,-1);
Kt = ttsv(Q3t, q, -2) + ttsv(Q4t,q,-2);
end




