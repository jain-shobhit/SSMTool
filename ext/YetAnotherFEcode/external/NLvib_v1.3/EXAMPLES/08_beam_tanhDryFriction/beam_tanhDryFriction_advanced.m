%========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of an Euler-Bernoulli beam
% with regularized Coulomb dry friction nonlinearity.
% 
% Two variants 'simple' and 'advanced' are available, which do essentially
% the same analysis, but are different in terms of degree of complexity of
% user interaction.
% The 'simple' version loads the system matrices and vectors and sets up
% the analysis in a rather simple way.
% The 'advanced' version makes use of the class 'FE_EulerBernoulliBeam',
% analyzes the linearized variants of the system (free and sticking
% friction contact), defines damping and initial guess accordingly, and
% also does a Shooting analysis (in addition to Harmonic Balance).
%========================================================================
% This file is part of NLvib.
% 
% If you use NLvib, please refer to the book:
%   M. Krack, J. Gross: Harmonic Balance for Nonlinear Vibration
%   Problems. Springer, 2019. https://doi.org/10.1007/978-3-030-14023-6.
% 
% COPYRIGHT AND LICENSING: 
% NLvib Version 1.3 Copyright (C) 2020  Malte Krack  
%										(malte.krack@ila.uni-stuttgart.de) 
%                     					Johann Gross 
%										(johann.gross@ila.uni-stuttgart.de)
%                     					University of Stuttgart
% This program comes with ABSOLUTELY NO WARRANTY. 
% NLvib is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%========================================================================
clearvars;
close all;
clc;
addpath('../../SRC');
addpath('../../SRC/MechanicalSystems');
%% Define system

% Properties of the beam
len = 2;                % length
height = .05*len;       % height in the bending direction
thickness = 3*height;   % thickness in the third dimension
E = 185e9;              % Young's modulus
rho = 7830;             % density
BCs = 'clamped-free';   % constraints

% Setup one-dimensional finite element model of an Euler-Bernoulli beam
n_nodes = 9;            % number of equidistant nodes along length
beam = FE_EulerBernoulliBeam(len,height,thickness,E,rho,...
    BCs,n_nodes);

% Apply elastic Coulomb element at node 4 in translational direction
inode = 4;
dir = 'trans';
eps = [];
muN = 1.5;
add_nonlinear_attachment(beam,inode,dir,'tanhDryFriction',...
    'eps',eps,'friction_limit_force',muN);

% Vectors recovering deflection at tip and nonlinearity location
T_nl = beam.nonlinear_elements{1}.force_direction';
T_tip = beam.L(end-1,:);
%% Modal analysis the linearized system

% Modes for free sliding contact
[PHI_free,OM2] = eig(beam.K,beam.M);
om_free = sqrt(diag(OM2));
% Sorting
[om_free,ind] = sort(om_free); PHI_free = PHI_free(:,ind);

% Modes for fixed contact
inl = find(T_nl); B = eye(length(beam.M)); B(:,inl) = [];
[PHI_fixed,OM2] = eig(B'*beam.K*B,B'*beam.M*B);
om_fixed = sqrt(diag(OM2));
% Sorting
[om_fixed,ind] = sort(om_fixed); PHI_fixed = B*PHI_fixed(:,ind);
%% Nonlinear frequency response analysis using Harmonic Balance
analysis = 'FRF';

% Analysis parameters
H = 13;                     % harmonic order
N = 2^7;                    % number of time samples per period
Om_s = 1.5*om_fixed(1);     % start frequency
Om_e = .9*om_free(1);       % end frequency

% Apply forcing to free end of beam in translational direction, with
% magnitude 'fex'
inode = n_nodes;
dir = 'trans';
fex = .2;
add_forcing(beam,inode,dir,fex);

% Specify stiffness proportional damping corresponding to D=1% at linear 
% at the first linearized resonance
beta   = 2*1e-2/om_fixed(1);
beam.D = beta*beam.K;

% Initial guess (from underlying linear system)
Fex1 = beam.Fex1;
Q1 = B*((B'*(-Om_s^2*beam.M + 1i*Om_s*beam.D + beam.K)*B)\(B'*Fex1));
Q1_free = (-Om_s^2*beam.M + 1i*Om_s*beam.D + beam.K)\Fex1;
qscl = mean(abs(Q1));
x0 = zeros((2*H+1)*size(Q1,1),1);
x0(size(Q1,1)+(1:2*size(Q1,1))) = [real(Q1);-imag(Q1)];

% Set tanh regularization parameter
beam.nonlinear_elements{1}.eps = om_free(1)*abs(T_nl*Q1_free);

% Solve and continue w.r.t. Om
ds = 5;
% Set options of path continuation
% flag = 0: no actual continuation, just sequential continuation
% stepadatp = 0: no automatic step size adjustment
Sopt = struct('flag',0,'stepadapt',0,...
    'Dscale',[1e0*qscl*ones(size(x0));Om_s]);
X_HB = solve_and_continue(x0,...
    @(X) HB_residual(X,beam,H,N,analysis),...
    Om_s,Om_e,ds,Sopt);

% Interpret solver output
OM_HB = X_HB(end,:);
Q_HB = X_HB(1:end-1,:);

% Determine harmonics and root-mean-square value of tip displacement
Qtip_HB = kron(eye(2*H+1),T_tip)*Q_HB;
Qtip_rms_HB = sqrt(sum(Qtip_HB.^2))/sqrt(2);

% Determine linear references
OM = linspace(min(OM_HB),max(OM_HB),1e3);
Qtip_fixed = zeros(size(OM));
Qtip_free = zeros(size(OM));
for i=1:length(OM)
    Q1tmp = B*((B'*(-OM(i)^2*beam.M + 1i*OM(i)*beam.D + beam.K)*B)\...
        (B'*Fex1));
    Qtip_fixed(i) = abs(T_tip*Q1tmp)/sqrt(2);
    Q1tmp = (-OM(i)^2*beam.M + 1i*OM(i)*beam.D + beam.K)\...
        (Fex1);
    Qtip_free(i) = abs(T_tip*Q1tmp)/sqrt(2);
end
%% Nonlinear frequency response analysis using Shooting

% Number of time samples per period
Ntd = 2^9;

% Initial guess (solution of linearized system)
ys = [real(Q1);-imag(Q1)];

% Solve and continue w.r.t. Om
ds = 5;
% Set options of path continuation, see above
Sopt = struct('flag',0,'stepadapt',0, ...
    'Dscale',[qscl*1e2*ones(size(ys));Om_s]);
Np = 1; % we seek period-one solutions
X_shoot = solve_and_continue(ys,...
    @(X) shooting_residual(X,beam,Ntd,Np,analysis,qscl),...
    Om_s,Om_e,ds,Sopt);

% Determine root-mean-square value of tip displacement
n = beam.n;
OM_shoot = X_shoot(end,:);
Qtip_rms_shoot = zeros(size(X_shoot,2),1);
for i=1:size(X_shoot,2)
    [~,~,~,Y] = shooting_residual(X_shoot(:,i),beam,Ntd,Np,analysis);
    Qtip_rms_shoot(i) = rms(Y(:,1:n)*T_tip');
end
%% Compare frequency response of both Harmonic Balance and Shooting methods
figure; hold on;
plot(OM,Qtip_fixed,'-','color',.75*[1 1 1]);
plot(OM_HB,Qtip_rms_HB,'g-');
plot(OM_shoot,Qtip_rms_shoot,'k--');
set(gca,'yscale','log');
xlabel('excitation frequency'); ylabel('tip displacement amplitude');
legend('fixed contact','HB','Shooting');
%% Comparison against direct time step integration

% Select resonance
[~,ind] = max(Qtip_rms_HB);
Om = OM_HB(ind);

% Synthesize time history of harmonic balance
tau = linspace(0,2*pi,N+1)';
DFT = [ones(N+1,1) zeros(N+1,2*H)];
DFT(:,2:2:end) = cos(tau*(1:H));
DFT(:,3:2:end) = sin(tau*(1:H));
q_HB = DFT*reshape(Q_HB(:,ind),beam.n,2*H+1)';
DER = zeros(2*H,1);
DER(2:2:end) = 1:H;
DER = ( diag(DER,1) - diag(DER,-1) )*Om;
u_HB = DFT*DER*reshape(Q_HB(:,ind),beam.n,2*H+1)';

% Time interval and discretization
ts = 0; te = 2*pi/Om*20;
maxstep = 2*pi/Om/20;%2*pi/om_sticking(end)/5;

% Initial values
q0 = zeros(beam.n,1);
u0 = zeros(beam.n,1);

% Translate properties
M = beam.M;
D = beam.D;
K = beam.K;
eps = beam.nonlinear_elements{1}.eps;
friction_limit_force = muN;

% Run simulation using ODE45
tic
[t,Y] = ode45(@(t,y) [y(n+1:2*n);M\(Fex1*cos(Om*t) - K*y(1:n) - ...
    D*y(n+1:2*n) - T_nl'*muN*tanh(T_nl*y(n+1:2*n)/eps) )],...
    [ts te],[q0;u0]);%,odeset('maxstep',maxstep));
disp(['Conventional direct time step integration for a single (!) ' ...
    'frequency from homogeneous initial conditions required ' ...
    num2str(toc) ' s.']);

% Select results of last few periods
ILP = t>(t(end)-5*2*pi/Om);

% Phase portrait of tip deflection
figure; hold on;
plot(T_tip*q_HB',T_tip*u_HB','g-','linewidth',2);
plot(T_tip*Y(ILP,1:n)',T_tip*Y(ILP,n+1:2*n)','k--');
xlabel('q_{tip}'); ylabel('u_{tip}');
legend('HB','time integration');

% Phase portrait of deflection at nonlinearity location
figure; hold on;
plot(T_nl*q_HB',T_nl*u_HB','g-','linewidth',2);
plot(T_nl*Y(ILP,1:n)',T_nl*Y(ILP,n+1:2*n)','k--');
xlabel('q_{fric}'); ylabel('u_{fric}');
legend('HB','time integration');