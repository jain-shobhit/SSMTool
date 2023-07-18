%========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of a two-degree-of-freedom (2DOF)
% oscillator with cubic spring, which is a common benchmark problem, see
% e.g. Kerschen et al. 2009.
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

% Parameters of the underlying linear system
mi = [1 1];         % masses
ki = [1 1 1];       % springs
di = 0*ki;          % no dampers

% Cubic spring applied to 1st DOF
nonlinear_elements = struct('type','cubicSpring',...
    'stiffness',.5,'force_direction',[1;0]);

% Define system as chain of oscillators
oscillator = ChainOfOscillators(mi,di,ki,...
    nonlinear_elements);

% Number of DOFs
n = oscillator.n;
%% Linear modal analysis
[PHI_lin,OM2] = eig(oscillator.K,oscillator.M);
[om_lin,ind] = sort(sqrt(diag(OM2)));
PHI_lin = PHI_lin(:,ind);
%% Nonlinear modal analysis of both modes
figure; hold on;
for imod = [2 1] % reverse order to plot invariant manifold of mode 1 later
%% Nonlinear modal analysis using Harmonic Balance
analysis = 'NMA';

% Analysis parameters
H = 7;              % harmonic order
N = 2^7;           % number of time samples per period
log10a_s = -1;      % start vibration level (log10 of modal mass)
log10a_e = 1.5;      % end vibration level (log10 of modal mass)
inorm = 1;          % coordinate for phase normalization

% Initial guess vector y0 = [Psi;om;del], where del is the modal
% damping ratio, estimate from underlying linear system
om = om_lin(imod); phi = PHI_lin(:,imod);
Psi = zeros((2*H+1)*n,1);
Psi(n+(1:n)) = phi;
x0 = [Psi;om;0];

% Solve and continue w.r.t. Om
ds = .02;
X_HB = solve_and_continue(x0,...
    @(X) HB_residual(X,oscillator,H,N,analysis,inorm),...
    log10a_s,log10a_e,ds);

% Interpret solver output
Psi_HB = X_HB(1:end-3,:);
om_HB = X_HB(end-2,:);
del_HB = X_HB(end-1,:);
log10a_HB = X_HB(end,:);
a_HB = 10.^log10a_HB;
Q_HB = Psi_HB.*repmat(a_HB,size(Psi_HB,1),1);

% Determine total energy in the system from the displacement and velocity
% at t=0
energy_HB = zeros(size(a_HB));
for i=1:size(X_HB,2)
    Qi = reshape(Q_HB(:,i),n,2*H+1);
    q0 = Qi(:,1) + sum(Qi(:,2:2:end),2);
    u0 = sum(Qi(:,3:2:end),2)*om_HB(i);
    energy_HB(i) = 1/2*u0'*oscillator.M*u0 + 1/2*q0'*oscillator.K*q0 + ...
        nonlinear_elements.stiffness*q0(1)^4/4;
end
%% Nonlinear modal analysis using Shooting

% Determine qs(inorm) for harmonic balance results
log10qsinorm_HB = log10(abs(sum(real(Q_HB(inorm:n:end,:)))));

% Limits of nonlinear modal analysis
log10qsinorm_s = -2;        % start vibration level (log10 of qs(inorm))
if imod==1
    log10qsinorm_e = .8;    % end vibration level (log10 of qs(inorm))
else
    log10qsinorm_e = 1;     % end vibration level (log10 of qs(inorm))
end
% Number of time samples per period
Ntd = 200;

% Initial guess vector x0 = [qs_/a;us_/(a*om);om;del], where del is the
% modal damping ratio and ()_ means all components except that of the
% coordiante 'inorm' which serves as normalization, estimate from 
% underlying linear system
nnorm = setdiff(1:n,inorm);
qs = phi/phi(inorm); us = 0*phi;
x0 = [qs(nnorm);us(nnorm);om_lin(imod);0];
qscl = mean(abs(qs(nnorm)));

% Solve and continue w.r.t. Om
ds = .02;
Np = 1; % we seek period-one solutions
X_shoot = solve_and_continue(x0,...
    @(X) shooting_residual(X,oscillator,Ntd,Np,analysis,qscl,[],inorm),...
    log10qsinorm_s,log10qsinorm_e,ds);

% Interpret solver output
om_shoot = X_shoot(end-2,:);
del_shoot = X_shoot(end-1,:);
log10qsinorm_shoot = X_shoot(end,:);

% Determine total energy in the system from the displacement and velocity
% at t=0
energy_shoot = zeros(size(om_shoot));
for i=1:size(X_shoot,2)
    q0 = zeros(n,1);
    q0(inorm) = 10^(log10qsinorm_shoot(i));
    q0(nnorm) = q0(inorm)*X_shoot(1:length(nnorm),i);
    u0 = zeros(n,1);
    u0(nnorm) = q0(inorm)*X_shoot(length(nnorm)+(1:length(nnorm)),i)*om_shoot(i);
    energy_shoot(i) = 1/2*u0'*oscillator.M*u0 + 1/2*q0'*oscillator.K*q0 + ...
        nonlinear_elements.stiffness*q0(1)^4/4;
end
%% Illustrate Frequency-Energy-Plot
plot(log10(energy_HB),om_HB,'g-');
plot(log10(energy_shoot),om_shoot,'k--');
end
xlabel('log10(energy)'); ylabel('\omega');
legend('HB','Shooting','LOCATION','NW');
set(gca,'xlim',[-2 3],'ylim',[0 5]);
%% Illustrate invariant manifold

% Select energy levels for illustration
ompl = linspace(om_HB(1), 1.3, 8);

figure; hold on;
for ipl=1:length(ompl)
    % Find energy level closest to desired frequency
    [~,i] = min(abs(om_HB-ompl(ipl)));
    
    % Synthesize time evolution of the coordinates
    tau = linspace(0,2*pi,2^8+1);
    Q1 = Q_HB(1:n:end,i); Q1 = [Q1(1);Q1(2:2:end)-1i*Q1(3:2:end)];
    Q2 = Q_HB(2:n:end,i); Q2 = [Q2(1);Q2(2:2:end)-1i*Q2(3:2:end)];
    q1 = real(exp(1i*tau(:)*(0:H))*Q1);
    u1 = real(exp(1i*tau(:)*(0:H))*(1i*om_HB(i)*(0:H)'.*Q1));
    q2 = real(exp(1i*tau(:)*(0:H))*Q2);
    plot3(q1,u1,q2);
end
view([42 -31]);
xlabel('q_1'); ylabel('u_1'); zlabel('u_2');