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
%% Nonlinear modal analysis using Harmonic Balance
analysis = 'NMA';

% Analysis parameters
H = 9;              % harmonic order
N = 2^7;            % number of time samples per period
imod = 1;           % mode to be analyzed
log10a_s = -1;      % start vibration level (log10 of modal mass)
log10a_e = pi;      % end vibration level (log10 of modal mass)
inorm = 2;          % coordinate for phase normalization

% Initial guess vector y0 = [Psi;om;del], where del is the modal
% damping ratio, estimate from underlying linear system
om = om_lin(imod); phi = PHI_lin(:,imod);
Psi = zeros((2*H+1)*n,1);
Psi(n+(1:n)) = phi;
x0 = [Psi;om;0];

% Solve and continue w.r.t. Om
% NOTE: Near the loops related to internal resonances, the 
% predictor-corrector continuation can be very sensitive to step length,
% scaling, parametrization and other parameters.
ds      = .013;
Sopt    = struct('Dscale',[1e-1*ones(size(x0,1)-2,1);1;1e-8;1],...
    'dynamicDscale',1);
[X_HB,Solinfo,Sol] = solve_and_continue(x0,...
    @(X) HB_residual(X,oscillator,H,N,analysis,inorm),...
    log10a_s,log10a_e,ds,Sopt);

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
%% Illustrate Frequency-Energy-Plot

% Overview
figure; hold on;
plot(log10(energy_HB),om_HB,'k-x');
xlabel('log10(energy)'); ylabel('\omega');

% Zoom into range near loops related to internal resonances
figure; hold on;
plot(log10(energy_HB),om_HB,'k-x');
xlabel('log10(energy)'); ylabel('\omega');
set(gca,'xlim',[2.4 6.1],'ylim',[1.375 sqrt(2)]);