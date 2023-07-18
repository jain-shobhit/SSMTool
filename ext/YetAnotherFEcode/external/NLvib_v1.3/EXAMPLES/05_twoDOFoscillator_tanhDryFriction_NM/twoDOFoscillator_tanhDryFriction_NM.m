%========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of a two-degree-of-freedom (2DOF)
% oscillator with regularized Coulomb dry friction nonlinearity, similar to
% that in Laxalde, Thouverez 2009.
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
mi = [.02 1];           % masses
ki = [0 40 600];        % springs
di = 0*ki;              % no dampers

% Elastic Coulomb element applied to 1st DOF
eps = .05;
muN = 5;
W = [1;0];
nonlinear_elements = struct('type','tanhDryFriction',...
    'eps',eps,'friction_limit_force',muN,'force_direction',W);

% Define system as chain of oscillators
oscillator = ChainOfOscillators(mi,di,ki,...
    nonlinear_elements);

% Number of DOFs
n = oscillator.n;
%% Modal analysis the linearized system

% Modes for free sliding contact
[PHI_free,OM2] = eig(oscillator.K,oscillator.M);
om_free = sqrt(diag(OM2));
% Sorting
[om_free,ind] = sort(om_free); PHI_free = PHI_free(:,ind);

% Modes for fixed contact
inl = find(W); B = eye(length(oscillator.M)); B(:,inl) = [];
[PHI_fixed,OM2] = eig(B'*oscillator.K*B,B'*oscillator.M*B);
om_fixed = sqrt(diag(OM2));
% Sorting
[om_fixed,ind] = sort(om_fixed); PHI_fixed = B*PHI_fixed(:,ind);
%% Nonlinear modal analysis using harmonic balance
analysis = 'NMA';

% Analysis parameters
H = 21;             % harmonic order
Ntd = 2^10;         % number of time samples per period
imod = 1;           % mode to be analyzed
log10a_s = -1.5;    % start vibration level (log10 of modal mass)
log10a_e = 1;       % end vibration level (log10 of modal mass)
inorm = 2;          % coordinate for phase normalization

% Initial guess vector x0 = [Psi;om;del], where del is the modal
% damping ratio, estimate from underlying linear system
om = om_fixed(imod); phi = PHI_fixed(:,imod);
Psi = zeros((2*H+1)*n,1);
Psi(n+(1:n)) = phi;
x0 = [Psi;om;0];

% Solve and continue w.r.t. Om
ds = .01;
fscl = mean(abs(oscillator.K*phi));
X_HB = solve_and_continue(x0,...
    @(X) HB_residual(X,oscillator,H,Ntd,analysis,inorm,fscl),...
    log10a_s,log10a_e,ds);

% Interpret solver output
Psi_HB = X_HB(1:end-3,:);
om_HB = X_HB(end-2,:);
del_HB = X_HB(end-1,:);
log10a_HB = X_HB(end,:);
a_HB = 10.^log10a_HB;
Q_HB = Psi_HB.*repmat(a_HB,size(Psi_HB,1),1);
%% Nonlinear modal analysis using shooting

% Determine qs(inorm) for harmonic balance results
log10qsinorm_HB = log10(sum(real(Q_HB(inorm:n:end,:))));

% Limits of nonlinear modal analysis
log10qsinorm_s = -1.5;      % start vibration level (log10 of qs(inorm))
log10qsinorm_e = 1;         % end vibration level (log10 of qs(inorm))

% Number of time samples per period
Ntd = 2^10;

% Initial guess vector x0 = [qs_/a;us_/(a*om);om;del], where del is the
% modal damping ratio and ()_ means all components except that of the
% coordiante 'inorm' which serves as normalization, estimate from 
% underlying linear system
nnorm = setdiff(1:n,inorm);
phi = PHI_free(:,imod);
qs = phi/phi(inorm); us = 0*phi;
x0 = [qs(nnorm);us(nnorm);om_fixed(imod);0];
qscl = max(abs(qs(nnorm)));

% Perform continuation on Om
ds = .02;
Np = 1; % we seek period-one solutions
X_shoot = solve_and_continue(x0,...
    @(X) shooting_residual(X,oscillator,Ntd,Np,analysis,qscl,[],inorm),...
    log10qsinorm_s,log10qsinorm_e,ds);

% Interpret solver output
om_shoot = X_shoot(end-2,:);
del_shoot = X_shoot(end-1,:);
log10qsinorm_shoot = X_shoot(end,:);
%% Compare modal characteristics for Shooting and Harmonic Balance methods

% Modal frequency vs. amplitude
figure; hold on;
plot(log10qsinorm_HB,om_HB/om_fixed(imod),'g-');
plot(log10qsinorm_shoot,om_shoot/om_fixed(imod),'k--');
xlabel('log10(q^s(i_{norm}))'); ylabel('\omega/\omega_0');
legend('HB','Shooting');

% Modal damping ratio vs. amplitude
figure; hold on;
plot(log10qsinorm_HB,del_HB*1e2,'g-');
plot(log10qsinorm_shoot,del_shoot*1e2,'k--');
xlabel('log10(q^s(i_{norm}))'); ylabel('modal damping ratio in %');
legend('HB','Shooting');