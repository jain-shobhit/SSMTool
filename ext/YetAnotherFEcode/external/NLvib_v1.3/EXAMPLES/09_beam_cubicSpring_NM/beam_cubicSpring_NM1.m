%========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of a cantilevered Euler-Bernoulli beam
% with cubic spring at its free end, in accordance with [F. Thouverez:
% Presentation of the ECL benchmark, MSSP 17 (1)(2003)195–202].
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
len = .7;               % length
height = .014;          % height in the bending direction
thickness = .014;       % thickness in the third dimension
E = 2.05e11;            % Young's modulus
rho = 7800;             % density
BCs = 'clamped-free';   % constraints

% Setup one-dimensional finite element model of an Euler-Bernoulli beam
n_nodes = 20;           % number of equidistant nodes along length
beam = FE_EulerBernoulliBeam(len,height,thickness,E,rho,...
    BCs,n_nodes);

% Apply cubic spring element at free end in translational direction
inode = n_nodes;
dir = 'trans';
knl = 6e9;
add_nonlinear_attachment(beam,inode,dir,'cubicSpring','stiffness',knl);

% Number of degrees of freedom
n = beam.n;
%% Modal analysis the linearized system
[PHI_lin,OM2] = eig(beam.K,beam.M);
om_lin = sqrt(diag(OM2));
[om_lin,ind] = sort(om_lin); PHI_lin = PHI_lin(:,ind);
%% Nonlinear modal analysis using harmonic balance
analysis = 'NMA';

% Analysis parameters
H = 5;              % harmonic order
N = 2^6;            % number of time samples per period
imod = 1;           % mode to be analyzed
log10a_s = -5;      % start vibration level (log10 of modal mass)
log10a_e = -3;      % end vibration level (log10 of modal mass)
inorm = n-1;        % coordinate for phase normalization

% Initial guess vector x0 = [Psi;om;del], where del is the modal
% damping ratio, estimate from underlying linear system
om = om_lin(imod); phi = PHI_lin(:,imod);
Psi = zeros((2*H+1)*n,1);
Psi(n+(1:n)) = phi;
x0 = [Psi;om;0];
qscl = max(abs(Psi));

% Solve and continue w.r.t. Om
ds = .6;
Sopt = struct('Dscale',[qscl*1e-2*ones(size(x0,1)-2,1);om;1e0;1e0],...
    'dynamicDscale',1);
fscl = mean(abs(beam.K*phi));
X_HB = solve_and_continue(x0,...
    @(X) HB_residual(X,beam,H,N,analysis,inorm,fscl),...
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
energy = zeros(size(a_HB));
for i=1:size(X_HB,2)
    Qi = reshape(Q_HB(:,i),n,2*H+1);
    q0 = Qi(:,1) + sum(Qi(:,2:2:end),2);
    u0 = sum(Qi(:,3:2:end),2)*om_HB(i);
    energy(i) = 1/2*u0'*beam.M*u0 + 1/2*q0'*beam.K*q0 + knl*q0(n-1)^4/4;
end
%% Illustrate nonlinear modal characteristics

% Modal frequency vs. amplitude
figure; hold on;
plot(log10(energy),om_HB/(2*pi),'k-o');
xlabel('log10(energy)'); ylabel('modal frequency in Hz');
set(gca,'ylim',[20 50]);

figure
plot(om_HB/(2*pi), log10(a_HB) , 'k-o');
xlabel('modal frequency in Hz');
ylabel('a_{HB}')
hold on
plot(xlim, log10a_s*[1 1],'r--')
plot(xlim, log10a_e*[1 1],'r--')
grid on
