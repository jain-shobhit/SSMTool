%========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of a mass connected via two orthogonal 
% springs undergoing large deflections. The equations of motion of the
% autonomous system with viscous linear damping read
% 
%   \ddot q1 + 2*zt1*om1*\dot q1 + om1^2*q1 + om1^2/2*(3*q1^2+q2^2) ...
%           + om2^2*q1*q2 + (om1^2+om2^2)/2*q1*(q1^2+q2^2) = 0
%   \ddot q2 + 2*zt2*om2*\dot q2 + om2^2*q2 + om2^2/2*(3*q2^2+q1^2) ...
%           + om1^2*q1*q2 + (om1^2+om2^2)/2*q2*(q1^2+q2^2) = 0.
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

% Fundamental parameters
om1 = 1.13;
om2 = 2;
zt1 = 1e-3;
zt2 = 5e-3;

% Properties of the underlying linear system
M = eye(2);
D = diag([2*zt1*om1 2*zt2*om2]);
K = diag([om1^2 om2^2]);

% Description of the polnomial nonlinearity, with the parameters p and E
% fnl = E'*z, where z(i) = q1^(p(i,1)) * q2^(p(i,2)) * q3^(p(i,3)) * ...
% 
% In our case, we want to describe 
% fnl(1)=om1^2/2*(3*q1^2+q2^2)+om2^2*q1*q2+(om1^2+om2^2)/2*q1*(q1^2+q2^2)
% fnl(2)=om2^2/2*(3*q2^2+q1^2)+om1^2*q1*q2+(om1^2+om2^2)/2*q2*(q1^2+q2^2)
% 
p = [2 0; 1 1; 0 2; 3 0; 2 1; 1 2; 0 3];
E = [3*om1^2/2 om2^2/2; om2^2 om1^2; om1^2/2 3*om2^2/2; ... % quadratic,
    (om1^2+om2^2)/2*[1 0; 0 1; 1 0; 0 1] ...                % cubic terms
    ];

% Fundamental harmonic of external forcing
Fex1 = [1;1];

% Define oscillator as system with polynomial stiffness nonlinearities
oscillator = System_with_PolynomialStiffnessNonlinearity(M,D,K,p,E,Fex1);

% Number of degrees of freedom
n = oscillator.n;
%% Nonlinear modal analysis using harmonic balance
analysis = 'NMA';

% Analysis parameters
H = 7;              % harmonic order
N = 2^10;           % number of time samples per period
imod = 1;           % mode to be analyzed
log10a_s = -6;      % start vibration level (log10 of modal mass)
log10a_e = -.15;    % end vibration level (log10 of modal mass)
inorm = 1;          % coordinate for phase normalization

% Initial guess vector x0 = [Psi;om;del], where del is the modal
% damping ratio, estimate from underlying linear system
om = sqrt(K(imod,imod)); phi = zeros(n,1); phi(imod) = 1;
Psi = zeros((2*H+1)*n,1);
Psi(n+(1:n)) = phi;
x0 = [Psi;om;0];

% Solve and continue w.r.t. Om
ds = .01;
Sopt.termination_criterion = {@(Y) Y(end-2)<1.0};   % stop if om<1.0
X = solve_and_continue(x0,...
    @(X) HB_residual(X,oscillator,H,N,analysis,inorm),...
    log10a_s,log10a_e,ds);

% Interpret solver output
Psi_NM = X(1:end-3,:);
om_NM = X(end-2,:);
del_NM = X(end-1,:);
log10a_NM = X(end,:);
a_NM = 10.^log10a_NM;
Q_NM = Psi_NM.*repmat(a_NM,size(Psi_NM,1),1);

% Redefine amplitude as magnitude of the fundamental harmonic of the 
% first coordinate's displacement
a_NM = sqrt(Q_NM(n+1,:).^2 + Q_NM(2*n+1,:).^2);
%% Compute frequency response using harmonic balance
analysis = 'FRF';

% Analysis parameters
Om_e = .8;      % start frequency
Om_s = 1.6;     % end frequency

% Excitation levels
exc_lev = [2e-4 5e-4 1e-3 3e-3];
Om = cell(size(exc_lev)); a = cell(size(exc_lev));
for iex=1:length(exc_lev)
    % Set excitation level
    oscillator.Fex1 = Fex1*exc_lev(iex);
    
    % Initial guess (solution of underlying linear system)
    Q1 = (-Om_s^2*M + 1i*Om_s*D + K)\oscillator.Fex1;
    y0 = zeros((2*H+1)*length(Q1),1);
    y0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)];
    qscl = max(abs((-om1^2*M + 1i*om1*D + K)\oscillator.Fex1));
    
    % Solve and continue w.r.t. Om
    ds = .005;
    Sopt = struct('Dscale',[qscl*ones(size(y0));Om_s]);
    X = solve_and_continue(y0,...
        @(X) HB_residual(X,oscillator,H,N,analysis),...
        Om_s,Om_e,ds,Sopt);
    
    % Interpret solver output
    Om{iex} = X(end,:);
    Q = X(1:end-1,:);
    
    % Define amplitude as magnitude of the fundamental harmonic of the 
    % first coordinate's displacement
    a{iex} = sqrt(Q(n+1,:).^2 + Q(2*n+1,:).^2);
end
%% Illustrate frequency responses
figure; hold on;
for iex=1:length(a)
    plot(Om{iex},a{iex},'g-');
end
plot(om_NM,a_NM,'k-');
plot([om1 om1],[0 1],'k--');
set(gca,'xlim',[1.085 1.15],'ylim',[0 .35]);
xlabel('excitation frequency');
ylabel('response amplitude |Q_{1,1}|');