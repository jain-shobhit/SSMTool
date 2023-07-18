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
%% Load system matrices (M,D,K), excitation vector (Fex1), force direction
% vector (W), tip displacement response vector (T_tip)
load('beam','M','D','K','Fex1','W','T_tip');

% Define nonlinear element (you can apply many of these)
muN = 1.5;      % friction limit force
eps = 6e-7;     % regularization parameter (tanh-approx. of signum funct.)
nonlinear_elements{1} = struct('type','tanhDryFriction',...
    'friction_limit_force',muN,'eps',eps,'force_direction',W);

% Define mechanical system
beam = MechanicalSystem(M,D,K,nonlinear_elements,Fex1);

% Analysis parameters
analysis = 'FRF';
H = 7;                      % harmonic order
N = 2^7;                    % number of time samples per period
Om_s = 370;                 % start frequency
Om_e = 110;                 % end frequency

% Initial guess (we here just provide a zero vector, although a good guess 
% could be derived from the underlying linear system)
x0 = zeros((2*H+1)*length(M),1);

% Solve and continue w.r.t. Om
ds = 5;
% Set options of path continuation
% flag = 0: no actual continuation, just sequential continuation
% stepadatp = 0: no automatic step size adjustment
Sopt = struct('flag',0,'stepadapt',0,...
    'Dscale',[1e-7*ones(size(x0));Om_s]);
X = solve_and_continue(x0,...
    @(X) HB_residual(X,beam,H,N,analysis),...
    Om_s,Om_e,ds,Sopt);

% Interpret solver output
OM = X(end,:);
Q = X(1:end-1,:);

% Determine harmonics and root-mean-square value of tip displacement
Qtip = kron(eye(2*H+1),T_tip)*Q;
Qtip_rms = sqrt(sum(Qtip.^2))/sqrt(2);

% Illustrate frequency response
figure; hold on;
plot(OM,Qtip_rms,'g-');
set(gca,'yscale','log');
xlabel('excitation frequency'); ylabel('tip displacement amplitude');