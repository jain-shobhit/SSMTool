%========================================================================
% DESCRIPTION: 
% Investigation of the forced and damped dynamics of a 
% chain of three oscillators with many nonlinear elements.
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
mi = [1 1 1];           % masses
ki = [1 1 1 1];         % springs
di = .02*ki;            % dampers
Fex1 = [0;1;0];         % fundamental excitation force harmonic

% Elastic dry friction elements (aka Jenkins elements)
knl = 20; muN = 1;
W1 = [1;0;0];
nonlinear_elements{1} = struct('type','elasticDryFriction',...
    'stiffness',knl,'friction_limit_force',muN,'ishysteretic',true,...
    'force_direction',W1);
W2 = [-1;1;0];
nonlinear_elements{2} = struct('type','elasticDryFriction',...
    'stiffness',knl,'friction_limit_force',muN,'ishysteretic',true,...
    'force_direction',W2);
dKfric = knl*(W1*W1'+W2*W2');

% Cubic spring elements
W3 = [0;1;0];
nonlinear_elements{3} = struct('type','cubicSpring',...
    'stiffness',1.,'force_direction',W3);
W4 = [0;-1;1];
nonlinear_elements{4} = struct('type','cubicSpring',...
    'stiffness',1.,'force_direction',W4);

% Unilateral spring element
W5 = [0;0;1];
nonlinear_elements{5} = struct('type','unilateralSpring',...
    'stiffness',1.,'gap',.25,'force_direction',W5);

% Define system as chain of oscillators
nonlinearMDOFoscillator = ChainOfOscillators(mi,di,ki,...
    nonlinear_elements,Fex1);
%% Compute frequency response using Harmonic Balance
analysis = 'FRF';

% Analysis parameters
H = 7;      % harmonic order
N = 2^10;   % number of time samples per period
Om_s = .5;  % start frequency
Om_e = 2;   % end frequency

% Initial guess (from underlying linear system)
Q1 = (-Om_s^2*nonlinearMDOFoscillator.M + ...
    1i*Om_s*nonlinearMDOFoscillator.D + ...
    nonlinearMDOFoscillator.K + dKfric)\Fex1;
x0 = zeros((2*H+1)*length(Q1),1);
x0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)];

% Solve and continue w.r.t. Om
ds = .005;
X = solve_and_continue(x0,...
    @(Y) HB_residual(Y,nonlinearMDOFoscillator,H,N,analysis),...
    Om_s,Om_e,ds);

% Interpret solver output
n = nonlinearMDOFoscillator.n;
OM = X(end,:);
Q1 = X(1:n:end-1,:);
Q2 = X(2:n:end-1,:); 
Q3 = X(3:n:end-1,:); 

% Define amplitude as root-mean-square value
Q1_rms = sqrt(sum(Q1.^2))/sqrt(2);
Q2_rms = sqrt(sum(Q2.^2))/sqrt(2);
Q3_rms = sqrt(sum(Q3.^2))/sqrt(2);

% Illustrate amplitude frequency curves
figure; hold on;
plot(OM,Q1_rms,'b-');
plot(OM,Q2_rms,'r-');
plot(OM,Q3_rms,'g-');
xlabel('excitation frequency');
ylabel('response amplitude (RMS)');
legend('q1','q2','q3');