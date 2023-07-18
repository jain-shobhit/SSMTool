%========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of a single-degree-of-freedom (SDOF)
% oscillator with cubic spring, i.e. the Duffing oscillator, under harmonic
% external forcing and light linear viscous damping. The dynamics are
% governed by the second-order ordinary diferential equation,
% 
%   mu \ddot q + zeta \dot q + kappa q + \gamma q^3 = P cos( Om * t ).
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
%% Parameters of the Duffing oscillator
mu = 1;
zeta = 0.05;
kappa = 1;
gamma = 0.1;
P = 0.18;
%% Compute frequency response using numerical implementation of the 
% harmonic balance method

% Analysis parameters
H = 7;      % harmonic order
N = 2^6;    % number of time samples per period
Om_s = .5;  % start frequency
Om_e = 1.6; % end frequency

% Initial guess (from underlying linear system)
Q = (-Om_s^2*mu+1i*Om_s*zeta+kappa)\P;
x0 = [0;real(Q);-imag(Q);zeros(2*(H-1),1)];

% Solve and continue w.r.t. Om
ds = .01;                       % Path continuation step size
Sopt = struct('jac','none');    % No analytical Jacobian provided here
X = solve_and_continue(x0,...
    @(X) HB_residual_Duffing(X,mu,zeta,kappa,gamma,P,H,N),...
    Om_s,Om_e,ds,Sopt);

% Determine excitation frequency and amplitude (magnitude of fundamental
% harmonic)
Om = X(end,:);
a = sqrt(X(2,:).^2 + X(3,:).^2);
%% Analytically calculate frequency response using single-term harmonic 
% balance.

% In this case, the excitation frequencies are determined depending on the
% amplitude. Hence, we specify the amplitude range and samples for which we
% determine the frequencies.
a_min = .1;
a_max = 3;
a_ana = linspace(a_min,a_max,50);

% For each amplitude, determine the two associated excitation frequencies
Om_ana = zeros(length(a_ana),2);
for i=1:length(a_ana)
    Om_ana(i,1) = sqrt(1-zeta^2/2+3*gamma*a_ana(i)^2/4 + ...
        sqrt(P^2/a_ana(i)^2+zeta^4/4-zeta^2-3*zeta^2*gamma*a_ana(i)^2/4));
    Om_ana(i,2) = sqrt(1-zeta^2/2+3*gamma*a_ana(i)^2/4 - ...
        sqrt(P^2/a_ana(i)^2+zeta^4/4-zeta^2-3*zeta^2*gamma*a_ana(i)^2/4));
end
% Only the real-valued solutions exist. Let us store the information which
% ones are valid.
valid_ana = imag(Om_ana(:,1))==0 & imag(Om_ana(:,2))==0;
%% Illustrate results
figure; hold on;
plot(Om,a,'k-');
plot(Om_ana(valid_ana,:),a_ana(valid_ana),'gx');
set(gca,'xlim',[Om_s Om_e]);
xlabel('excitation frequency \Omega'); ylabel('response amplitude a');
legend(['numerical HB, H=' num2str(H)],'analytical HB, H=1',...
    'LOCATION','NW');