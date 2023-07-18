%========================================================================
% DESCRIPTION: 
% Investigation of the forced and damped dynamics of a 
% two-degree-of-freedom (2DOF) oscillator with cubic spring.
% Besides the well-known occurrence of turning points in the
% amplitue-frequency curve, the system exhibits Neimark-Sacker 
% bifurcations, which is highlighted in this example.
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
mi = [1 .05];           % masses
ki = [1 .0453 0];       % springs
di = [.002 .013 0];     % dampers
Fex1 = [.11;0];         % fundamental excitation force harmonic

% Unilateral spring applied to 1st DOF
nonlinear_elements{1} = struct('type','cubicSpring',...
    'stiffness',1,'force_direction',[1;0]);
nonlinear_elements{2} = struct('type','cubicSpring',...
    'stiffness',.0042,'force_direction',[1;-1]);

% Define system as chain of oscillators
oscillator = ChainOfOscillators(mi,di,ki,...
    nonlinear_elements,Fex1);
n = oscillator.n;   % number of degrees of freedom
%% Compute frequency response using harmonic balance
analysis = 'FRF';

% Analysis parameters
H = 7;         % harmonic order
N = 2^8;       % number of time samples per period
Om_s = .8;     % start frequency
Om_e = 1.4;    % end frequency

% Initial guess (solution of underlying linear system)
Q1 = (-Om_s^2*oscillator.M + 1i*Om_s*oscillator.D + oscillator.K)\Fex1;
y0 = zeros((2*H+1)*length(Q1),1);
y0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)];

% Solve and continue w.r.t. Om
ds = .02;
Sopt = struct('Dscale',[1e-0*ones(size(y0));Om_s]);
[X_HB,Solinfo_HB] = solve_and_continue(y0,...
    @(X) HB_residual(X,oscillator,H,N,analysis),...
    Om_s,Om_e,ds,Sopt);

% Interpret solver output
Om_HB = X_HB(end,:);
Q_HB = X_HB(1:end-1,:);

% Define amplitude as magnitude of the fundamental harmonic of the second
% coordinate's displacement
a_HB = sqrt(Q_HB(n+2,:).^2 + Q_HB(2*n+2,:).^2);
a_rms_HB = sqrt(sum(Q_HB(1:2:end,:).^2))/sqrt(2);
%% Compute frequency response using shooting method

% Number of time samples per period
Ntd = 2^8;

% Initial guess (solution of underlying linear system)
ys = [real(Q1);-Om_s*imag(Q1)];

% Solve and continue w.r.t. Om
ds = .02;
Sopt = struct('Dscale',[1e0*ones(size(ys));Om_s]);%,'dynamicDscale',1);
Np = 1; % we seek period-one solutions
[X_shoot,Solinfo] = ...
    solve_and_continue(ys,...
    @(X) shooting_residual(X,oscillator,Ntd,Np,analysis),...
    Om_s,Om_e,ds,Sopt);

% Interpret solver output
Om_shoot = X_shoot(end,:);
Ys = X_shoot(1:end-1,:);
%% Determine and compare amplitude-frequency curves for both harmonic 
% balance and shooting method

% Calculate amplitude also for the results of the shooting method, and
% determine the asymptotic stability according to Floquet theory
a_shoot = zeros(size(Om_shoot));
a_rms_shoot = zeros(size(Om_shoot));
stable = zeros(size(Om_shoot));
mucrit = zeros(size(Om_shoot));
BPs(1:length(Om_shoot)) = struct('Om',[],'a',[],'a_rms',[],'type','');
nBP = 0;
tic;
for i=1:length(Om_shoot)    
    % Evaluate solution and monodromy matrix
    [~,~,~,Y,dye_dys] = shooting_residual(X_shoot(:,i),...
        oscillator,Ntd,Np,analysis);
    
    % Determine fundamental harmonic magnitude
    Qc = fft(Y(:,1))/Ntd;
    a_shoot(i) = 2*abs(Qc(2));
    a_rms_shoot(i) = sqrt(sum(abs(Qc(1:H+1)).^2))/sqrt(2)*2;
    
    % Determine stability in accordance with Floquet theory: a periodic
    % solution is stable, if all eigenvalues of the monodromy matrix remain
    % within the unit circle in the complex plane
    mucrit(i) = eigs(dye_dys,1,'lm');   % leading Floquet multiplier
    stable(i) = abs(mucrit(i))<=1; % allow for some tolerance
    
    % Analyze type of stability loss
    if i>1 && stable(i)~=stable(i-1)
        % Increment number of bifurcation points (BPs)
        nBP = nBP+1;
        
        % Store location of BP
        BPs(nBP).Om = (Om_shoot(i)+Om_shoot(i-1))/2;
        BPs(nBP).a = (a_shoot(i)+a_shoot(i-1))/2;
        BPs(nBP).a_rms = (a_rms_shoot(i)+a_rms_shoot(i-1))/2;
        
        % Identify and store type of BP
        if abs(angle(mucrit(i)))<1e-3
            % Critical Floquet multiplier appears to leave unit circle 
            % via +1. This indicates a turning point or branching
            % bifuraction. In this case example it is a turning point (TP).
            BPs(nBP).type = 'TP';
        elseif abs(abs(angle(mucrit(i)))-pi)<1e-3
            % Critical Floquet multiplier appears to leave unit circle 
            % via -1. This indicates a period doubling bifurcation (PD).
            BPs(nBP).type = 'PD';
        else
            % Critical Floquet multiplier appears to leave unit circle as
            % complex-valued pair. This indicates a Neimark-Sacker
            % bifurcation (NS).
            BPs(nBP).type = 'NS';
        end
    end
end
BPs(nBP+1:end) = [];
disp(['A posteriori Floquet stability analysis required ' ...
    num2str(toc) ' s.']);

% Plot results
figure; hold on;
plot(Om_HB,a_rms_HB,'g-');
plot(Om_shoot,a_rms_shoot,'k--.','markersize',5);
plot(Om_shoot(~stable),a_rms_shoot(~stable),'rx');
xlabel('excitation frequency'); ylabel('response amplitude');
legend('HB','Shooting','unstable','LOCATION','NE');
% Plot BPs
for i=1:nBP
    switch BPs(i).type
        case 'TP'
            styl = {'o','markersize',10,'markerfacecolor','r'};
        case 'NS'
            styl = {'^','markersize',10,'markerfacecolor','c'};
        case 'PD'
            styl = {'s','markersize',10,'markerfacecolor','m'};
    end
    plot(BPs(i).Om,BPs(i).a_rms,styl{:});
end
set(gca,'xlim',sort([Om_s Om_e]));