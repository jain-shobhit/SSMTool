%========================================================================
% DESCRIPTION: 
% Investigation of the forced and damped dynamics of a 
% two-degree-of-freedom (2DOF) oscillator with unilateral spring.
% Besides the well-known occurrence of turning points in the
% amplitue-frequency curve, the system exhibits period 
% doubling and Neimark-Sacker bifurcations, which is highlighted in this
% example.
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
mi = [1 1];             % masses
ki = [0 1 1];           % springs
di = .03*ki;            % dampers
Fex1 = [0;.1];          % fundamental excitation force harmonic

% Unilateral spring applied to 1st DOF
nonlinear_elements = struct('type','unilateralSpring',...
    'stiffness',1e2,'gap',1,'force_direction',[1;0]);

% Define system as chain of oscillators
oscillator = ChainOfOscillators(mi,di,ki,...
    nonlinear_elements,Fex1);
n = oscillator.n;   % number of degrees of freedom
%% Compute frequency response using harmonic balance
analysis = 'FRF';

% Analysis parameters
H = 21;         % harmonic order
N = 2^10;       % number of time samples per period
Om_e = .5;      % start frequency
Om_s = .8;      % end frequency

% Initial guess (solution of underlying linear system)
Q1 = (-Om_s^2*oscillator.M + 1i*Om_s*oscillator.D + oscillator.K)\Fex1;
y0 = zeros((2*H+1)*length(Q1),1);
y0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)];

% Solve and continue w.r.t. Om
ds = .0085;
Sopt = struct('Dscale',[1e-0*ones(size(y0));Om_s]);
[X_HB,Solinfo_HB] = solve_and_continue(y0,...
    @(X) HB_residual(X,oscillator,H,N,analysis),...
    Om_s,Om_e,ds,Sopt);
% return
% Interpret solver output
Om_HB = X_HB(end,:);
Q_HB = X_HB(1:end-1,:);
%% Compute frequency response using shooting method

% Number of time samples per period
Ntd = 2^12;

% Initial guess (solution of underlying linear system)
ys = [real(Q1);-Om_s*imag(Q1)];

% Solve and continue w.r.t. Om
ds = .02;
Sopt = struct('Dscale',[1e0*ones(size(ys));Om_s],'dynamicDscale',1);
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

% Define amplitude as magnitude of the fundamental harmonic of the second
% coordinate's displacement
a_HB = sqrt(Q_HB(n+2,:).^2 + Q_HB(2*n+2,:).^2);

% Calculate amplitude also for the results of the shooting method, and
% determine the asymptotic stability according to Floquet theory
a_shoot = zeros(size(Om_shoot));
stable = zeros(size(Om_shoot));
mucrit = zeros(size(Om_shoot));
tic;
for i=1:length(Om_shoot)    
    % Evaluate solution and monodromy matrix
    [~,~,~,Y,dye_dys] = shooting_residual(X_shoot(:,i),...
        oscillator,Ntd,Np,analysis);
    
    % Determine fundamental harmonic magnitude
    Qc = fft(Y(:,2))/Ntd;
    a_shoot(i) = 2*abs(Qc(2));
    
    % Determine stability in accordance with Floquet theory: a periodic
    % solution is stable, if all eigenvalues of the monodromy matrix remain
    % within the unit circle in the complex plane
    mucrit(i) = eigs(dye_dys,1,'lm');   % leading Floquet multiplier
    stable(i) = abs(mucrit(i))<=1; % allow for some tolerance
end
disp(['A posteriori Floquet stability analysis required ' ...
    num2str(toc) ' s.']);

% Select interesting points for further investigation:
label{1} = 'Point in regime beyond period doubling bifurcation';
[~,ind(1)] = max(a_shoot);
label{2} = 'Point in regime of stable period-one solutions';
[~,ind(2)] =  min(abs(Om_shoot-.635));
label{3} = 'Point in regime beyond Neimark-Sacker bifurcation';
[~,ind(3)] = min(abs(Om_shoot-.6));

%% Plot results
figure; hold on;
plot(Om_HB,a_HB,'g-');
plot(Om_shoot,a_shoot,'k.','markersize',10);
% NOTE: Some actually unstable points on the overhanging branch might be
% identified as stable, e.g. because the time discretization is still
% comparatively coarse for a reliable stability analysis.
plot(Om_shoot(~stable),a_shoot(~stable),'rx');
plot(Om_shoot(ind),a_shoot(ind),'bo');
xlabel('excitation frequency'); ylabel('response amplitude');
legend({'HB','Shooting','unstable','selected'},'Location','NorthWest');
%% Investigate the selected response points
for i=1:length(ind)
    Om = Om_shoot(ind(i));
    ys = Ys(:,ind(i));
    %% Simulate a few periods using direct time step integration
    tic
    [t,Y] = ode45(@(t,y) [y(n+1:2*n); oscillator.M\( Fex1*cos(Om*t) -...
        oscillator.K*y(1:n) - oscillator.D*y(n+1:2*n) - ...
        oscillator.nonlinear_elements{1}.force_direction*...
        oscillator.nonlinear_elements{1}.stiffness*(...
        oscillator.nonlinear_elements{1}.force_direction'*y(1:n)-...
        oscillator.nonlinear_elements{1}.gap)*double(...
        oscillator.nonlinear_elements{1}.force_direction'*y(1:n)-...
        oscillator.nonlinear_elements{1}.gap > 0) )],[0 2*pi/Om*1e2],...
        ys,odeset('maxstep',2*pi/Om/2e2));
    disp(['Conventional forward numerical integration for a single (!) ' ...
        'frequency from initial conditions of Shooting required ' ...
    num2str(toc) ' s.']);
    %% Illustrate results in phase projection
    figure; hold on; title(label{i});
    
    % Plot computed periodic solution (stable or unstable).
    % We take here the results from the shooting method, but could of
    % course also synthesize the results from the harmonic balance method
    [~,~,~,Y_shoot] = shooting_residual([ys;Om],...
        oscillator,Ntd,Np,analysis);
    plot(Y_shoot(:,1),Om*Y_shoot(:,n+1),'g-','linewidth',2);
    
    % Select last few periods of direct time step integration
    ILP = t>(t(end)-5*2*pi/Om);
    plot(Y(ILP,1),Y(ILP,n+1),'k--');
    xlabel('q_1'); ylabel('u_1');
    legend('Shooting','forward num. integration');
end