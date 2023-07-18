% nlvib_ShootingPostprocess
%
% Syntax:
% s = nlvib_ShootingPostprocess(X, MechSystem, Ns, Np, analysis, H, inorm)
%
% Description: this function computes the time series and spectra
% associated to the solution X of the shooting method, which contains only
% the set of initial conditions for the periodic orbit. The outputs are
% organized in a structure which should be easy to navigate.
% INPUTS:
%   X: solution from the shooting method
%   MechSystem: the FE_system for NLvib
%   Ns: number of time samples for integration over the periodic orbit
%   Np: X was sought for an Np-periodic orbit
%   analysis: 'FRF' or 'NMA'
%   H: desired output harmonics
%   inorm (optional): coordinate wrt which normalization was done during X 
%   computation in the case of NMA analysis (not required for FRF).
% OUTPUTS:
%   s: struct variable containing the results. The fields are
%       .harmonics.(Qabs/Qphase/rms/max).(displ/vel)
%       .time.(displ/vel/Nsamples)
%       .floquet.(mucrit/stable/ctime)
%
% Author: Jacopo Marconi, Politecnico di Milano
% Created: 26/04/2021
% Last modified: 27/04/2021

function s = nlvib_ShootingPostprocess(X, MechSystem, Ns, Np, analysis, H, inorm)

if nargin < 7
    inorm = [];
end

n = MechSystem.n; % number of DOFs

% Calculate amplitudes also for the results of the shooting method, and
% determine the asymptotic stability according to Floquet theory

if strcmpi(analysis, 'NMA')
    om  = X(end-2,:);                   % rad/s
    del = X(end-1,:);                   % modal damping ratio
    log10qsinorm_shoot = X(end,:);
    a = 10.^log10qsinorm_shoot;
elseif strcmpi(analysis, 'FRF')
    om  = X(end,:);                     % rad/s
end

nw = length(om);

% sizes (frequency): [nDOFs * nOmega * (1+H)]
%       (time):      [nDOFs * nOmega * Nsamples]
Qabs    = zeros(2*n, nw, 1+H);	% amplitude
Qphase  = zeros(2*n, nw, 1+H);	% phase
RMS = zeros(2*n, nw);           % rms
MAX = zeros(2*n, nw);           % max
stable = zeros(nw);
mucrit = zeros(nw);
YY = zeros(n*2, nw, Ns);

t0 = tic;
for i = 1 : nw
    % Evaluate solution and monodromy matrix
    [~,~,~,Y,dye_dys] = shooting_residual(X(:, i),...
        MechSystem, Ns, Np, analysis, [], [], inorm);
    
    % Determine first h harmonic magnitudes
    Qc = fft(Y(:,:))/Ns;
    Qabs(:, i, 1:1+H) = 2*abs( Qc(1:1+H, :) ).';
    Qabs(:, i, 1) = 0.5*Qabs(:, i, 1);
    Qphase(:, i, 1:1+H) = angle( Qc(1:1+H, :) ).';
    
    RMS(:, i) = rms(Y,1).'; % root mean square values
    MAX(:, i) = max(Y).';   % maximum values
    
    YY(:, i, :) = Y.';      % time series
    
    % Determine stability in accordance with Floquet theory: a periodic
    % solution is stable, if all eigenvalues of the monodromy matrix remain
    % within the unit circle in the complex plane
    mucrit(i) = eigs(dye_dys,1,'lm');	% leading Floquet multiplier
    if isnan(mucrit(i))
        % some times eigs do not converge is asked only 1 eigenvalue. Here
        % we ask to compute 20 (this number can be changed if needed).
        mutemp = eigs(dye_dys, min([size(dye_dys,1) 20]), 'lm');
        mucrit(i) = mutemp(1);
        fprintf(' \b\b Adjusted mucrit = %.2f + %.2f i\n', real(mucrit(i)), imag(mucrit(i)))
    end
    stable(i) = abs(mucrit(i))<=1;    	% allow for some tolerance
end
t1 = toc(t0);

s.omega = om;
if strcmpi(analysis, 'NMA')
    s.omega = om;
    s.delta = del;
    s.a = a;
end
s.harmonics.Qabs.displ   = Qabs(n+1:end, :, :);
s.harmonics.Qabs.vel     = Qabs(1:n, :, :);
s.harmonics.Qphase.displ = Qphase(n+1:end, :, :);
s.harmonics.Qphase.vel   = Qphase(1:n, :, :);
s.harmonics.rms.displ = RMS(n+1:end, :);
s.harmonics.rms.vel   = RMS(1:n, :);
s.harmonics.max.displ = MAX(n+1:end, :);
s.harmonics.max.vel   = MAX(1:n, :);
s.harmonics.info = 'sizes: [nDOFs, nOmega, 1+H]';
s.time.displ = YY(n+1:end, :, :);
s.time.vel   = YY(1:n, :, :);
s.time.Nsamples = Ns;
s.time.info = 'sizes: [nDOFs, nOmega, nTimeSamples]';
s.floquet.mucrit = mucrit;
s.floquet.stable = stable;
s.floquet.ctime = t1; % computational time

