%========================================================================
% DESCRIPTION: 
% Matlab function setting up the frequency-domain residual vector 'R' and 
% its derivatives for given frequency 'Om' and vector of harmonics of the 
% generalized coordiantes 'Q'. The corresponding model is a single DOF
% oscillator with cubic spring nonlinearity, i.e. the Duffing oscillator,
% governed by the time-domain equation of motion
% 
%       mu * \ddot q + zeta * \dot q + kappa * q + gamma * q^3 = P * cos(Om*t).
% 
% Compared to the more general variant 'HB_residual', we
%       - consider only a single-DOF oscillator, instead of multi-DOF
%       systems
%       - consider only a cubic spring nonlinearity, instead of various
%       nonlinear elements
%       - carry out only frequency response analysis, instead of e.g. 
%       nonlinear modal analysis, and
%       - do not calculate analytical derivatives.
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
function R = HB_residual_Duffing(X,mu,zeta,kappa,gamma,P,H,N)
% Conversion of sine-cosine to complex-exponential representation
Q_ce = [flipud(X(2:2:end -1)+1i*X(3:2:end -1))/2; ...
    X(1); ...
    (X(2:2:end -1) -1i*X(3:2:end -1))/2];

% Excitation frequency
Om = X(end);

% P is the magnitude of the cosine forcing
Fex_ce = [zeros(H-1,1);P/2;0;P/2; zeros(H-1 ,1)];

% Determine inverse discrete Fourier transform matrix
E_NH = exp(1i*2*pi/N*(0:N-1)'*(-H:H));

% Apply inverse discrete Fourier transform
q = real(E_NH*Q_ce);

% Evaluate nonlinear force in the time domain
fnl = gamma*q.^3;

% Apply discrete Fourier transform
Fnl_ce = E_NH'/N*fnl;

% Dynamic force equilibrium
R_ce = ( -((-H:H)'*Om).^2 * mu + 1i*(-H:H)'*Om * zeta + kappa ).* Q_ce+...
    Fnl_ce-Fex_ce;

% Conversion of complex-exponential to sine-cosine representation
R = [real(R_ce(H+1));real(R_ce(H+2: end));-imag(R_ce(H+2: end))];