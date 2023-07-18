function stability = check_stability(rho,psi,gamma,lambda,epsilon,R1)
% CHECK_STABILITY
% This function checks the stability of periodic orbits computed from 2D
% polar reduced dynamics by evaluation of the Jacobian of the hyperbolic
% fixedpoints of their polar odes.
%
% stability = CHECK_STABILITY(rho,psi,gamma,lambda,epsilon,R1)
%
% rho:      parametrisation coordinate polar amplitude
% psi:      parametrisatino coordinate polar phase omega*t
% gamma:    reduced dynamics coefficients
% lambda:   master mode eigenvalue
% epsilon:  excitation amplitude
% R1:       non-autonomous reduced dynamics coefficients
%
% stability:array that contains stability information for all periodic
%           orbits that are input
%
% See also: FRC_JACOBIAN, FRC_LEVEL_SET

nRho = length(rho);
stability = zeros(size(rho));
for j = 1:nRho
%% 
% The stability is assesed by the eigenvalues of the Jacobian 
% 
% $$J(\rho)=\left[\begin{array}{cc}\partial_{\rho}a(\rho) & -\rho\left[b(\rho)-m\Omega\right]\\\partial_{\rho}b(\rho)+\frac{\left[b(\rho)-m\Omega\right]}{\rho} 
% & \frac{a(\rho)}{\rho}\end{array}\right]$$
    J = frc_Jacobian(rho(j),psi(j),gamma,lambda,epsilon,R1);
    trJ = trace(J);
    detJ = det(J);
%% 
% Routh Stability criterion
    if detJ>0 && trJ<0 % asymptotically stable
        stability(j) = 1;
    else % not asymptotically stable
        stability(j) = 0;
    end
end
stability = logical(stability);
end