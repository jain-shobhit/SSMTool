function [a, b] = frc_ab(rho,omega,gamma,lambda)
% FRC_AB
% This function evaluates the autonomous reduced dynamics on the
% parametrisation coordinates, for a 2D SSM
%
% [a, b] = FRC_AB(rho,omega,gamma,lambda)
%
% rho:      parametrisation coordinate polar amplitude
% omega:    forcing frequency
% gamma:    reduced dynamics coefficients
% lambda:   master mode eigenvalue
%
% a:        contribution to the amplitude ODE of the reduced dynamics
% b:        contribution to the angle ODE of the reduced dynamics
%
% See also: COMPUTE_REDUCED_DYNAMICS_2D_POLAR

a = rho * real(lambda);
b = rho * (imag(lambda) - omega);

for j = 1:length(gamma)
    a = a + real(gamma(j))* rho.^(2*j+1);    
    b = b + imag(gamma(j)) * (rho.^(2*j+1));
end

end
