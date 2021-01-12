function [J] = frc_Jacobian(rho,psi,gamma,lambda,epsilon,f)

c = epsilon*(real(f)*cos(psi) + imag(f) * sin(psi));
d = epsilon*(-real(f)*sin(psi) + imag(f) * cos(psi));

J = [real(lambda),  d;
        -d/(rho^2), -c/rho];

for ell = 1:length(gamma)
    J(1,1)  = J(1,1) + real(gamma(ell)) * (2*ell + 1) * rho^(2*ell);
    J(2,1)  = J(2,1) + imag(gamma(ell)) *  2*ell * rho^(2*ell-1);
end

end