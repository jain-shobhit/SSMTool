function [a, b] = frc_ab(rho,omega,gamma,lambda)
a = rho * real(lambda);
b = rho * (imag(lambda) - omega);

for j = 1:length(gamma)
    a = a + real(gamma(j))* rho.^(2*j+1);    
    b = b + imag(gamma(j)) * (rho.^(2*j+1));
end

end
