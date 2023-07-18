function psi = frc_psi(rho,omega,gamma,lambda,f)
[a,b] = frc_ab(rho,omega,gamma,lambda);

psi = atan2( (rho.*b * real(f) - a*imag(f)), (-a*real(f) - rho.*b*imag(f)) );

end