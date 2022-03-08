function dz = lorenz(z,sigma,rho,beta)
% vector field of Lorenz system
dz = [sigma*(z(2)-z(1));
    rho*z(1)-z(2)-z(1)*z(3);
    -beta*z(3)+z(1)*z(2)];

end