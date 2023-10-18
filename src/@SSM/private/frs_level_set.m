function y = frs_level_set(rho,omega,epsilon,gamma,lambda,r,isbaseForce)
% FRS_LEVEL_SET This function returns the level set for forced response
% surface defined via 2D SSM

[a,b] = frc_ab(rho,omega,gamma,lambda);
if isbaseForce
    y = a.^2+b.^2-epsilon.^2.*omega.^2*r^2;
else
    y = a.^2+b.^2-epsilon.^2*r^2;
end

end