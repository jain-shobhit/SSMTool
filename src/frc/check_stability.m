function stability = check_stability(rho,psi,gamma,lambda,epsilon,f)
nRho = length(rho);
stability = zeros(size(rho));
for j = 1:nRho
%% 
% The stability is assesed by the eigenvalues of the Jacobian 
% 
% $$J(\rho)=\left[\begin{array}{cc}\partial_{\rho}a(\rho) & -\rho\left[b(\rho)-m\Omega\right]\\\partial_{\rho}b(\rho)+\frac{\left[b(\rho)-m\Omega\right]}{\rho} 
% & \frac{a(\rho)}{\rho}\end{array}\right]$$
    J = frc_Jacobian(rho(j),psi(j),gamma,lambda,epsilon,f);
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