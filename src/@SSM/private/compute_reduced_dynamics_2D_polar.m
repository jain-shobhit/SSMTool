function [rhodot, rhopsidot, eta] = compute_reduced_dynamics_2D_polar(RHO,PSI, lambda, gamma, R_1,Omega,epsilon)
%  COMPUTE_REDUCED_DYNAMICS_2D_POLAR: This functions computed the reduced
%  dynamics on 2D-manifolds over a polar grid (RHO,PSI). The equations are
%  as follows:
% $$\left[\begin{array}{c}\dot{\rho}\\\rho\dot{\psi}\end{array}\right]	=\mathbf{r}(\rho,\psi,\Omega):=\left[\begin{array}{c}a(\rho)\\b(\rho,\Omega)\end{array}\right]+\left[\begin{array}{cc}\cos\psi 
% & \sin\psi\\-\sin\psi & \cos\psi\end{array}\right]\left[\begin{array}{c}\mathrm{Re}\left(f\right)\\\mathrm{Im}\left(f\right)\end{array}\right],\dot{\phi}	
% =\Omega$$
% 
% where 
% 
% $$a(\rho)	=\mathrm{Re}\left(\rho\lambda+\sum_{\ell\in\mathbb{N}}\gamma_{\ell}\rho^{2\ell+1}\right),\\b(\rho,\Omega)	
% =\mathrm{Im}\left(\rho\lambda+\sum_{\ell\in\mathbb{N}}\gamma_{\ell}\rho^{2\ell+1}\right)-\eta\rho\Omega,\\f	
% =\epsilon\mathbf{u}^{\star}\mathbf{F}_{\eta}^{ext},\\\psi	=\theta-\eta\phi$$
% 
% 
R_10 = R_1{1}.coeffs;
f = epsilon*R_10(1,2);
eta = R_1{1}.kappas(2); % forcing harmonic present in the normal form
if eta < 0 % ensure the positive multiple
    eta = -eta;
    f = conj(f);
end
% fprintf(' %d times the forcing frequency %.4d is nearly resonant with the eigenvalue %.4d + i%.4d \n', eta, Omega, real(lambda),imag(lambda))
% Autonomous components of the reduced dynamics
[a,b] = frc_ab(RHO, eta*Omega, gamma, lambda);
% Reduced dynamics
rhodot = a + cos(PSI) * real(f) + sin(PSI) * imag(f);
rhopsidot = b + (cos(PSI) * imag(f) - sin(PSI) * real(f));
    
end