function rho = compute_rho_grid(omegaRange,nOmega,rhoScale,gamma,lambda,nRho)
omega = linspace(omegaRange(1),omegaRange(end), nOmega);
%%
% *Explicit quadratic approximation of the backbone curve*
%
% $$\rho_{backbone} = \sqrt{\frac{\Omega-\Im(\lambda)}{\Im(\gamma_1)}}$$
rho_bb = real(sqrt((omega - imag(lambda))/imag(gamma(1))));
rhomax = rhoScale * max(rho_bb);
rho = (rhomax/nRho) * (1:nRho);
end