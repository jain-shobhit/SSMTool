function [J] = frc_Jacobian(rho,psi,gamma,lambda,epsilon,R1)
% FRC_JACOBIAN
% This function evaluates the jacobian of polar reduced dynamics wrt the
% parametrisation coordinates, for a 2D SSM
%
% [J] = FRC_JACOBIAN(rho,psi,gamma,lambda,epsilon,R1)
%
% rho:      parametrisation coordinate polar amplitude
% psi:      parametrisatino coordinate polar phase omega*t
% gamma:    reduced dynamics coefficients
% lambda:   master mode eigenvalue
% epsilon:  excitation amplitude
% R1:       non-autonomous reduced dynamics coefficients
%
% J:        Jacobian of 2D polar (nonaut) reduced dynamics evaluated on (rho,psi)
%
% See also: FRC_LEVEL_SET, CHECK_STABILITY

f = nonzeros(R1(1).R(1).coeffs);

if isempty(f)
    f = 0;
end
kappa = R1(1).kappa;

if kappa <0
    f = conj(f);
    kappa = -kappa;
end

kappa0 = kappa;


c = epsilon*(real(f)*cos(psi) + imag(f) * sin(psi));
d = epsilon*(-real(f)*sin(psi) + imag(f) * cos(psi));

J = [real(lambda),  d;
        -d/(rho^2), -c/rho];

for ell = 1:length(gamma)
    J(1,1)  = J(1,1) + real(gamma(ell)) * (2*ell + 1) * rho^(2*ell);
    J(2,1)  = J(2,1) + imag(gamma(ell)) *  2*ell * rho^(2*ell-1);
end


% Higher order non-autonomous contributions
J = J_nonaut(epsilon,J,R1, rho, psi, kappa0);
end


function [J] = J_nonaut(epsilon,J,R_1, RHO, PSI, kappa0)

% assuming reduced dynamics of form (rhodot,psidot)
J_11 = J(1,1);
J_21 = J(2,1);
J_12 = J(1,2);
J_22 = J(2,2);

num_kappa = numel(R_1);
for i = 1:num_kappa
    

    kappa = R_1(i).kappa/kappa0;
    cos_k = cos(kappa*PSI);
    sin_k = sin(kappa*PSI);
    
    for k = 1:(numel(R_1(i).R)-1)
        Rk = R_1(i).R(k+1);

        [~,col,r] = find(Rk.coeffs(1,:));
        if any(col)
        midx = sum(Rk.ind(col,:),2); %exponent of spatial component
        
        run_idx = 1;
        for col_j = col
            rj = r(run_idx); % reduced dynamics coefficient in column col_j
            
            J_11 = J_11 + epsilon*RHO.^(midx(run_idx)-1).*midx(run_idx).*... % drhodot/drho
                ( real(rj)* cos_k + imag(rj)*sin_k );            
                        
            J_21 = J_21 + epsilon*RHO.^(midx(run_idx)-2).*(midx(run_idx)-1).*... % dpsidot/drho
                (imag(rj)*cos_k  - real(rj)*sin_k  );  
            
            J_12 = J_12 + epsilon*RHO.^(midx(run_idx)).*kappa.*... % drhodot/dpsi
                ( - real(rj)* sin_k + imag(rj)*cos_k);
            
            J_22 = J_22 - epsilon*RHO.^(midx(run_idx)-1).*kappa.*...  % dpsidot/dpsi
                ( imag(rj)*sin_k  + real(rj)*cos_k );
                
            run_idx = run_idx+1;
        end
        end
    end
end    
J = [J_11,J_12; J_21, J_22];
end
