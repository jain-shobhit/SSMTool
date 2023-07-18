function [J] = ode_2DSSM_cartesian_fixROM_DFDX(t, x, p, data)
% ODE_2DSSM_CARTESIAN_FIXROM_DFDX
% This function computes the jacobian of reduced dynamics in its time dependent normal
% form in cartesian coordinates for a 2D SSM with respect to the parametrisation coordinates. 
% This function is not vectorised as of now.
% It assumes weak dependence of the non-autonomous RD coefficients upon
% changing omega and thus does not recompute the ROM when omega is changed.
%
% [y] = ODE_2DSSM_CARTESIAN_FIXROM_DFDX(t, x, p, data,obj)
%
% t:             time variable
% x:             phase space coordinates
% p:             array containing continuation parameters omega and epsilon
% data contains: order of RD
%                R - autonomous reduced dynamics coefficients
%                W - autonomous SSM coefficients
% obj            SSM object
%
% J:             result of evaluating Jacobian of ode J_x(t,x,p)
%
% See also: ODE_2DSSM_CARTESIAN_FIXROM_DFDP, ODE_2DSSM_CARTESIAN_FIXROM 

order = data.order;
om = p(1,:);



J = aut_J_dx(x,data.R);
J = nonaut_J_dx(t,x,p,data.S,J);
end



function [J] = aut_J_dx(x,R)


%parametrisation coordinates
q1 = x(1,:) + 1i * x(2,:);
q2 = x(1,:) - 1i * x(2,:); %conjugate of q1

order = numel(R);
n_gamma = floor((order-1)/2);

% leading order Autonomous part of the jacobian
[~,~,r] = find(R(1).coeffs(1,:));
J_11 = real(r); %dx(1,:)dotdx(1,:)
J_21 = imag(r); %dx(2,:)dotdx(1,:)
J_12 = real(1i*r); %dx(1,:)dotdx(2,:)
J_22 = imag(1i*r); %dx(2,:)dotdx(2,:)


for j = 1:n_gamma
    if ~isempty(R(2*j+1).ind)
        [~, loc] = ismember([j+1,j],R(2*j+1).ind,'rows'); %in my notation j = m2, m1-1 = m2
        gamma = R(2*j+1).coeffs(1,loc);
        
        % Derivative of spatial part using chain rule
        
        deriv_x1 =    (j+1)*q1.^(j).*q2.^(j) +    j*q1.^(j+1).*q2.^(j-1);
        deriv_x2 = 1i*(j+1)*q1.^(j).*q2.^(j) - 1i*j*q1.^(j+1).*q2.^(j-1);
        
        
        J_11 = J_11  + real(gamma .* deriv_x1 ); %dx(1,:)dotdx(1,:)
        J_21 = J_21  + imag(gamma .* deriv_x1 ); %dx(2,:)dotdx(1,:)
        
        J_12 = J_12  + real(gamma .* deriv_x2 ); %dx(1,:)dotdx(1,:)
        J_22 = J_22  + imag(gamma .* deriv_x2 ); %dx(2,:)dotdx(1,:)
    end
end

J = [ J_11, J_12 ...
    ; J_21, J_22];
end


function [J] = nonaut_J_dx(t,x,p,S,J)
% Nonautonomous part of the reduced dynamics
J_11 = J(1,1); %dx(1,:)dotdx(1,:)
J_21 = J(2,1); %dx(2,:)dotdx(1,:)
J_12 = J(1,2); %dx(1,:)dotdx(2,:)
J_22 = J(2,2); %dx(2,:)dotdx(2,:)

%parametrisation coordinates
q1 = x(1,:) + 1i * x(2,:);
q2 = x(1,:) - 1i * x(2,:); %conjugate of q1

om  = p(1,:);
eps = p(2,:);

num_kappa = numel(S);
for i = 1:num_kappa %each harmonic
    kappa = S(i).kappa;
    exp_kap = exp(1i* kappa * om.*t);
    for k = 1:(numel(S(i).R)-1) %every spatial expansion order
        Sk = S(i).R(k+1); %order k coefficients
        
        [~,col,s] = find(Sk.coeffs(1,:));
        if any(col)
            m = Sk.ind(col,:); %exponent of spatial component, multiindices in m
            run_idx = 1;
            for s_j = s
                % Spatial contribution 
                m1 = m(run_idx,1);
                m2 = m(run_idx,2);
                % m1 = 0 = m2 is leading order which is zero and not
                % considered
                
                % avoid divergences at q = 0 where m_i= 0 
                if m2 == 0                    
                    deriv_x1 =        m1*q1.^(m1-1);
                    deriv_x2 =   1i * m1*q1.^(m1-1); 

                elseif m1 == 0                     
                    deriv_x1 =        m2 .*q2.^(m2-1);
                    deriv_x2 = - 1i * m2 .*q2.^(m2-1);
                elseif m1 == 0 && m2 == 0
                    deriv_x1 = 0; 
                    deriv_x2 = 0;
                else                    
                    deriv_x1 =     m1*q1.^(m1-1).*q2.^(m2) +      m2 * q1.^(m1).*q2.^(m2-1);
                    deriv_x2 = 1i* m1*q1.^(m1-1).*q2.^(m2) - 1i * m2 * q1.^(m1).*q2.^(m2-1); 
                    
                end

                % Frequency contribution with coeffs
                freq_part    = s_j*exp_kap;
                
                J_11 = J_11 + eps.*(real(deriv_x1.* freq_part));  %dx(1,:)dotdx(1,:)
                J_21 = J_21 + eps.*(imag(deriv_x1.* freq_part));  %dx(2,:)dotdx(1,:)
                
                J_12 = J_12 + eps.*(real(deriv_x2.* freq_part));  %dx(1,:)dotdx(2,:)
                J_22 = J_22 + eps.*(imag(deriv_x2.* freq_part));  %dx(2,:)dotdx(2,:)
                run_idx = run_idx+1;
            end
        end
    end
end
J = [ J_11, J_12 ...
    ; J_21, J_22];

end