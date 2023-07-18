function [J] = ode_2DSSM_cartesian_DFDP(t, x, p, data,obj)
% ODE_2DSSM_CARTESIAN_DFDP
% This function computes the jacobian of reduced dynamics in its time dependent normal
% form in cartesian coordinates for a 2D SSM with respect to the parameters omega and epsilon. 
% This function is not vectorised as of now.
%
% [y] = ODE_2DSSM_CARTESIAN_DFDP(t, x, p, data,obj)
%
% t:             time variable
% x:             phase space coordinates
% p:             array containing continuation parameters omega and epsilon
% data contains: order of RD
%                R - autonomous reduced dynamics coefficients
%                W - autonomous SSM coefficients
% obj            SSM object
%
% J:             result of evaluating Jacobian of ode J_p(t,x,p)
%
% See also: ODE_2DSSM_CARTESIAN_DFDX, ODE_2DSSM_CARTESIAN             


order = data.order;
om = p(1,:);

%% Omega dependence of nonautonomous coefficients is taken into account recompute coefficients for each new frequency
% set correct omega
if  isempty(obj.R_1)|| om(1) ~=  obj.R_1(1).Omega
    % As nonautonomous reduced dynamics implicitly depend on omega
    % compute them again
    obj.System.Omega = om(1);
    [X, S] = obj.compute_perturbed_whisker(order-1,data.W,data.R);
    
    obj.R_1 = S; % Nonautonomous reduced dynamics
    obj.W_1 = X; % Nonautonomous SSM coefficients
end

% Autonomous part
        
J = aut_J_dp();


% Nonautonomous part
if obj.Options.contribNonAuto
    %{
    if obj.FRCOptions.omDepNonAuto
        %% Omega dependence of nonautonomous coefficients is taken into account recompute coefficients for each new frequency
        % set correct omega
        % Compute sensitivity coefficients
        if isempty(obj.DR_1) || om(1) ~=  obj.DR_1(1).Omega
            
            [~,DS] = obj.compute_sensitivity_coefficients(order-1,data.W,data.R,obj.W_1);
            obj.DR_1 = DS;
            
        end
    else
        %% Weak dependence of non-autonomous coefficients is assumed, only compute them once for exact resonance        
        if  isempty(obj.DR_1)
            [~,DS] = obj.compute_sensitivity_coefficients(order-1,data.W,data.R,obj.W_1);
            obj.DR_1 = DS;
        end
    end
    %}
    
    J = nonaut_J_dp(t,x,p,obj.R_1,J);

else
    J = lead_nonaut_J_dp(t,p,obj.R_1,J);
end

end



function [J] = aut_J_dp()
% Autonomous part of the jacobian
% no dependence on omega and epsilon
J = zeros(2);
end


function [J] = nonaut_J_dp(t,x,p,S,J)
% Nonautonomous part of the reduced dynamics including leading order
J_11 = J(1,1); %dx(1,:)dotdOmega
J_21 = J(2,1); %dx(2,:)dotdOmega
J_12 = J(1,2); %dx(1,:)dotdepsilon
J_22 = J(2,2); %dx(2,:)dotdepsilon

%parametrisation coordinates
q1 = x(1,:) + 1i * x(2,:);
q2 = x(1,:) - 1i * x(2,:); %conjugate of q1

om  = p(1,:);
eps = p(2,:);

%% Leading order
% Reduced Dynamics Coefficients
[~,~,s] = find(S(1).R(1).coeffs(1,:));
if ~isempty(s)
exp_kap = exp(1i* S(1).kappa * om.*t);

J_11 = J_11 + eps.*(real(1i*S(1).kappa.*t.* s .* exp_kap)); 
J_21 = J_21 + eps.*(imag(1i*S(1).kappa.*t.*s .* exp_kap));
J_12 = J_12 + (real(s .* exp_kap)); 
J_22 = J_22 + (imag(s .* exp_kap));
end

%% Higher orders
% {
num_kappa = numel(S);
for i = 1:num_kappa %each harmonic
    kappa = S(i).kappa;
    exp_kap = exp(1i* kappa * om.*t);
    for k = 1:(numel(S(i).R)-1) %every spatial expansion order
        Sk = S(i).R(k+1); %order k coefficients
        
        [~,col,s] = find(Sk.coeffs(1,:));
        
        if any(col)
            exp_k = Sk.ind(col,:); %exponent of spatial component
            run_idx = 1;
            for s_j = s
                
                % Spatial contribution 
                m1 = exp_k(run_idx,1);
                m2 = exp_k(run_idx,2);
                spatial_part = q1.^m1.*q2.^m2;

                % Frequency contribution with coeffs
                freq_part    = s_j*exp_kap;
                i_kappa_t    = 1i*kappa*t; % derivative coefficient form exp_kap
                
                J_11 = J_11 + eps.*(real(spatial_part .* i_kappa_t.* freq_part));  %dx(1,:)dotdOmega
                J_21 = J_21 + eps.*(imag(spatial_part .* i_kappa_t.* freq_part));  %dx(2,:)dotdOmega
                
                J_12 = J_12 +      real((spatial_part .* freq_part)); %dx(1,:)dotdepsilon
                J_22 = J_22 +      imag((spatial_part .* freq_part)); %dx(2,:)dotdepsilon
                               
                run_idx = run_idx+1;
            end
        end
        
        %{
        % Sensitivity coefficients
        DSk = DS(i).R(k+1); %order k coefficients
        [~,dcol,ds] = find(DSk.coeffs(1,:));
        
        if any(dcol)
            exp_k = DSk.ind(dcol,:); %exponent of spatial component
            run_idx = 1;
            for ds_j = ds
                % Spatial contribution 
                m1 = exp_k(run_idx,1);
                m2 = exp_k(run_idx,2);
                spatial_part = q1.^m1.*q2.^m2;

                % Frequency contribution with coeffs
                freq_der     = ds_j *exp_kap;

                % Derivatives of the coefficients
                J_11 = J_11 + eps.*(real(spatial_part .* freq_der));  %dx(1,:)dotdOmega
                J_21 = J_21 + eps.*(imag(spatial_part .* freq_der));  %dx(2,:)dotdOmega
                
                run_idx = run_idx+1;
            end
        end
        %}
    end
end
%}
J = [ J_11, J_12 ...
    ; J_21, J_22];
end


function [J] = lead_nonaut_J_dp(t,p,S,J)
% Nonautonomous part of the reduced dynamics including only leading order
J_11 = J(1,1); %dx(1,:)dotdOmega
J_21 = J(2,1); %dx(2,:)dotdOmega
J_12 = J(1,2); %dx(1,:)dotdepsilon
J_22 = J(2,2); %dx(2,:)dotdepsilon


om  = p(1,:);
eps = p(2,:);

%% Leading order
% Reduced Dynamics Coefficients
[~,~,s] = find(S(1).R(1).coeffs(1,:));
if ~isempty(s)
exp_kap = exp(1i* S(1).kappa * om.*t);

J_11 = J_11 + eps.*(real(1i*S(1).kappa.*t.* s .* exp_kap)); 
J_21 = J_21 + eps.*(imag(1i*S(1).kappa.*t.*s .* exp_kap));
J_12 = J_12 + (real(s .* exp_kap)); 
J_22 = J_22 + (imag(s .* exp_kap));
end

J = [ J_11, J_12 ...
    ; J_21, J_22];
end
