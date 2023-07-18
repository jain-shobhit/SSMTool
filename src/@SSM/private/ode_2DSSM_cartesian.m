function [y] = ode_2DSSM_cartesian(t, x, p, data,obj)
% ODE_2DSSM_CARTESIAN
% This function computes the reduced dynamics in its time dependent normal
% form in cartesian coordinates for a 2D SSM. This function is not
% vectorised as of now.
%
% [y] = ODE_2DSSM_CARTESIAN(t, x, p, data,obj)
%
% t:             time variable
% x:             phase space coordinates
% p:             array containing continuation parameters omega and epsilon
% data contains: order of RD
%                R - autonomous reduced dynamics coefficients
%                W - autonomous SSM coefficients
% obj            SSM object
%
% y:             result of applying ode: y = f(t,x,p)
%
% See also: ODE_2DSSM_CARTESIAN_DFDX, ODE_2DSSM_CARTESIAN_DFDP


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

y = aut_RD(x,data.R);
y = nonaut_RD(t,x,p,obj.R_1,y);
end



function [y] = aut_RD(x,R)
% Autonomous part of the reduced dynamics

%parametrisation coordinates
q1 = x(1,:) + 1i * x(2,:);
q2 = x(1,:) - 1i * x(2,:); %conjugate of q1

y = zeros(2,size(q1,2));
order = numel(R);
n_gamma = floor((order-1)/2);

% leading order part
[~,~,r] = find(R(1).coeffs(1,:));
y(1,:) = real(r*q1);
y(2,:) = imag(r*q1);

% higher orders
for j = 1:n_gamma
    if ~isempty(R(2*j+1).ind)
        [~, loc] = ismember([j+1,j],R(2*j+1).ind,'rows');
        gamma = R(2*j+1).coeffs(1,loc);
        
        spatial_part = q1 .^(j+1).*q2.^(j);
        y(1,:) = y(1,:) + real(gamma .* spatial_part);
        y(2,:) = y(2,:) + imag(gamma .* spatial_part); %in my notation j = m2, m1-1 = m2
    end
end
end


function [y] = nonaut_RD(t,x,p,S,y)
% Nonautonomous part of the reduced dynamics including leading order

%parametrisation coordinates
q1 = x(1,:) + 1i * x(2,:);
q2 = x(1,:) - 1i * x(2,:); %conjugate of q1
om = p(1,:);
eps = p(2,:);

% Leading order
if ~isempty(S(1).R(1).coeffs)

[~,~,s] = find(S(1).R(1).coeffs(1,:));
if ~isempty(s)
    exp_kap = exp(1i* S(1).kappa * om.*t);
    y(1,:) = y(1,:) + eps.*(real(s .* exp_kap)); 
    y(2,:) = y(2,:) + eps.*(imag(s .* exp_kap));
end
end
% {
% Higher orders
num_kappa = numel(S);
for i = 1:num_kappa %each harmonic
    kappa = S(i).kappa;
    exp_kap = exp(1i* kappa * om.*t);
    for k = 1:(numel(S(i).R)-1) %every spatial expansion order
        Sk = S(i).R(k+1); %order k-1 coefficients
        if ~isempty(Sk.coeffs)
        [~,col,s] = find(Sk.coeffs(1,:));
        if any(col)
            m = Sk.ind(col,:); %exponents of spatial component in multindices
            run_idx = 1;
            for s_j = s
                % Spatial contribution 
                m1 = m(run_idx,1);
                m2 = m(run_idx,2);
                spatial_part = q1.^m1.*q2.^m2;
                
                % Frequency contribution with coeffs
                freq_part    = s_j*exp_kap;

                y(1,:) = y(1,:) + eps.*(real(spatial_part .* freq_part)); 
                y(2,:) = y(2,:) + eps.*(imag(spatial_part .* freq_part));
                run_idx = run_idx+1;
            end
        end
        end
    end
end
%}
end
