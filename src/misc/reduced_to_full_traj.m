function [z] = reduced_to_full_traj(t,p,W0,varargin)

if nargin == 3
    N = size(W0(1).coeffs ,1);
    z = zeros(N,1);
elseif nargin ==6
    W1 = varargin{1};
    epsilon = varargin{2};
    om = varargin{3};
    if ~isempty(W1)        
        % Nonautonomous first order time contributions
        num_kappa = numel(W1);
        phi = om*t; % assuming single periodic frequency
        
        N = size(W0(1).coeffs ,1);
        z = zeros(N,1);
        
        for i = 1:num_kappa
            order = numel(W1(i).W);
            % zeroth order
            W10 = W1(i).W(1);
            z = z + epsilon * real( W10.coeffs * exp(1i * W1(i).kappa * phi));    % Higher order time contributions
            
            for j = 1:order-1 %array starting at 0
                Wij = W1(i).W(j+1);
                z = z + epsilon * real(expand_multiindex(Wij,p) .* exp(1i * W1(i).kappa * phi));
            end
            
        end
        
    else
        N = size(W0(1).coeffs ,1);
        z = zeros(N,1);
    end
else
    error('Incorrect number of arguments: Check input')
end

for j = 1:length(W0)
    z =  z + real(expand_multiindex(W0(j),p)); 
end

end