function [z] = reduced_to_full_traj(t,p,W_0,varargin)

if nargin == 3
    N = size(W_0{1}.coeffs ,1);
    z = zeros(N,1);
elseif nargin ==6
    W_1 = varargin{1};
    epsilon = varargin{2};
    om = varargin{3};
    if ~isempty(W_1)
        W_10 = W_1{1};
        phi = om*t; % assuming single periodic frequency
        z =  epsilon*real( W_10.coeffs * exp(1i * W_10.kappas * phi));
    else
        N = size(W_0{1}.coeffs ,1);
        z = zeros(N,1);
    end
else
    error('Incorrect number of arguments: Check input')
end

for j = 1:length(W_0)
    z =  z + real(expand_multiindex(W_0{j},p)); 
end

end