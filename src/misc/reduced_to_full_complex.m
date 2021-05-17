function [z] = reduced_to_full_complex(p,W0,W1,epsilon)

nt = size(p,2);

if isempty(W1) || epsilon == 0
    N = size(W0{1}.coeffs ,1);
    z = zeros(N,1);
else
    W10 = W1{1};
    phi = linspace(0,2*pi,nt); % assuming single periodic frequency
    z =  epsilon * ( W10.coeffs * exp(1i * W10.kappas * phi));
end

for j = 1:length(W0)
    z =  z + expand_multiindex(W0{j},p); 
end

end