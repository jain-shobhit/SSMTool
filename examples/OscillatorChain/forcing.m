function [data] = forcing(n,f_0)
% Cosine type forcing of the oscillatory chain with forcing amplitude
% vector f_0

data(1).kappa = 1;
data(2).kappa = -1;
data(1).f_n_k(1).coeffs = f_0/2;
data(1).f_n_k(1).ind    = zeros(1,n);
data(2).f_n_k(1).coeffs = f_0/2;
data(2).f_n_k(1).ind    = zeros(1,n);

end