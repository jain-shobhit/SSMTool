function [a,b,da,db] = cal_ab_dab(rho,gamma,lambda)
% As indicated by the name of this function, it returns the results of a, b
% and their derivatives w.r.t. rho

a  = rho * real(lambda);
b  = imag(lambda);
da = real(lambda);
db = 0;

for j = 1:length(gamma)
    a = a + real(gamma(j))* rho.^(2*j+1);    
    b = b + imag(gamma(j)) * (rho.^(2*j));
    da = da + (2*j+1)*real(gamma(j))* rho.^(2*j);
    db = db + (2*j)*imag(gamma(j)) * (rho.^(2*j-1));
end

end