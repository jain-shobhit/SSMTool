function [A,B,F] = build_model(sigma,rho,beta)
% linear part
A = [-sigma sigma 0;rho -1 0;0 0 -beta];
B = eye(3);
% quadratic part
F2 = sptensor([3 3 3]);
F2(2,1,3) = -1; % -xz term
F2(3,1,2) = 1;  % xy term
F = {F2};
end