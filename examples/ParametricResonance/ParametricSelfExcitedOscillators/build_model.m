function[M_matr,C,K,fnl,fext] = build_model_order2()

[alpha,beta,gamma,mu,q,rho,M,lambda] = deal( 0.01,... %alpha
                                             0.05,... % beta
                                             0.1,...  %gamma
                                             0.2,...  %mu
                                             0.2,...  %q
                                             0.2,...  % rho
                                             0.5,...  % M
                                             4);      % lambda 

M_matr = eye(2);
C = [0,0;...
     0,-alpha];
K = [  lambda+1 , -lambda  ;...
     -lambda*M,   lambda*M];
 
 
 
 fnl(2).coeffs = [gamma, 0; ...
                  0, beta; ];
                
              
 fnl(2).ind    = [ 3,0,0,0;...
                   0,0,0,3]; 
 
 n = 2;    
%lambda=0;
 fext.data = set_forcing(q,lambda,M,mu,n);
 
 end
 
 function [data] = set_forcing(q,lambda,M,mu,n)
 % Set external forcing coefficients
% Multi-indices are stored in rows, coefficients in columns

f_ord = 2; %nonlinearity x^3 including zeroth order


idle = repmat(struct('coeffs',[],'ind',[]),f_ord,1);
data = repmat(struct('kappa',[],'f_n_k',idle),4,1);

% External forcing

data(1).kappa = 1;
data(1).f_n_k(1).coeffs = q/2 * [1;0]; 
data(1).f_n_k(1).ind    = zeros(1,2);

data(2).kappa = -1;
data(2).f_n_k(1).coeffs = q/2 * [1;0]; 
data(2).f_n_k(1).ind    = zeros(1,2);

% Parametric Excitation
data(3).kappa = 2;
data(3).f_n_k(2).coeffs = [ lambda*mu/2, -lambda*mu/2;...
                            -lambda*M*mu/2, lambda*M*mu/2 ];
data(3).f_n_k(2).ind = [1,0;...
                        0,1];

data(4).kappa = -2;
data(4).f_n_k(2).coeffs = [ lambda*mu/2, -lambda*mu/2;...
                            -lambda*M*mu/2, lambda*M*mu/2];
data(4).f_n_k(2).ind = [1,0;...
                        0,1];
 end