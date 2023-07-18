function [M,C,K,fext,fnl] = build_model(w_0sq,c,lam3,gam3, dampnl)
M = 1;
C = c;
K = w_0sq;

%% External forcing
fext.data = set_forcing_coefficients(w_0sq,gam3);


%% Nonlinearity

fnl(2).coeffs = [lam3, dampnl];
fnl(2).ind    = [3,0; 0,3];
end


function [data] = set_forcing_coefficients(w_0sq,gam3)
% Set external forcing coefficients
% Multi-indices are stored in rows, coefficients in columns

f_ord = 2; %nonlinearity x^3 including zeroth order


idle = repmat(struct('coeffs',[],'ind',[]),f_ord,1);
data = repmat(struct('kappa',[],'F_n_k',idle),2,1);

% kappa_1
data(1).kappa = 1;
data(1).f_n_k(2).coeffs = [w_0sq/2 ];
data(1).f_n_k(2).ind    = [1];

% kappa_2
data(2).kappa = -1;
data(2).f_n_k(2).coeffs = [w_0sq/2];
data(2).f_n_k(2).ind    = [1];
% {
% cubic nonlinearity 
% kappa_1
data(1).kappa = 1;
data(1).f_n_k(4).coeffs = [gam3/2];
data(1).f_n_k(4).ind    = [3];

% kappa_2
data(2).kappa = -1;
data(2).f_n_k(4).coeffs = [gam3/2];
data(2).f_n_k(4).ind    = [3];
%}
end