function [M,C,K,fnl,fext] = build_model(psi,mu)
% Build model for 2 coupled nonlinear parametric amplifiers
c  = -0.02 ; 
c1 = -1*c ;
c2 =  0 ;
c3 = -1*c ;

k1 = 1;
k2 = 1;
k3 = 1;

m1 = 1;
m2 = 1;

M = [m1,0;0,m2];
C = -[-c1 , 0 ; 0 , -c3];
K = -[(-k1-k2) , k2 ; k2 , (-k2 -k3)];


%% Nonlinearity

kap1 = 0.1;
kap2 = 0;
kap3 = 0.1;

gam1 = 0.02;
gam2 = 0.0;
gam3 = 0.02;
% Stiffness
fnl(2).coeffs =    [ -kap1    , 0 ;...
                     0    , - kap3   ];

fnl(2).ind = [3 ,  0;...
              0  , 3; ...
              zeros(2,2)];

          % Damping 
fnl(2).coeffs = -[ fnl(2).coeffs ,  [ -gam1,0 ;...
                                    0, - gam3   ] ];

fnl(2).ind =  [ fnl(2).ind , [ zeros(2,2); ...
                                3  , 0;...
                                0  , 3    ] ].';         
                            
                            
                            
%% (Parametric) excitation
q = 0.2;
f0 = q * [ 0  ; 1];
Q = [( -k1) , 0 ; 0 , (-k3 )];
% Stiffness

 
fext.data = external_excitation(f0,Q, mu ,psi);
end


function [data] = external_excitation(f0,Q,mu1,psi)

f_ord = 2; %nonlinearity x including zeroth order
nKappa = 2;

idle = repmat(struct('coeffs',[],'ind',[]),f_ord,1);
data = repmat(struct('kappa',[],'F_n_k',idle),nKappa,1);

n = size(f0,1);

% External excitation
% kappa_1
data(1).kappa = 1;
data(1).f_n_k(1).coeffs = f0/2* exp(1j*psi);
data(1).f_n_k(1).ind    = sparse(1,n);

% kappa_1
data(2).kappa = -1;
data(2).f_n_k(1).coeffs = f0/2* exp(-1j*psi);
data(2).f_n_k(1).ind    = sparse(1,n);

% {
% Linear Parametric Excitation
% kappa_1
data(3).kappa = 2;
data(3).f_n_k(2).coeffs = mu1*Q/2;
data(3).f_n_k(2).ind    = speye(n);

% kappa_1
data(4).kappa = -2;
data(4).f_n_k(2).coeffs = mu1*Q/2;
data(4).f_n_k(2).ind    = speye(n);

 
%}
end
