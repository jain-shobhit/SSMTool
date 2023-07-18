function [M,C,K,fnl,fext] = build_model(mu1,mu3)
% Build model for 2 coupled parametric duffing oscillators
c  = -0.1 ; % 0.05 2D - mode 3,4,  mode 1
c1 = -1*c ;
c2 =  1*c ;
c3 = -1*c ;

k1 = 1;
k2 = 1;
k3 = 1;

m1 = 1;
m2 = 1;

M = [m1,0;0,m2];
C = -[(-c1-c2) , c2 ; c2 , (-c2 -c3)];
K = -[(-k1-k2) , k2 ; k2 , (-k2 -k3)];


%% Nonlinearity

kap1 = 0.01;
kap2 = 0.01;
kap3 = 0.01;

gam1 = 0;
gam2 = 0.001;
gam3 = 0;
% Stiffness
fnl(2).coeffs =    [ -kap1 - kap2 , 3*kap2 , -3*kap2 , kap2 ;...
                     kap2    , -3*kap2, 3*kap2  , -kap2 - kap3   ];

fnl(2).ind = [3 ,2 ,1, 0;...
              0 ,1 ,2, 3; ...
              zeros(2,4)];
% Damping 
fnl(2).coeffs = -[ fnl(2).coeffs ,  [ -gam1 - gam2 , 3*gam2 , -3*gam2 , gam2 ;...
                                    gam2    , -3*gam2, 3*gam2  , -gam2 - gam3   ] ];

fnl(2).ind =  [ fnl(2).ind , [ zeros(2,4); ...
                                3 ,2 ,1, 0;...
                                0 ,1 ,2, 3    ] ].';         
                            
                            
                            
%% (Parametric) excitation
f0 = [ 0  ; 1];
Q = [( -k2) , k2 ; k2 , (-k2 )];
% Stiffness
f_par(2).coeffs = [  - kap2 , 3*kap2 , -3*kap2 , kap2 ;...
                     kap2    , -3*kap2, 3*kap2  , -kap2    ];

f_par(2).ind = [3 ,2 ,1, 0;...
              0 ,1 ,2, 3; ...
              zeros(2,4)].';
% mu tunes magnitude of linear to cubic parametric excitation
fext.data = external_excitation(f0,Q,f_par,mu1,mu3);
end


function [data] = external_excitation(f0,Q,fnl,mu1,mu3)

f_ord = 2; %nonlinearity x including zeroth order
nKappa = 2;

idle = repmat(struct('coeffs',[],'ind',[]),f_ord,1);
data = repmat(struct('kappa',[],'F_n_k',idle),nKappa,1);

n = size(f0,1);

% External excitation
% kappa_1
data(1).kappa = 1;
data(1).f_n_k(1).coeffs = f0/2;
data(1).f_n_k(1).ind    = sparse(1,n);

% kappa_1
data(2).kappa = -1;
data(2).f_n_k(1).coeffs = f0/2;
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

% Cubic Parametric Excitation

% kappa_1
data(3).f_n_k(4).coeffs = mu3 * fnl(2).coeffs/2;
data(3).f_n_k(4).ind    = fnl(2).ind;

% kappa_2
data(4).f_n_k(4).coeffs = mu3 * fnl(2).coeffs/2;
data(4).f_n_k(4).ind    = fnl(2).ind;
%}
end
