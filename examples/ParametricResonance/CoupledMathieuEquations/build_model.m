function [M,C,K,fnl,fext] = build_model_phys()

c = 0.05 ; % 0.05 2D - mode 3,4,  mode 1
c1 = 1*c;
c2 = 1*c;
c3 = 1*c;

k1 = 1;
k2 = 1;
k3 = 1;

m1 = 1;
m2 = 1;

g1 = 0.1;
g2 = 0.1;
g3 = 0.1;

M = [m1,0;0,m2];
C = -[(-c1-c2) , c2 ; c2 , (-c2 -c3) ];
K = -[(-k1-k2) , k2 ; k2 , (-k2 -k3)];

Q = [ -k2  , k2 ; k2 ,  -k2  ];


fnl(2).coeffs = -[ -g1 - g2 , 3*g2 , -3*g2 , g2 ;...
                     g2    , -3*g2, 3*g2  , -g2 - g3   ];

fnl(2).ind = [3 ,2 ,1, 0;...
              0 ,1 ,2, 3].';

fext.data = external_excitation(Q);
end


function [data] = external_excitation(Q)

f_ord = 2; %nonlinearity x including zeroth order
nKappa = 2;

idle = repmat(struct('coeffs',[],'ind',[]),f_ord,1);
data = repmat(struct('kappa',[],'F_n_k',idle),nKappa,1);

n = size(Q,1);
% kappa_1
data(1).kappa = 1;
data(1).f_n_k(2).coeffs = Q/2;
data(1).f_n_k(2).ind    = speye(n);

% kappa_1
data(2).kappa = -1;
data(2).f_n_k(2).coeffs = Q/2;
data(2).f_n_k(2).ind    = speye(n);

end
