function [M,C,K,fnl,fext] = build_model_phys_nD(n)

c = 0.05; % 0.05 2D - mode 3,4


k = 1;

m = 1;

g = 0.1;

M = m * speye(n);



Cel = -[-c , c  ; c  , -c ];
Kel = -[-k , k ; k , -k];

C = sparse(n,n); 
K = sparse(n,n);

for i = 1:n-1
    C(i:i+1,i:i+1) = C(i:i+1,i:i+1) +Cel;
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) +Kel;
end

% Parametric excitation matrix
Q = K;

% Copuling to wall
C(1,1) = C(1,1) + c; C(n,n) = C(n,n) +c;
K(1,1) = K(1,1) + k; K(n,n) = K(1,1) +k;

% Nonlinearity of springs


nl_elem = -[ - g  , 3*g  , -3*g  , g  ;...
                     g     , -3*g , 3*g   , -g    ];
ind_elem = [3 ,2 ,1, 0;...
              0 ,1 ,2, 3];
nl = sparse(n,3*(n-1)+1);
ind = sparse(n,3*(n-1)+1);
for i = 1:n-1
   idx = 1+ 3*(i-1);
   nl(i:i+1,idx:(idx+3)) = nl(i:i+1,idx:(idx+3)) + nl_elem;
   ind(i:i+1, idx:(idx+3) ) =  ind_elem;
end

% Coupling to wall
nl(1,1) = nl(1,1) + g;
nl(n,end) = nl(n,end) + g;

fnl(2).coeffs = nl;
fnl(2).ind =  ind.';

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