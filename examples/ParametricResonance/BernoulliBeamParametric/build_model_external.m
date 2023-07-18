function [M,C,K,fnl,fext] = build_model(kappa, gamma, nElements,mu)

[M,C,K]=L_Bernoulli_Beam_Model(nElements);

n = length(M);
forcing_dof = n-1;

%fext = sparse(n,1);
%fext(forcing_dof) = 1;


fext.data = forcing(forcing_dof,n,mu);

nDOF = 2*nElements;
subs3 = [[nDOF-1,2*nDOF- 1,2*nDOF- 1,2*nDOF- 1];
         [nDOF-1,nDOF-1,nDOF-1,nDOF-1]];
vals3 = [gamma;kappa];
f3 = sptensor(subs3, vals3, [nDOF,2*nDOF,2*nDOF,2*nDOF]);


%% 
% first-order tensors

f2 = sptensor([n,n,n]);
% f3 = sptensor([n,n,n,n]);
% %% 
% % Adding cubic and damper to the middle and end node of the beam 
% 
% dof = [ n/2-1, n-1]; % only even number of dofs allowed
% for j = dof
%     f3(j,j,j,j) = kappa;
% end

fnl = {f2,f3};
end


function [data] = forcing(dof,n,mu)

f_0 = sparse(n,1);
f_0(dof) = 1;

%% External Excitation
% kappa_1
data(1).kappa = 1;
data(1).f_n_k(1).coeffs = f_0/2;
data(1).f_n_k(1).ind    = sparse(1,n);

% kappa_1
data(2).kappa = -1;
data(2).f_n_k(1).coeffs = f_0/2;
data(2).f_n_k(1).ind    = sparse(1,n);

%% Linear parametric Excitation
% kappa_3
data(3).kappa = 2;
data(3).f_n_k(2).coeffs = mu*f_0/2;
data(3).f_n_k(2).ind    = f_0.';

% kappa_4
data(4).kappa = -2;
data(4).f_n_k(2).coeffs = mu*f_0/2;
data(4).f_n_k(2).ind    = f_0.';

end