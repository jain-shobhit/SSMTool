function [M,C,K,fnl,f_0] = build_model(kappa, gamma, nElements)

[M,C,K]=L_Bernoulli_Beam_Model(nElements);

n = length(M);
forcing_dof = n-1;

f_0 = sparse(n,1);
f_0(forcing_dof) = 1;


nDOF = 2*nElements;
subs3 = [[nDOF-1,2*nDOF- 1,2*nDOF- 1,2*nDOF- 1];
         [nDOF-1,nDOF-1,nDOF-1,nDOF-1]];
vals3 = [gamma;kappa];
f3 = sptensor(subs3, vals3, [2*nDOF,2*nDOF,2*nDOF,2*nDOF]);


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
