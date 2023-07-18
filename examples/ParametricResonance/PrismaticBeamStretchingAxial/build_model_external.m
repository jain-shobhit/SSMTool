function [mass,damp,stiff,fnl,fext] = build_model_external(c,rLsq,n,mu,f_0)

l = 1.7;
disp('Getting nonlinearity coefficients')
fileName = ['tensors_',num2str(n),'_',num2str(l),'.mat'];

try
    load(fileName,'alpha','omega','nu');   
    disp('Loaded coefficients from storage');
catch
    disp('Calculating coefficients');
    [alpha, omega, nu] = cal_parameters(n,l);
    disp('Saving coefficients');
    save(fileName,'alpha','omega','nu','-v7.3')
end

paramfileName = ['paramtensors_',num2str(n),'_',num2str(l),'.mat'];

try 
    load(paramfileName,'aji');   
    disp('Loaded coefficients from storage');
catch
    disp('Building parametric coefficients');
    aji = axial_forcing(n,l);
    disp('Saving parametric coefficients');
    save(paramfileName,'aji','-v7.3')
end


mass = eye(n);
damp = 2*c*eye(n)*rLsq;
stiff = diag(omega.^2);

f3 = -rLsq*nu*sptensor(alpha);
fnl = {[],f3};


fext.data = forcing(n,aji,mu,f_0);

end


function [data] = forcing(n,aji,mu,f_0)

ind_1 = eye(n);


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
data(3).f_n_k(2).coeffs = mu*aji/2;
data(3).f_n_k(2).ind    = ind_1;

% kappa_4
data(4).kappa = -2;
data(4).f_n_k(2).coeffs = mu*aji/2;
data(4).f_n_k(2).ind    = ind_1;

end