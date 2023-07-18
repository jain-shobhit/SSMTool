function [mass,damp,stiff,fnl,fext] = build_model_parametric(c,rLsq,n)

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

%% Construct coefficients of parametric excitation
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


fext.data = forcing(n,aji);

end


function [data] = forcing(n,aji)

ind_1 = eye(n);


%% External Excitation
% kappa_1
data(1).kappa = 1;

% kappa_1
data(2).kappa = -1;

%% Linear parametric Excitation
% kappa_1
data(1).f_n_k(2).coeffs = aji/2;
data(1).f_n_k(2).ind    = ind_1;

% kappa_1
data(2).f_n_k(2).coeffs = aji/2;
data(2).f_n_k(2).ind    = ind_1;

end