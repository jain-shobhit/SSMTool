function [mass,damp,stiff,fnl,fext] = build_model(n,u,beta,miu,Gamma,alpha,BC,varargin)

l = 1;
disp('Getting linear and nonlinearity coefficients')
fileName = [BC,'_tensors_',num2str(n),'_',num2str(1),'.mat'];
try 
    load(fileName,'a_ijkl','b_ijkl','c_ijkl','g_ijkl','a','b','delta','delta1','intphi');   
    disp('Loaded coefficients from storage');
catch
    disp('Calculating coefficients');
    [a,b,delta,delta1,a_ijkl,b_ijkl,c_ijkl,g_ijkl,intphi] = cal_parameters(n,l,BC);
    disp('Saving coefficients');
    save(fileName,'a_ijkl','b_ijkl','c_ijkl','g_ijkl','a','b','delta','delta1','intphi','-v7.3')
end

mass  = delta;
stiff = delta1+(u^2-Gamma)*b;
alpha_ijkl = u^2*c_ijkl+a_ijkl;
beta_ijkl  = 2*u*sqrt(beta)*b_ijkl;
gamma_ijkl = g_ijkl;
if numel(varargin)>0
    damp_type = varargin{1};
    if strcmp(damp_type, 'nonlinear_damp')
      % visoelastic damping
          damp = alpha*delta1+2*u*sqrt(beta)*a;
          tmp1 = sptensor(alpha_ijkl);
          tmp2 = sptensor(beta_ijkl);
          tmp3 = sptensor(gamma_ijkl);
          subs1 = tmp1.subs; 
          subs2 = tmp2.subs;
          subs3 = tmp3.subs;
          if ~isempty(subs2)
            subs2(:,4) = subs2(:,4)+n; %dq_l
          end
          subs3(:,3) = subs3(:,3)+n; %dq_k
          subs3(:,4) = subs3(:,4)+n; %dq_l
          subs = [subs1; subs2; subs3];
          vals = [tmp1.vals; tmp2.vals; tmp3.vals];
          f3 = sptensor(subs, vals, [n,2*n,2*n,2*n]); 
    else
        % viscoelastic damping but ignore nonlinear part of q*dq*q
        damp = lamda*delta1+2*u*sqrt(beta)*a;
        f3 = sptensor(alpha_ijkl);
    end
else
    % propotional damping c\dot{w}
    damp = 2*u*sqrt(beta)*a;
    f3 = -0.5*miu*sptensor(alpha);
end
        
fnl = {[],f3};
fext = intphi;
end