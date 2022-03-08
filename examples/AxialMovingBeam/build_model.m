function [mass,damp,gyro,stiff,fnl,fext] = build_model(n,varargin)

A = 0.04*0.03; % m^2
I = 0.04*0.03^3/12;
rho = 7680;
E = 30e9;
% eta = 1e-5*E;
eta = 1e-4*E;
L = 1;
P = 67.5e3;

kf = sqrt(E*I/(P*L^2));
k1 = sqrt(E*A/P);
alpha = I*eta/(L^3*sqrt(rho*A*P));
gamma = 0.5128;
mass  = eye(n);
damp  = zeros(n);
gyro  = zeros(n);
stiff = zeros(n);
fext  = zeros(n,1);
cubic_coeff = zeros(n,n,n,n);
for i=1:n
    damp(i,i)  = alpha*(i*pi)^4;
    stiff(i,i) = kf^2*(i*pi)^4-(gamma^2-1)*(i*pi)^2;
    fext(i)    = (1-(-1)^i)/(i*pi);
    for j=1:n
        if j~=i
            gyro(i,j) = 4*gamma*i*j/(i^2-j^2)*(1-(-1)^(i+j));
        end
        cubic_coeff(i,j,j,i) = i^2*j^2;
    end
end

disp('the first four eigenvalues for undamped system');
lamd = eigs([zeros(n) eye(n);-mass\stiff -mass\(gyro)],4,'smallestabs')

if numel(varargin)>0 && strcmp(varargin{1}, 'nonlinear_damp')
    % visoelastic damping
    tmp  = sptensor(cubic_coeff);
    subs1 = tmp.subs;
    subs2 = subs1;
    subs2(:,3) = subs1(:,3)+n;
    subs = [subs1; subs2];
    vals = [0.25*k1^2*pi^4*tmp.vals; 0.5*alpha*k1^2/kf^2*pi^4*tmp.vals];
    f3 = sptensor(subs, vals, [n,2*n,2*n,2*n]); 
else
    % viscoelastic damping but ignore nonlinear part
    f3 = 0.25*k1^2*pi^4*sptensor(cubic_coeff);
end
fnl = {[],f3};

end