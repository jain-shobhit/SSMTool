function [aji] = axial_forcing(N,l)

syms x
% N = 2; % number of modes
% l = 2;

% construct modal functions
phis   = [];
dphis  = [];
ddphis = [];
alphal = zeros(N,1);
alphal(1:4) = [3.927 7.069 10.210 13.352]';
alphal(5:N) = (4*(5:N)'+1)*pi/4;
alpha  = alphal/l;
R = sin(alphal)./sinh(alphal);
E = (0.5*l*(1-R.^2)+(R.^2.*sinh(2*alphal)-sin(2*alphal))/(4*alpha)).^(-0.5);
for n=1:N
   phin   = E(n)*(sin(alpha(n)*x)-R(n)*sinh(alpha(n)*x));
   dphin  = diff(phin,x);
   ddphin = diff(dphin,x);
   phis   = [phis; phin];
   dphis  = [dphis; dphin];
   ddphis = [ddphis; ddphin];
end


%% Build Parametric Excitation coefficients
aji = zeros(N,N);

for i = 1:N
    for j = 1:N
        aji(j,i) = int( phis(j)* ddphis(i),x, 0,l);
    end
end

%% Check orthonormality
%{
for n = 1:N
   spphi = int(phis(1)*phis(n), x, 0, l);
   fprintf('Scalarproduct <phi%d phi%d> is %d\n',1, n, spphi);
end
%}
end
    
                
                