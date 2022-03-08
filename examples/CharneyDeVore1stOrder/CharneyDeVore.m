function y = CharneyDeVore(x,x1star)

% parameter set
gamma = 1;
beta  = 1.25;
b     = 0.5;
C     = 0.1;
m = 1:2;
alpham    = 8*sqrt(2)/pi*m.^2./(4*m.^2-1).*(b^2+m.^2-1)./(b^2+m.^2);
betam     = beta*b^2./(b^2+m.^2);
deltam    = 64*sqrt(2)/(15*pi)*(b^2-m.^2+1)./(b^2+m.^2);
gammaastm = gamma*4*m*sqrt(2)*b./(4*m.^2-1)/pi;
gammam    = gamma*4*m.^3./(4*m.^2-1)*sqrt(2)*b./pi./(b^2+m.^2);
epsilon   = 16*sqrt(2)/(5*pi);

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
x4 = x(4,:);
x5 = x(5,:);
x6 = x(6,:);
x4star = x1star*(-0.4);

y = zeros(6,numel(x1));
y(1,:) = gammaastm(1)*x3-C*(x1-x1star);
y(2,:) = -(alpham(1)*x1-betam(1)).*x3-C*x2-deltam(1)*x4.*x6;
y(3,:) = (alpham(1)*x1-betam(1)).*x2-gammam(1)*x1-C*x3+deltam(1)*x4.*x5;
y(4,:) = gammaastm(2)*x6-C*(x4-x4star)+epsilon*(x2.*x6-x3.*x5);
y(5,:) = -(alpham(2)*x1-betam(2)).*x6-C*x5-deltam(2)*x4.*x3;
y(6,:) = (alpham(2)*x1-betam(2)).*x5-gammam(2)*x4-C*x6+deltam(2)*x4.*x2;

end