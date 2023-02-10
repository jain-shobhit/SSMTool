function y = slide_pend_ode(t,x,p,c1,c2,k1,k2)

% parameters
m1 = 1;
m2 = 1;
len = 1;
J2 = m2*len^2/12;
g = 9.8;
% c1 = 0.17;
% c2 = 0.02;
% k1 = 7.45;
% k2 = 1;

% states and parameters
x1  = x(1,:);
phi = x(2,:);
v1  = x(3,:);
phiDot = x(4,:);

epf = p(1,:);
om  = p(2,:);

% ODE formulation
detm = (m1+m2)*(J2+0.25*m2*len^2)-0.25*m2^2*len^2*(cos(phi)).^2;
fx1  = -c1*v1-k1*x1+epf.*cos(om.*t)+0.5*m2*len*sin(phi).*phiDot.^2;
fphi = -c2*phiDot-k2*phi-0.5*len*m2*g*sin(phi);
y(1,:) = v1;
y(2,:) = phiDot;
y(3,:) = (J2+0.25*m2*len^2)*fx1-0.5*m2*len*cos(phi).*fphi;
y(4,:) = -0.5*m2*len*cos(phi).*fx1+(m1+m2)*fphi;
y(3,:) = y(3,:)./detm;
y(4,:) = y(4,:)./detm;

end