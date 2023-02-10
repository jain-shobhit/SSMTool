function [lamd2,lamd3] = cal_reaction_force(t,x,epf,om,c1,c2,k1,k2)

% parameters
m1 = 1;
m2 = 1;
len = 1;
J2 = m2*len^2/12;
g = 9.8;

% states and parameters
x1  = x(:,1);
phi = x(:,2);
v1  = x(:,3);
phiDot = x(:,4);


% ODE formulation
detm = (m1+m2)*(J2+0.25*m2*len^2)-0.25*m2^2*len^2*(cos(phi)).^2;
fx1  = -c1*v1-k1*x1+epf.*cos(om.*t)+0.5*m2*len*sin(phi).*phiDot.^2;
fphi = -c2*phiDot-k2*phi-0.5*len*m2*g*sin(phi);

ax1  = (J2+0.25*m2*len^2)*fx1-0.5*m2*len*cos(phi).*fphi;
aphi = -0.5*m2*len*cos(phi).*fx1+(m1+m2)*fphi;
ax1  = ax1./detm;
aphi = aphi./detm;

ax2 = ax1+0.5*len*(cos(phi).*aphi-sin(phi).*phiDot.^2);
ay2 = 0.5*len*(-sin(phi).*aphi-cos(phi).*phiDot.^2);
lamd2 = -m2*ax2;
lamd3 = m2*g-m2*ay2;



% cos(phi)*phiDot
% -sin(phi)*phiDot



end