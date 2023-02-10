function y = pend_chain_ode(t,x,p,m1,m,c1,c,k1,k,len,n)

% m1 = 1;
% m = 1;
% len = 1;
J = m*len^2/12;
g = 9.8;

ndof  = 2+3*(n-1);
qidx  = 1:ndof;
vidx  = ndof+1:2*ndof;
xidx  = qidx([1 3:3:ndof-2]);
yidx  = qidx([2 4:3:ndof-1]);
phi_idx  = [0 qidx(5:3:end)]; % add 0 for (d)phi1 (coding purpose)
dphi_idx = [0 vidx(5:3:end)];

% inverse of mass matrix
Minv = diag([1/m1; 1/m1; repmat([1/m; 1/m; 1/J],[n-1,1])]);
% g vector, G matrix and Gdot matrix
gconstraints    = zeros(2*n-1,1);
gconstraints(1) = x(yidx(1));
gconstraints(2) = x(xidx(2))-0.5*len*sin(x(phi_idx(2)))-x(xidx(1));
gconstraints(3) = x(yidx(2))-0.5*len*cos(x(phi_idx(2)))-x(yidx(1));
for i=2:n-1
    gconstraints(2*i)   = x(xidx(i+1))-0.5*len*sin(x(phi_idx(i+1)))-x(xidx(i))-0.5*len*sin(x(phi_idx(i)));
    gconstraints(2*i+1) = x(yidx(i+1))-0.5*len*cos(x(phi_idx(i+1)))-x(yidx(i))-0.5*len*cos(x(phi_idx(i)));
end

G      = zeros(2*n-1,ndof);
G(1,2) = 1;
G(2,1) = -1; G(2,3) = 1; G(2,5) = -0.5*len*cos(x(phi_idx(2)));
G(3,2) = -1; G(3,4) = 1; G(3,5) = 0.5*len*sin(x(phi_idx(2)));
for i=2:n-1
    cosphi = cos(x(phi_idx(i:i+1)))';
    sinphi = sin(x(phi_idx(i:i+1)))';
    G(2*i,[xidx(i:i+1) phi_idx(i:i+1)])   = [-1 1 -0.5*len*cosphi];
    G(2*i+1,[yidx(i:i+1) phi_idx(i:i+1)]) = [-1 1 0.5*len*sinphi];
end

Gdot      = zeros(2*n-1,ndof);
Gdot(2,5) = 0.5*len*sin(x(phi_idx(2)))*x(dphi_idx(2));
Gdot(3,5) = 0.5*len*cos(x(phi_idx(2)))*x(dphi_idx(2));
for i=2:n-1
    cosphi = cos(x(phi_idx(i:i+1)))';
    sinphi = sin(x(phi_idx(i:i+1)))';
    phidot = x(dphi_idx(i:i+1))';
    Gdot(2*i,phi_idx(i:i+1))   = 0.5*len*sinphi.*phidot;
    Gdot(2*i+1,phi_idx(i:i+1)) = 0.5*len*cosphi.*phidot;
end

% F vector
F    = zeros(ndof,1);
F(1) = p(1)*cos(p(2)*t)-k1*x(xidx(1))-c1*x(vidx(1));
F(2) = m1*g;
F(4:3:end-1) = m*g;
F(5) = -k*(2*x(phi_idx(2))-x(phi_idx(3)))-c*(2*x(dphi_idx(2))-x(dphi_idx(3)));
F(8:3:end-3) = -k*(2*x(phi_idx(3:n-1))-x(phi_idx(2:n-2))-x(phi_idx(4:n)))...
    -c*(2*x(dphi_idx(3:n-1))-x(dphi_idx(2:n-2))-x(dphi_idx(4:n)));
F(end) = -k*(x(phi_idx(n))-x(phi_idx(n-1)))-c*(x(dphi_idx(n))-x(dphi_idx(n-1)));

% Equation of motion
alpha = 1; beta = 1; % stabilization
gdotconstraints = G*x(vidx);
c = -Gdot*x(vidx)-alpha*gconstraints-beta*gdotconstraints;
y = zeros(2*ndof,1);
y(1:ndof)        = x(vidx,:); % dot(q)=v
y(ndof+1:2*ndof) = Minv*(F+G'*((G*Minv*G')\(c-G*Minv*F)));

end