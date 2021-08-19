function Jbc = torus_bc_DFDX(data, T0, T, x0, x1, p)
%TORUS_BC_DFDX   Jacobian of torus boundary conditions.
%
% Trajectory end points lie on a curve on the invariant torus. The return
% map corresponds to identical times-of-flight and describes a rigid
% rotation.
% Appropriate phase conditions are imposed as well to yield unique solution

om2    = p(end-1);
varrho = p(end);
nt = numel(T);
nx = numel(x0);
np = numel(p);

Th  = (1:data.N)*2*pi*varrho;
SIN = [ zeros(size(Th)) ; sin(Th) ];
R   = diag([1 kron(cos(Th), [1 1])]);
R   = R  + diag(SIN(:), +1)- diag(SIN(:), -1);
RF  = kron(R*data.Fs, eye(data.dim));

% Derivative w.r.t varrho
DSIN = [ zeros(size(Th)) ; cos(Th).*(1:data.N)*2*pi ];
DR   = diag([0 kron(-sin(Th).*(1:data.N)*2*pi, [1 1])]);
DR   = DR  + diag(DSIN(:), +1)- diag(DSIN(:), -1);
DRF  = kron(DR*data.Fs, eye(data.dim));


Jbc1 = [
  eye(nt), zeros(nt,nt+2*nx+np);
  zeros(nt), eye(nt), zeros(nt,2*nx), zeros(nt,np-2),2*pi/om2^2*ones(nt,1), zeros(nt,1);
  zeros(nx,2*nt), -RF, data.F, zeros(nx,np-1), -DRF*x0;
  zeros(1,2*nt+2*nx), zeros(1,np-3), 1, -varrho, -om2];
if data.autonomous
    Jbc2 = [zeros(1,2*nt), data.f00', zeros(1,nx-data.dim+nx+np);
        zeros(1,2*nt), data.f0', zeros(1, nx-data.dim+nx+np)];
else
    Jbc2 = [zeros(1,2*nt+2*nx+np);
        zeros(1,2*nt), data.f00', zeros(1, nx-data.dim+nx+np)];
    Jbc2(1,2*nt+2*nx+data.Om2idx) = 1;
    Jbc2(1,2*nt+2*nx+np-1) = -1;
end

Jbc = [Jbc1; Jbc2];

end