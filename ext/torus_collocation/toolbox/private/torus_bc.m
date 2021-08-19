function fbc = torus_bc(data, T0, T, x0, x1, p)
%TORUS_BC   Torus boundary conditions.
%
% Trajectory end points lie on a curve on the invariant torus. The return
% map corresponds to identical times-of-flight and describes a rigid
% rotation.
% Appropriate phase conditions are imposed as well to yield unique solution

om1    = p(end-2);
om2    = p(end-1);
varrho = p(end);
Th  = (1:data.N)*2*pi*varrho;
SIN = [ zeros(size(Th)) ; sin(Th) ];
R   = diag([1 kron(cos(Th), [1 1])]);
R   = R  + diag(SIN(:), +1)- diag(SIN(:), -1);
RF  = kron(R*data.Fs, eye(data.dim));

fbc1 = [T0; T-2*pi/om2; data.F*x1-RF*x0; om1-varrho*om2];

if data.autonomous
    phase1 = data.f00'*(x0(1:data.dim)-data.x00);
    phase2 = data.f0'*(x0(1:data.dim)-data.x0);
    fbc2   = [phase1; phase2];
else
    Om2 = p(data.Om2idx);
    par_coup = Om2-om2;
    % Poincare condition: <v^ast_\phi, v(0,0)-v^\ast>=0 
    % Here v(0,0) is the evaluation of the first segment at t=0. v^\ast
    % is the same evaluation at the solution of previous continuation
    % step.
    phase = data.f00'*(x0(1:data.dim)-data.x00);
    fbc2  = [par_coup; phase];
end

fbc = [fbc1; fbc2];

end
