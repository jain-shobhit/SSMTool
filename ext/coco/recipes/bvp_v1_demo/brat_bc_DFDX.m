function Jbc = brat_bc_DFDX(T, x0, x1, p)
%BRAT_BC_DFDX   'bvp'-compatible encoding of Jacobian of bratu boundary conditions.

Jbc = zeros(3,6);
Jbc(1,1) = 1;
Jbc(2,2) = 1;
Jbc(3,4) = 1;

end
