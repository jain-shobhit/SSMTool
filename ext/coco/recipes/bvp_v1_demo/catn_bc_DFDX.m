function Jbc = catn_bc_DFDX(T, x0, x1, p)
%CATN_BC_DFDX   'bvp'-compatible encoding of Jacobian of catenary boundary conditions.

Jbc = zeros(3,6);
Jbc(1,1) = 1;
Jbc(2,2) = 1;
Jbc(3,4) = 1;
Jbc(3,6) = -1;

end
