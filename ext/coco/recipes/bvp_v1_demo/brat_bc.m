function fbc = brat_bc(T, x0, x1, p)
%BRAT_BC   'bvp'-compatible encoding of bratu boundary conditions.
  fbc = [T-1; x0(1); x1(1)];
end
