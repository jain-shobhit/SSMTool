function fbc = caty_bc(T, x0, x1, p)
%CATY_BC   'bvp'-compatible encoding of catenary boundary conditions.
  fbc = [T-1; x0(1)-1; x1(1)-p];
end
