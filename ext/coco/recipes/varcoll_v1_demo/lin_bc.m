function fbc = lin_bc(T, x0, x1, p)
%LIN_BC   'bvp'-compatible encoding of boundary conditions for forced linear oscillator.
  fbc = [x1(1:2)-x0(1:2); x1(3)-x0(3)-2*pi; x0(1)];
end
