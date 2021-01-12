function fbc = lin_bc(T, x0, x1, p)
%LIN_BC   'bvp'-compatible encoding of boundary conditions for forced linear oscillator.
  fbc = [x0(3); x1(1)-1; x1(2)];
end
