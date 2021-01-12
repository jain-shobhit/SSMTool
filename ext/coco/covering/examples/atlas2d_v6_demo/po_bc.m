function fbc = po_bc(T, x0, x1, p) %#ok<INUSD,INUSL>
%PO_BC   'bvp-compatible encoding of periodic boundary conditions on Poincare section.
  fbc = [x1-x0; x0(2)];
end
