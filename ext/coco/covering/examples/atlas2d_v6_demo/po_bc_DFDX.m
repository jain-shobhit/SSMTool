function Jbc = po_bc_DFDX(T, x0, x1, p) %#ok<INUSL>
%PO_BC_DFDX   'bvp-compatible encoding of Jacobian of periodic boundary conditions on Poincare section.

Jbc = [zeros(3,1) , -eye(3) , eye(3) , zeros(3,numel(p))
       0           [0 1 0]   [0 0 0]   zeros(1,numel(p))];

end
