function Jbc = po_bc_DFDX(~, T, x0, x1, p) %#ok<INUSL>
Jbc = [zeros(3,1) , -eye(3) , eye(3) , zeros(3,numel(p))
       0           [0 0 1]   [0 0 0]   zeros(1,numel(p))];
end
