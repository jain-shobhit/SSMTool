function fbc = po_bc(~, T, x0, x1, p) %#ok<INUSD,INUSL>
  fbc = [x1-x0; x0(3)-0.7];
end
