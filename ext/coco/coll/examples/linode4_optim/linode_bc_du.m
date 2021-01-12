function J = linode_bc_du(~,T,x0,x1,p) %#ok<INUSD>
J = [0 -1 0 0 1 0 0 0; 0 0 -1 0 0 1 0 0; 0 0 1 0 0 0 0 0; 0 0 0 -1 0 0 1 0];
end
