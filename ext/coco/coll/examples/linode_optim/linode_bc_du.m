function [data, J] = linode_bc_du(prob, data, u) %#ok<INUSD,INUSL>

J = [-1 0 1 0 0 0; 0 -1 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];

end
