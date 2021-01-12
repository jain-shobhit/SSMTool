function Jt = linode_DFDT(t, x, p)
%LINODE_DFDT   'coll'-compatible encoding of Jacobian of 'linode' vector field w.r.t. time.
%
% Encoding is of a non-autonomous vector field.

Jt = zeros(2,numel(t));
Jt(2,:) = -sin(t);

end
