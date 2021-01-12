function y = vdp(t, x, p)
%VDP   'coll'-compatible encoding of langford vector field.

x1 = x(1,:);
x2 = x(2,:);
om = p(1,:);
c  = p(2,:);
a  = p(3,:);

y(1,:) = x2;
y(2,:) = c.*x2.*(1-x1.^2)-x1+a.*cos(om.*t);

end
