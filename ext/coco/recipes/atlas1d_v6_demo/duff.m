function y = duff(x,p)
%DUFF   'coll'-compatible encoding of the duffing vector field.

la = p(1,:);
al = p(2,:);
ep = p(3,:);
A  = p(4,:);
om = p(5,:);

y(1,:) = x(2,:);
y(2,:) = A.*x(4,:)-la.*x(2,:)-al.*x(1,:)-ep.*x(1,:).^3;

ss     = x(3,:).*x(3,:)+x(4,:).*x(4,:);
y(3,:) =      x(3,:) + om.*x(4,:) - x(3,:).*ss;
y(4,:) = -om.*x(3,:) +     x(4,:) - x(4,:).*ss;

end
