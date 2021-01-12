function y = linode(x, p)
%LINODE   'coll'-compatible encoding of harmonically excited linear oscillator vector field.

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
p1 = p(1,:);

y(1,:) = x2;
y(2,:) = -x2-p1.*x1+cos(x3);
y(3,:) =  1;

end
