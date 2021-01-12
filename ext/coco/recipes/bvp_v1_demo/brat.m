function y = brat(x, p)
%BRAT   'coll'-compatible encoding of bratu vector field.

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);

y(1,:) = x2;
y(2,:) = -p1.*exp(x1);

end
