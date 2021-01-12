function f = duff(t,x,p, cdata)
% p = [A om la al eps]
% x'' + la*x' + al*x + eps*x^3 = A*cos(om*t)

y      = x(1,:);
yr     = cdata.yr(t);
P      = y - yr;
D      = (P-x(3,:))/cdata.scale;
u      = cdata.k1*P + cdata.k2*D;

f(1,:) = x(2,:);
f(2,:) = p(1)*cos(p(2)*t) - p(3)*x(2,:) - p(4)*x(1,:) - p(5)*x(1,:).^3 + u;
f(3,:) = D;

end
