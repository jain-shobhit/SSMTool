function y = chemosc(x, p)
%CHEMOSC   'dft'-compatible encoding of vector field.

p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
p4 = p(4,:);
p5 = p(5,:);
p6 = p(6,:);
p7 = p(7,:);
p8 = p(8,:);
p9 = p(9,:);

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
x4 = x(4,:);

T1 = p1.*x1.*x2.*x3;
T2 = p3.*x1.*x2.*x4;
T3 = 2*p2.*x3.*x3;

y(1,:) = -T1-T2+p7-p9.*x1;
y(2,:) = -T1-T2+p8;
y(3,:) =  T1-T3+2*T2-p4.*x3+p6;
y(4,:) = -T2+T3-p5.*x4;

end
