function [data, y] = map(prob, data, u)

x0 = u(1);
y0 = u(2);
x1 = u(3);
y1 = u(4);
om = u(5);
fr = u(6);
T  = u(7);

ex = exp(T);
co = cos(fr*T);
si = sin(fr*T);
r0 = sqrt(x0^2+y0^2);
de = 1+om^2-r0*om^2+ex*r0*(1+om^2-cos(om*T)-om*sin(om*T));

y = [x1 - ex*(1+om^2)*(x0*co-y0*si)/de; y1 - ex*(1+om^2)*(y0*co+x0*si)/de];

end
