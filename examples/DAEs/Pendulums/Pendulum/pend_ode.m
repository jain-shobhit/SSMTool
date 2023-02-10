function y = pend_ode(t,x,p)

c = 0.1;
q1 = x(1,:);
v1 = x(2,:);
epf = p(1,:);
om  = p(2,:);

y(1,:) = v1;
y(2,:) = epf.*cos(om.*t)-c*v1-sin(q1);
end
