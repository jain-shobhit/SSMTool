function y = pend_ode_auto(x,p)

q1 = x(1,:);
v1 = x(2,:);
c  = p(1,:);

y(1,:) = v1;
y(2,:) = -c.*v1-sin(q1);
end