function y = vanderpol(x, p)
%VANDERPOL   'coll'-compatible encoding of Van der Pol vector field

y(1,:) = p(2,:).*(p(1,:)-x(2,:));
y(2,:) = x(1,:)+x(2,:)-x(2,:).^3/3;

end
