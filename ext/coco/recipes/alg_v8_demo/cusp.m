function y = cusp(x, p)
%CUSP   'alg'-compatible encoding of cusp normal form

x  = x(1);
ka = p(1);
la = p(2);

y = ka-x*(la-x^2);

end
