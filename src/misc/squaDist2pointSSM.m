function y = squaDist2pointSSM(z0,u,W_0,data)
% SQUADIST2POINTSSM This function compute the squared distance between a point z0 and a point u on
% SSM. z0 here is in full physicial space. u is a point gives coordinates
% of this point on SSM charaterized by W_0. In addition, u here is
% rearranged as real coordinates. Such a mapping is given by data


% mapping u to (complex conjugate) state
realx = data.realx;
compx = data.compx;
realz = realx;
compz = compx;

x_real = u(realz);
x_comp = u(compz(1:2:end-1))+1i*u(compz(2:2:end));

x  = zeros(data.dim,1); % state
x(realx) = x_real;
x(compx(1:2:end-1)) = x_comp;
x(compx(2:2:end))   = conj(x_comp);

% locating corresponding z in physical space
z = reduced_to_full(x,W_0);

% calculating squared distance between z and z0
y = sum((z-z0).^2);

end