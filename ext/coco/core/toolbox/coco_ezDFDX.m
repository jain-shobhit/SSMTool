function [varargout] = coco_ezDFDX(key, varargin)
%COCO_NUM_DFDX  Numerical differentiation wrt. state.
%
%   J        = COCO_EZDFDX(KEY, ...)
%   The routines COCO_NUM_DFDX,
%   COCO_NUM_DFDXV, COCO_NUM_DFDP, and COCO_NUM_DFDPV use mid-point
%   finite-difference scheme for numerically approximating the Jacobian of
%   the right-hand side with respect to the state vector at a number of
%   distinct points and with respect to the parameter vector at a number of
%   distinct points. Each function includes the possibility of
%   differentiation of an algorithm rather than a straightforward function.
%
%   NOTE: These functions are inefficient for functions F with sparse
%   Jacobian J, for example, algorithms. Use the vectorised versions
%   whenever possible for differentiating functions provided by the user.
%   In this case the variable opts.pdat.vectorised is set to 'true', which
%   is the default.
%
%   See also: gcont_num_DFDP, gcont_num_DFDPv, gcont_num_DFDXv
%

switch key
  %case 'f(x)v'
  %case 'f(x)'
  case 'f(x,p)v'
    [varargout{1:nargout}] = coco_num_DFDXv__A(varargin{:});
  case 'f(x,p)'
    [varargout{1:nargout}] = coco_num_DFDX__A(varargin{:});
  case 'f(o,x)v'
    [varargout{1:nargout}] = coco_num_DFDXv__B(varargin{:});
  case 'f(o,x)'
    [varargout{1:nargout}] = coco_num_DFDX__B(varargin{:});
  case 'f(o,d,x)v'
    [varargout{1:nargout}] = coco_num_DFDXv__C(varargin{:});
  case 'f(o,d,x)'
    [varargout{1:nargout}] = coco_num_DFDX__C(varargin{:});
  case 'f(o,x,p)v'
    [varargout{1:nargout}] = coco_num_DFDXv__D(varargin{:});
  case 'f(o,x,p)'
    [varargout{1:nargout}] = coco_num_DFDX__D(varargin{:});
  case 'f(o,d,x,p)v'
    [varargout{1:nargout}] = coco_num_DFDXv__E(varargin{:});
  case 'f(o,d,x,p)'
    [varargout{1:nargout}] = coco_num_DFDX__E(varargin{:});
  case 'f(o,d,u,l,v)v'
    [varargout{1:nargout}] = coco_num_DFDXv__F(varargin{:});
  case 'f(o,d,u,l,v)'
    [varargout{1:nargout}] = coco_num_DFDX__F(varargin{:});
  otherwise
    error('%s: key ''%s'' unrecognised', mfilename, key);
end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function J = coco_num_DFDXv__A(F, x, p, varargin)

[m, n] = size(x);
idx    = repmat(1:n, [m 1]);
x0     = x(:,idx);
p0     = p(:,idx);

idx = repmat(1:m, [1 n]);
idx = sub2ind([m m*n], idx, 1:m*n);

h = 1.0e-8*(1.0 + abs(x0(idx)));
x = x0;

x(idx) = x0(idx)+h;
fr     = F(x, p0, varargin{:});
x(idx) = x0(idx)-h;
fl     = F(x, p0, varargin{:});

l  = size(fr);
hi = reshape(repmat(0.5./h, [prod(l(1:end-1)),1]), l);
J  = reshape(hi.*(fr-fl), [l(1:end-1) m n]);

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function J = coco_num_DFDX__A(F, x, p, varargin)

[m, n] = size(x);

fr = F(x(:,1), p(:,1), varargin{:});
l  = size(fr);

J  = zeros(prod(l),m,n);

for j=1:n
	x0 = x(:,j);
	h  = 1.0e-8*( 1.0 + abs(x0) );
	hi = 0.5./h;
	for i=1:m
		xx       = x0;
		xx(i)    = x0(i)+h(i);
		fr       = F(xx, p(:,j), varargin{:});
		xx(i)    = x0(i)-h(i);
		fl       = F(xx, p(:,j), varargin{:});
		J(:,i,j) = hi(i)*(fr(:)-fl(:));
	end
end

if numel(l)==2 && l(2)==1
  J = reshape(J, [l(1) m n]);
else
  J = reshape(J, [l m n]);
end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [opts J] = coco_num_DFDXv__B(opts, F, x, varargin)

[m n] = size(x);
h     = reshape(1.0e-8*( 1.0 + abs(x) ), 1, m*n);

idx1 = kron( 1:n , ones(1,m) );
x0   = x(:,idx1);

idx1 = repmat(1:m, [1 n]);
idx2 = 1:m*n;
idx  = sub2ind([m m*n], idx1, idx2);

x    = x0;

x(idx)    = x0(idx)+h;
[opts fr] = F(opts, x,varargin{:});
x(idx)    = x0(idx)-h;
[opts fl] = F(opts, x,varargin{:});

l  = size(fr, 1);
hi = repmat(0.5./h, l, 1);
J  = hi.*(fr-fl);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [opts J] = coco_num_DFDX__B(opts, F, x, varargin)

[m n] = size(x);

[opts fr] = F(opts, x(:,1),varargin{:});
l  = size(fr,1);

J  = zeros(l,m,n);

for j=1:n
	x0 = x(:,j);
	h  = 1.0e-8*( 1.0 + abs(x0) );
	hi = 0.5./h;
	for i=1:m
		xx        = x0;
		xx(i)     = x0(i)+h(i);
		[opts fr] = F(opts, xx,varargin{:});
		xx(i)     = x0(i)-h(i);
		[opts fl] = F(opts, xx,varargin{:});
		J(:,i,j)  = hi(i)*(fr-fl);
	end
end
J = reshape(J, [l m*n]);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [data, J] = coco_num_DFDXv__C(opts, data, F, x, varargin)

[m, n] = size(x);
idx   = repmat(1:n, [m 1]);
x0    = x(:,idx);

idx = repmat(1:m, [1 n]);
idx = sub2ind([m m*n], idx, 1:m*n);

h = 1.0e-8*(1.0 + abs(x0(idx)));
x = x0;

x(idx)     = x0(idx)+h;
[data, fr] = F(opts, data, x, varargin{:});
x(idx)     = x0(idx)-h;
[data, fl] = F(opts, data, x, varargin{:});

l  = size(fr);
hi = reshape(repmat(0.5./h, [prod(l(1:end-1)),1]), l);
J  = reshape(hi.*(fr-fl), [l(1:end-1) m n]);

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [data, J] = coco_num_DFDX__C(opts, data, F, x, varargin)

[m, n] = size(x);

[data, fr] = F(opts, data, x(:,1),varargin{:});
l  = size(fr);

J  = zeros(prod(l),m,n);

for j=1:n
	x0 = x(:,j);
	h  = 1.0e-8*( 1.0 + abs(x0) );
	hi = 0.5./h;
	for i=1:m
		xx         = x0;
		xx(i)      = x0(i)+h(i);
		[data, fr] = F(opts, data, xx, varargin{:});
		xx(i)      = x0(i)-h(i);
		[data, fl] = F(opts, data, xx, varargin{:});
    J(:,i,j)   = hi(i)*(fr(:)-fl(:));
	end
end

if numel(l)==2 && l(2)==1
  J = reshape(J, [l(1) m n]);
else
  J = reshape(J, [l m n]);
end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [opts J] = coco_num_DFDXv__D(opts, F, x, p, varargin)

[m n] = size(x);
h     = reshape(1.0e-8*( 1.0 + abs(x) ), 1, m*n);

idx1 = kron( 1:n , ones(1,m) );
x0   = x(:,idx1);

idx1 = repmat(1:m, [1 n]);
idx2 = 1:m*n;
idx  = sub2ind([m m*n], idx1, idx2);

x    = x0;

x(idx)    = x0(idx)+h;
[opts fr] = F(opts, x, p,varargin{:});
x(idx)    = x0(idx)-h;
[opts fl] = F(opts, x, p,varargin{:});

l  = size(fr, 1);
hi = repmat(0.5./h, l, 1);
J  = hi.*(fr-fl);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [opts J] = coco_num_DFDX__D(opts, F, x, p, varargin)

[m n] = size(x);

[opts fr] = F(opts, x(:,1),varargin{:});
l  = size(fr,1);

J  = zeros(l,m,n);

for j=1:n
	x0 = x(:,j);
	h  = 1.0e-8*( 1.0 + abs(x0) );
	hi = 0.5./h;
	for i=1:m
		xx        = x0;
		xx(i)     = x0(i)+h(i);
		[opts fr] = F(opts, xx, p,varargin{:});
		xx(i)     = x0(i)-h(i);
		[opts fl] = F(opts, xx, p,varargin{:});
		J(:,i,j)  = hi(i)*(fr-fl);
	end
end
J = reshape(J, [l m*n]);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [data J] = coco_num_DFDXv__E(opts, data, F, x, p, varargin)

[m n] = size(x);
h     = reshape(1.0e-8*( 1.0 + abs(x) ), 1, m*n);

idx1 = kron( 1:n , ones(1,m) );
x0   = x(:,idx1);

idx1 = repmat(1:m, [1 n]);
idx2 = 1:m*n;
idx  = sub2ind([m m*n], idx1, idx2);

x    = x0;

x(idx)    = x0(idx)+h;
[data fr] = F(opts, data, x, p,varargin{:});
x(idx)    = x0(idx)-h;
[data fl] = F(opts, data, x, p,varargin{:});

l  = size(fr, 1);
hi = repmat(0.5./h, l, 1);
J  = hi.*(fr-fl);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [data J] = coco_num_DFDX__E(opts, data, F, x, p, varargin)

[m n] = size(x);

[data fr] = F(opts, data, x(:,1),varargin{:});
l  = size(fr,1);

J  = zeros(l,m,n);

for j=1:n
	x0 = x(:,j);
	h  = 1.0e-8*( 1.0 + abs(x0) );
	hi = 0.5./h;
	for i=1:m
		xx        = x0;
		xx(i)     = x0(i)+h(i);
		[data fr] = F(opts, data, xx, p,varargin{:});
		xx(i)     = x0(i)-h(i);
		[data fl] = F(opts, data, xx, p,varargin{:});
		J(:,i,j)  = hi(i)*(fr-fl);
	end
end
J = reshape(J, [l m*n]);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [data, J] = coco_num_DFDXv__F(opts, data, F, u, l, v, varargin)

% w.r.t. u
if isempty(u)
  Ju = [];
else
  [m, n] = size(u);
  idx   = repmat(1:n, [m 1]);
  u0    = u(:,idx);
  
  idx = repmat(1:m, [1 n]);
  idx = sub2ind([m m*n], idx, 1:m*n);
  
  h = 1.0e-8*(1.0 + abs(u0(idx)));
  u = u0;
  
  u(idx)     = u0(idx)+h;
  [data, fr] = F(opts, data, u, l, v, varargin{:});
  u(idx)     = u0(idx)-h;
  [data, fl] = F(opts, data, u, l, v, varargin{:});
  
  k  = size(fr);
  hi = reshape(repmat(0.5./h, [prod(k(1:end-1)),1]), k);
  Ju  = reshape(hi.*(fr-fl), [k(1:end-1) m n]);
end

% w.r.t. l
if isempty(l)
  Jl = [];
else
  [m, n] = size(l);
  idx   = repmat(1:n, [m 1]);
  l0    = l(:,idx);
  
  idx = repmat(1:m, [1 n]);
  idx = sub2ind([m m*n], idx, 1:m*n);
  
  h = 1.0e-8*(1.0 + abs(l0(idx)));
  l = l0;
  
  l(idx)     = l0(idx)+h;
  [data, fr] = F(opts, data, u, l, v, varargin{:});
  l(idx)     = l0(idx)-h;
  [data, fl] = F(opts, data, u, l, v, varargin{:});
  
  k  = size(fr);
  hi = reshape(repmat(0.5./h, [prod(k(1:end-1)),1]), k);
  Jl  = reshape(hi.*(fr-fl), [k(1:end-1) m n]);
end

% w.r.t. v
if isempty(v)
  Jv = [];
else
  [m, n] = size(v);
  idx   = repmat(1:n, [m 1]);
  v0    = v(:,idx);
  
  idx = repmat(1:m, [1 n]);
  idx = sub2ind([m m*n], idx, 1:m*n);
  
  h = 1.0e-8*(1.0 + abs(v0(idx)));
  v = v0;
  
  v(idx)     = v0(idx)+h;
  [data, fr] = F(opts, data, u, l, v, varargin{:});
  v(idx)     = v0(idx)-h;
  [data, fl] = F(opts, data, u, l, v, varargin{:});
  
  k  = size(fr);
  hi = reshape(repmat(0.5./h, [prod(k(1:end-1)),1]), k);
  Jv  = reshape(hi.*(fr-fl), [k(1:end-1) m n]);
end

J = {Ju, Jl, Jv};

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [data, J] = coco_num_DFDX__F(opts, data, F, u, l, v, varargin)

% w.r.t. u
if isempty(u)
  Ju = [];
else
  [m, n] = size(u);
  
  [data, fr] = F(opts, data, u(:,1), l, v, varargin{:});
  k  = size(fr);
  
  Ju  = zeros(prod(k),m,n);
  
  for j=1:n
    u0 = u(:,j);
    h  = 1.0e-8*( 1.0 + abs(u0) );
    hi = 0.5./h;
    for i=1:m
      uu         = u0;
      uu(i)      = u0(i)+h(i);
      [data, fr] = F(opts, data, uu, l, v, varargin{:});
      uu(i)      = u0(i)-h(i);
      [data, fl] = F(opts, data, uu, l, v, varargin{:});
      Ju(:,i,j)   = hi(i)*(fr(:)-fl(:));
    end
  end
  
  if numel(k)==2 && k(2)==1
    Ju = reshape(Ju, [k(1) m n]);
  else
    Ju = reshape(Ju, [k m n]);
  end
end
% w.r.t. l
if isempty(l)
  Jl = [];
else
  [m, n] = size(l);
  
  [data, fr] = F(opts, data, u, l(:,1), v, varargin{:});
  k  = size(fr);
  
  Jl  = zeros(prod(k),m,n);
  
  for j=1:n
    l0 = l(:,j);
    h  = 1.0e-8*( 1.0 + abs(l0) );
    hi = 0.5./h;
    for i=1:m
      ll         = l0;
      ll(i)      = l0(i)+h(i);
      [data, fr] = F(opts, data, u, ll, v, varargin{:});
      ll(i)      = l0(i)-h(i);
      [data, fl] = F(opts, data, u, ll, v, varargin{:});
      Jl(:,i,j)   = hi(i)*(fr(:)-fl(:));
    end
  end
  
  if numel(k)==2 && k(2)==1
    Jl = reshape(Jl, [k(1) m n]);
  else
    Jl = reshape(Jl, [k m n]);
  end
end

% w.r.t. v
if isempty(v)
  Jv = [];
else
  [m, n] = size(v);
  
  [data, fr] = F(opts, data, u, l, v(:,1),varargin{:});
  k  = size(fr);
  
  Jv  = zeros(prod(k),m,n);
  
  for j=1:n
    v0 = v(:,j);
    h  = 1.0e-8*( 1.0 + abs(v0) );
    hi = 0.5./h;
    for i=1:m
      vv         = v0;
      vv(i)      = v0(i)+h(i);
      [data, fr] = F(opts, data, u, l, vv, varargin{:});
      vv(i)      = v0(i)-h(i);
      [data, fl] = F(opts, data, u, l, vv, varargin{:});
      Jv(:,i,j)   = hi(i)*(fr(:)-fl(:));
    end
  end
  
  if numel(k)==2 && k(2)==1
    Jv = reshape(Jv, [k(1) m n]);
  else
    Jv = reshape(Jv, [k m n]);
  end
end

J = {Ju, Jl, Jv};

end
