function Jx = alg_fhan_DFDX(data, x, p)
%ALG_FHAN_DFDX   Compute x-derivative of DATA.FHAN(X,P).
%
% When dfdxhan is empty, an approximate Jacobian is obtained using
% numerical differentiation.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_fhan_DFDX.m 2839 2015-03-05 17:09:01Z fschild $

if isempty(data.dfdxhan)
  Jx = coco_ezDFDX('f(x,p)', data.fhan, x, p);
else
  Jx = data.dfdxhan(x, p);
end

end
