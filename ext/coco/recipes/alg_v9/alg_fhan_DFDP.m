function Jp = alg_fhan_DFDP(data, x, p)
%ALG_FHAN_DFDP   Compute p-derivative of DATA.FHAN(X,P).
%
% When dfdphan is empty, an approximate Jacobian is obtained using
% numerical differentiation.
%
% Identical to alg_v7.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_fhan_DFDP.m 2839 2015-03-05 17:09:01Z fschild $

if isempty(data.dfdphan)
  Jp = coco_ezDFDP('f(x,p)', data.fhan, x, p);
else
  Jp = data.dfdphan(x, p);
end

end
