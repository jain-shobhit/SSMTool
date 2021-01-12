function [data J] = calcvar_DFDU(prob, data, u)
%CALCVAR_DFDU   Linearization of zero function wrapper.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: calcvar_DFDU.m 2839 2015-03-05 17:09:01Z fschild $

x   = u(data.x_idx); % Extract basepoint values
p   = u(data.p_idx); % Extract problem parameter

f  = data.W  * x;                 % Collocation node values
fp = (2*data.NTST) * data.Wp * x; % Collocation node derivatives
pp  = repmat(p, [data.NTST*data.NCOL 1]);

d2Ldfdf   = data.fhan(f, fp, pp, 'd2Ldfdf');
d2Ldfdf   = spdiags(d2Ldfdf,0,data.NTST*data.NCOL,data.NTST*data.NCOL);
d2Ldfpdf  = (2 * data.NTST) * data.fhan(f, fp, pp, 'd2Ldfpdf');
d2Ldfpdf  = spdiags(d2Ldfpdf,0,data.NTST*data.NCOL,data.NTST*data.NCOL);
d2Ldfpdfp = (2 * data.NTST)^2 * data.fhan(f, fp, pp, 'd2Ldfpdfp');
d2Ldfpdfp = spdiags(d2Ldfpdfp,0,data.NTST*data.NCOL,data.NTST*data.NCOL);

Jint = (0.5 / data.NTST) *...
  (data.W'*data.wt*(d2Ldfdf*data.W + d2Ldfpdf*data.Wp) +...
  data.Wp'*data.wt*(d2Ldfpdf*data.W + d2Ldfpdfp*data.Wp));
Jint(data.fint1_idx,:) = Jint(data.fint1_idx,:) + Jint(data.fint2_idx,:); % Eliminate Lagrange multipliers
Jint = Jint(data.fint3_idx,:); % W.r.t. unknown basepoint values

[rows cols vals] = find(Jint);

off  = data.NTST*data.NCOL-1; % Row offset

[r c v] = find(data.Q); % Continuity conditions w.r.t. basepoint values
rows = [rows ; off + r];
cols = [cols ; c];
vals = [vals ; v];

off = off + data.NTST - 1;% Row offset

rows = [rows ; off + [1; 2]];
cols = [cols ; 1; data.NTST*(data.NCOL+1)];
vals = [vals ; 1 ; 1];

J1 = sparse(rows, cols, vals);

J2 = [zeros(data.NTST*(data.NCOL+1)-2,1); 0; -1]; % W.r.t. problem parameter

J = sparse([J1 J2]);

end
% Use the following command to check the explicit Jacobian:
% [data Jt] = coco_ezDFDX('f(o,d,x)', prob, data, @calcvar_F, u);
