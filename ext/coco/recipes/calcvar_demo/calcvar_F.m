function [data y] = calcvar_F(prob, data, u)
%CALCVAR_F   Zero function wrapper to definition of functional.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: calcvar_F.m 2839 2015-03-05 17:09:01Z fschild $

x   = u(data.x_idx); % Extract basepoint values
p   = u(data.p_idx); % Extract problem parameter

f  = data.W  * x;                 % Collocation node values
fp = (2*data.NTST) * data.Wp * x; % Collocation node derivatives
pp = repmat(p, [data.NTST*data.NCOL 1]);

dLdf  = data.fhan(f, fp, pp, 'dLdf'); % Partial w.r.t. f
dLdfp = (2 * data.NTST) * data.fhan(f, fp, pp, 'dLdfp'); % Partial w.r.t. f'

fint  = (0.5 / data.NTST) * (data.W' * data.wt * dLdf...
  + data.Wp' * data.wt * dLdfp); % Jacobian of quadrature approximation
fint(data.fint1_idx) = fint(data.fint1_idx) + fint(data.fint2_idx); % Eliminate Lagrange multipliers
fint  = fint(data.fint3_idx);                        % Quadrature conditions
fcont = data.Q * x;                                  % Continuity conditions
fbound = [x(1) - 1; x(data.NTST*(data.NCOL+1)) - p]; % Boundary conditions

y = [fint; fcont; fbound];

end
