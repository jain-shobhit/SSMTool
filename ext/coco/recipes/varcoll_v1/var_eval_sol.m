function M = var_eval_sol(data, u)
%VAR_EVAL_SOL   Evaluate Jacobian of time-T flow
%
% Solve the variational equations with initial condition given by the
% identity matrix and extract fundamental solution at t=1.

x = u(data.xbp_idx); % Extract basepoint values
T = u(data.T_idx);   % Extract interval length
p = u(data.p_idx);   % Extract problem parameters

xx = reshape(data.W*x, data.x_shp); % Values at collocation nodes
pp = repmat(p, data.p_rep);

if isempty(data.dfdxhan)
  dxode = coco_ezDFDX('f(x,p)v', data.fhan, xx, pp);
else
  dxode = data.dfdxhan(xx, pp);
end
dxode = sparse(data.dxrows, data.dxcols, dxode(:));
dxode = (0.5*T/data.coll.NTST)*dxode*data.W-data.Wp;    % Variational equation

dim    = data.dim;                                      % State-space dimension
M1_idx = data.xbp_idx(end-dim+(1:dim));                 % Index of Jacobian of time-T flow
row    = [eye(dim), zeros(dim, data.xbp_idx(end)-dim)]; % Initial condition and rhs

M = [row; dxode; data.Q]\row'; % Fundamental solution
M = M(M1_idx,:); 

end
