function [prob, data] = ep_construct_eqn(prob, data, sol, CacheJacobian)
%EP_CONSTRUCT_EQN   Add EP zero problem.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_construct_eqn.m 2872 2015-08-06 20:12:06Z hdankowicz $

data = init_data(data);
eqn  = data.ep_eqn;
fid  = eqn.fid;

if CacheJacobian
  data.sh.ep_Jx = [];
  data.sh.ep_Jp = [];
  prob = coco_add_func(prob, fid, @FDF1, data, 'zero', ...
    'u0', sol.u0, 't0', sol.t0, 'F+DF');
  data = data.protect('ep_Jx', 'ep_Jp');
  data.no_save = [ data.no_save { 'ep_Jx' 'ep_Jp' } ];
else
  prob = coco_add_func(prob, fid, @FDF2, data, 'zero', ...
    'u0', sol.u0, 't0', sol.t0, 'F+DF');
end

if ~isempty(data.pnames)
  uidx = coco_get_func_data(prob, fid, 'uidx');
  pfid = coco_get_id(fid, 'pars');
  prob = coco_add_pars(prob, pfid, uidx(eqn.p_idx), data.pnames);
end

end

function data = init_data(data)
%INIT_DATA   Initialize data for EP zero problem.

xdim      = data.xdim;
pdim      = data.pdim;
eqn.x_idx = (1:xdim)';
eqn.p_idx = xdim+(1:pdim)';
eqn.fid   = coco_get_id(data.oid, 'ep');

data.ep_eqn  = eqn;

end

function [data, y, J] = FDF1(~, data, u)
%FDF1   EP zero problem with cacheing of derivatives.

pr  = data.pr;
eqn = pr.ep_eqn;

x = u(eqn.x_idx);
p = u(eqn.p_idx);

y  = pr.ode_F(pr, [], x, p);
Jx = pr.ode_DFDX(pr, [], x, p);
Jp = pr.ode_DFDP(pr, [], x, p);

J = [Jx Jp];

data.ep_Jx = Jx;
data.ep_Jp = Jp;
end

function [data, y, J] = FDF2(~, data, u)
%FDF2   EP zero problem without cacheing of derivatives.

pr  = data.pr;
eqn = pr.ep_eqn;

x = u(eqn.x_idx);
p = u(eqn.p_idx);

y = pr.ode_F(pr, [], x, p);

if nargout>=3
  Jx = pr.ode_DFDX(pr, [], x, p);
  Jp = pr.ode_DFDP(pr, [], x, p);
  J  = [Jx Jp];
end

end
