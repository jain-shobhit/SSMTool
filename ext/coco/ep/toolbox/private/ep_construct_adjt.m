function [prob, data] = ep_construct_adjt(prob, data, sol)
%EP_CONSTRUCT_ADJT   Add EP adjoint problem.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ep_construct_adjt.m 2872 2015-08-06 20:12:06Z hdankowicz $

data = init_data(data);
opt  = data.ep_opt;
fid  = opt.fid;

prob   = coco_add_adjt(prob, fid, @adj, @adj_DU, data, 'l0', sol.l0, ...
  'tl0', sol.tl0);
if ~isempty(data.pnames)
  pfid   = coco_get_id(fid, 'pars');
  dnames = coco_get_id('d', data.pnames);
  axidx  = coco_get_adjt_data(prob, fid, 'axidx');
  prob   = coco_add_adjt(prob, pfid, dnames, 'aidx', axidx(opt.p_idx), ...
    'l0', sol.pars_l0, 'tl0', sol.pars_tl0);
end

end

function data = init_data(data)
%INIT_DATA   Initialize data for EP adjoint problem.

xdim      = data.xdim;
pdim      = data.pdim;
opt.x_idx = (1:xdim)';
opt.p_idx = xdim+(1:pdim)';

opt.dJrows = xdim;
opt.dJcols = (xdim+pdim)^2;
opt.dfdxdxrows = repmat(1:xdim, [1, xdim^2]);
cols = repmat(1:xdim, [xdim, 1]);
opt.dfdxdxcols = repmat(cols(:), [1, xdim]) + ...
  repmat(0:xdim+pdim:(xdim-1)*(xdim+pdim), [xdim^2, 1]);

opt.dfdxdprows = repmat(1:xdim, [1, xdim*pdim]);
cols = repmat(1:xdim, [xdim, 1]);
opt.dfdxdpcols = repmat(cols(:), [1, pdim]) + ...
  repmat(xdim*(xdim+pdim):xdim+pdim:(xdim+pdim-1)*(xdim+pdim), ...
  [xdim^2, 1]);

opt.dfdpdxrows = repmat(1:xdim, [1, pdim*xdim]);
cols = repmat(1:xdim+pdim:xdim*(xdim+pdim), [xdim, 1]);
opt.dfdpdxcols = repmat(cols(:), [1, pdim]) + ...
  repmat(xdim:xdim+pdim-1, [xdim^2, 1]);

opt.dfdpdprows = repmat(1:xdim, [1, pdim*pdim]);
cols = repmat(1:pdim, [xdim, 1]);
opt.dfdpdpcols = repmat(cols(:), [1, pdim]) + ...
  repmat(xdim*(xdim+pdim+1):xdim+pdim:(xdim+pdim-1)*(xdim+pdim)+xdim, ...
  [xdim*pdim, 1]);

opt.fid   = coco_get_id(data.oid, 'ep');

data.ep_opt  = opt;

end

function [data, J] = adj(prob, data, u) %#ok<INUSL>

pr  = data.pr;
opt = pr.ep_opt;

x = u(opt.x_idx);
p = u(opt.p_idx);

Jx = pr.ode_DFDX(pr, [], x, p);
Jp = pr.ode_DFDP(pr, [], x, p);

J = [Jx Jp];

end

function [data, dJ] = adj_DU(prob, data, u) %#ok<INUSL>

pr  = data.pr;
opt = pr.ep_opt;

x = u(opt.x_idx);
p = u(opt.p_idx);

dfdxdx = pr.ode_DFDXDX(pr, [], x, p);
dfdxdp = pr.ode_DFDXDP(pr, [], x, p);
dfdpdp = pr.ode_DFDPDP(pr, [], x, p);

dJrows = opt.dJrows;
dJcols = opt.dJcols;

dJ = ...
  sparse(opt.dfdxdxrows, opt.dfdxdxcols, dfdxdx(:), dJrows, dJcols) + ...
  sparse(opt.dfdxdprows, opt.dfdxdpcols, dfdxdp(:), dJrows, dJcols) + ...
  sparse(opt.dfdpdxrows, opt.dfdpdxcols, dfdxdp(:), dJrows, dJcols) + ...
  sparse(opt.dfdpdprows, opt.dfdpdpcols, dfdpdp(:), dJrows, dJcols);
  
end
