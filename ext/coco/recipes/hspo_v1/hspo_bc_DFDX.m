function J = hspo_bc_DFDX(data, T, x0, x1, p)
%HSPO_BC_DFDX   Linearization of multi-segment periodic boundary conditions.
%
% Trajectory segments terminate on event surface and connect in circular
% order through resets.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: hspo_bc_DFDX.m 2839 2015-03-05 17:09:01Z fschild $

cdim  = data.cdim;  % Total state-space dimension
pdim  = data.pdim;  % Number of problem parameters
nsegs = data.nsegs; % Number of segments
vals  = [];
for i=1:nsegs
  if ~isempty(data.dhdxhan) % events w.r.t. x
    dhdx = data.dhdxhan(x1(data.x1_idx{i}), p, data.events{i});
  else
    dhdx = coco_ezDFDX('f(x,p)v', data.hhan, x1(data.x1_idx{i}), ...
      p, data.events{i});
  end
  if ~isempty(data.dhdphan) % events w.r.t. p
    dhdp = data.dhdphan(x1(data.x1_idx{i}), p, data.events{i});
  else
    dhdp = coco_ezDFDP('f(x,p)v', data.hhan, x1(data.x1_idx{i}), ...
      p, data.events{i});
  end
  if ~isempty(data.dgdxhan) % resets w.r.t. x
    dgdx = data.dgdxhan(x1(data.x1_idx{i}), p, data.resets{i});
  else
    dgdx = coco_ezDFDX('f(x,p)v', data.ghan, x1(data.x1_idx{i}), ...
      p, data.resets{i});
  end
  if ~isempty(data.dgdphan) % resets w.r.t. p
    dgdp = data.dgdphan(x1(data.x1_idx{i}), p, data.resets{i});
  else
    dgdp = coco_ezDFDP('f(x,p)v', data.ghan, x1(data.x1_idx{i}), ...
      p, data.resets{i});
  end
  vals = [vals; dhdx(:); dhdp(:); ...
    ones(data.dim(mod(i,data.nsegs)+1),1); -dgdx(:); -dgdp(:)];
end
J = sparse(data.rows, data.cols, vals, nsegs+cdim, nsegs+2*cdim+pdim);

end
