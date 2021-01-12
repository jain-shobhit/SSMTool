function [data y] = var_evs(prob, data, u)
%VAR_EVS   COCO-compatible encoding of eigenvector conditions.

fdata = coco_get_func_data(prob, data.tbid, 'data'); % 'varcoll' instance data

M  = reshape(u(1:end-4), fdata.u_shp);
M1 = M(fdata.M1_idx,:); % Jacobian of time-T flow

vec = u(end-3:end-1);   % Extract eigenvector
lam = u(end);           % Extract eigenvalue

y = [M1*vec-lam*vec; vec'*vec-1]; % Unit eigenvector

end
