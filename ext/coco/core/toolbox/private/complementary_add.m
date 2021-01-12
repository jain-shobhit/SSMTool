function prob = complementary_add(prob, data, args)
%COMPLEMENTARY_ADD Add complementary zero/monitor function wrapper to continuation problem

data = init_comp_data(prob, data);
uidx = [ data.u_idx data.lu_idx data.vu_idx(1:end-data.v_num) ];
prob = coco_add_func(prob, data.identifyer, @coco_complementary, data, ...
  args{:}, 'uidx', uidx, 'f+df', 'remesh', @remesh_comp);
prob = coco_add_slot(prob, data.identifyer, @coco_save_data, ...
  data, 'save_full');

end

function data = init_comp_data(prob, data)

data.uidx = [];
data.lidx = [];
data.vidx = [];
adata = coco_get_func_data(prob, 'coco_adjoint', 'data');

if ~isempty(data.u_idx)
  data.uidx = 1:numel(data.u_idx);
end
if ~isempty(data.l_idx)
  data.lidx = numel(data.u_idx)+(1:numel(data.l_idx));
end
if ~isempty(data.v_idx)
  data.vidx = numel(data.u_idx)+numel(data.l_idx)+(1:numel(data.v_idx));
end
data.lu_idx = adata.x_dim + data.l_idx;
data.vu_idx = adata.x_dim + adata.l_dim + data.v_idx;

data.axdim = adata.x_dim;
data.aldim = adata.l_dim;

end

function [data, y, J] = coco_complementary(prob, data, u)

u0 = u(data.uidx);
l0 = u(data.lidx);
v0 = u(data.vidx);

switch data.call_mode
  case 1
    [data.data, y] = complementary_call_F(prob, data.data, data, u0, l0, v0);
    if nargout==3
      [data.data, J] = complementary_call_DF(prob, data.data, data, u0, l0, v0);
      J = cell2mat(J);
    end
  case 2
    [data.data, y, J] = complementary_call_FDF(prob, data.data, data, u0, l0, v0);
    J = cell2mat(J);
end

end

function [prob, stat, xtr] = remesh_comp(prob, data, chart, ub, Vb)

efunc         = prob.efunc;
complementary = prob.complementary;

prob.complementary.cfidx = ...
  find(strcmpi(data.identifyer, complementary.identifyers),1);

x0beg = numel(efunc.x0);
xtr   = efunc.xtr;

if isempty(data.remesh)
  adata = coco_get_func_data(prob, 'coco_adjoint', 'data');

  uidx = xtr(data.u_idx(data.u_idx<=numel(xtr)));
  lidx = xtr(data.lu_idx(data.lu_idx<=numel(xtr)))-adata.x_dim;
  vidx = xtr(data.vu_idx(data.vu_idx<=numel(xtr)))-adata.x_dim-adata.l_dim;
  if any(uidx<=0) || any(lidx<=0) || any(vidx<=0)
    error('%s: compulsory variables have been removed, cannot relocate ''%s''', ...
      mfilename, data.identifyer);
  end
  v0idx = data.vidx(data.vu_idx>numel(xtr));
  prob  = coco_change_comp(prob, data, 'uidx', uidx, 'lidx' , lidx, ...
    'vidx', vidx, 'v0', ub(v0idx), 'vecs', Vb(v0idx,:));
  xtr  = [ xtr ; x0beg+(1:numel(v0idx))' ];
  stat = 'success';
else
  v0idx = data.vidx;
  [prob, stat, tr] = data.remesh(prob, data, chart, ub(v0idx), ...
    Vb(v0idx,:));
  xbeg = x0beg*(tr~=0);
  xtr  = [ xtr ; xbeg(:)+tr(:) ];
end

data = prob.complementary.funcs(prob.complementary.cfidx);
data = init_comp_data(prob, data);
uidx = [data.u_idx data.lu_idx data.vu_idx(1:end-data.v_num)];
prob = coco_change_func(prob, data, 'uidx', uidx, ...
  'u0', ub(end-data.v_num+1:end), 'vecs', Vb(end-data.v_num+1:end,:));

end
