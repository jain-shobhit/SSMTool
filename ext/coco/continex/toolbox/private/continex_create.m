function prob = continex_create(prob, oid, data, sol)

tbid = coco_get_id(oid, 'continex');

% use continex specific atlas and corrector classes
if ~coco_exist('atlas', 'class_prop', prob, 'cont', '-no-inherit')
	prob = coco_set(prob, 'cont', 'atlas', 'continex');
end
if ~coco_exist('corrector', 'class_prop', prob, 'cont', '-no-inherit')
  prob = coco_set(prob, 'cont', 'corrector', 'continex');
end

% get toolbox settings
data = continex_get_settings(prob, tbid, data);

% add zero problem
data.scaleF  = data.continex.scale;
data.scale   = data.scaleF(sol.x0, sol.p0);
data.jac     = true;
data.NJac    = data.continex.NJac;
data.DFCount = data.continex.NJac;
data.JacH    = data.continex.JacH;
data.oid     = oid;
data.tbid    = tbid;
data         = coco_func_data(data);

prob = coco_add_chart_data(prob, data.tbid, [], []);
prob = coco_add_func(prob, data.tbid, @func_FDF, data, 'zero', 'f+df', ...
  'x0', [sol.x0 ; sol.p0], 't0', sol.t0, 'passChart');

fid  = coco_get_id(data.tbid, 'update');
prob = coco_add_slot(prob, fid, @continex_update, data, 'update');
fid  = coco_get_id(data.tbid, 'set_jac');
prob = coco_add_slot(prob, fid, @continex_set_jac, data, 'corr_sample');
fid  = coco_get_id(data.tbid, 'force_FDM');
prob = coco_add_slot(prob, fid, @continex_force_FDM, data, 'continex_force_FDM');

% add external parameters
if ~isempty(data.pnames)
  fid  = coco_get_id(data.tbid, 'pars');
  xidx = coco_get_func_data(prob, data.tbid, 'xidx');
  prob = coco_add_pars(prob, fid, xidx(data.p_idx), data.pnames);
end

% initialise output slots
fid  = coco_get_id(data.tbid, 'bddat');
prob = coco_add_slot(prob, fid, @continex_bddat, data, 'bddat');
fid  = coco_get_id(data.tbid, 'save');
prob = coco_add_slot(prob, fid, @coco_save_data, data, 'save_full');
fid  = coco_get_id(data.tbid, 'corr_print');
prob = coco_add_slot(prob, fid, @continex_print, data, 'corr_print');
fid  = coco_get_id(data.tbid, 'cont_print');
prob = coco_add_slot(prob, fid, @continex_print, data, 'cont_print');

fid  = coco_get_id(data.tbid, 'plot');
prob = coco_add_slot(prob, fid, @continex_plot, data, 'FSM_state_co_flush_end');

end

function data = continex_get_settings(prob, tbid, data) %#ok<INUSL>

defaults.NJac  = 0 ; % recompute Jacobian using FDM every NJac steps, 0 = off
defaults.JacH  = [0.01 0.01]; % step sizes for FDM
defaults.xpar  = ''; % name of x-axis
defaults.ypar  = ''; % name of y-axis
defaults.phan  = [] ; % plot handle
defaults.pfunc = @(x) sqrt(sum(x.^2,1)); % plot function
defaults.scale = @(x,p) 1/(1+norm(x)); % scaling of residuum

data.continex = coco_merge(defaults, coco_get(prob, 'continex'));

if numel(data.continex.JacH)<2
  data.continex.JacH = [data.continex.JacH data.continex.JacH];
end

end

function [data chart F J] = func_FDF(prob, data, chart, xp)
% continex zero problem

% F
x = xp(data.x_idx,1);
p = xp(data.p_idx,1);

if isfield(data, 'F_sv_xp') && all( xp-data.F_sv_xp == 0 )
  F = data.F_sv_F;
else
  if isempty(data.fpars)
    F = data.F(x, p);
  else
    % bug: can be simplified if multiple substitutions are implemented in
    % coco_func_data.subsref
    fpars = data.fpars;
    [data.fpars{1} F] = data.F(fpars{1}, x, p, fpars{2:end});
  end
  data.F_sv_xp  = xp;
  data.F_sv_F   = F;
end

if nargout<=3
  F = data.scale*F;
  return
end

% DFDX
data.DFCount = data.DFCount+1;

if data.NJac>0 && data.jac && data.DFCount>data.NJac
  % compute FDM approximation of Jacobian
  data.DFCount = 1;
  J = zeros( numel(F), numel(xp) );
  for i=1:numel(xp)
    xx = xp;
    if ismember(i, data.x_idx)
      h  = data.JacH(1)*(1+abs(xx(i)));
    else
      h  = data.JacH(2)*(1+abs(xx(i)));
    end
    xx(i) = xx(i)+h;
    x1 = xx(data.x_idx,1);
    p1 = xx(data.p_idx,1);
    if isempty(data.fpars)
      F1 = data.F(x1, p1);
    else
      % bug: can be simplified if multiple substitutions are implemented in
      % coco_func_data.subsref
      fpars = data.fpars;
      [data.fpars{1} F1] = data.F(fpars{1}, x1, p1, fpars{2:end});
    end
    J(:,i) = (F1-F)/h;
  end
  data.sv_J  = J;
  data.sv_F  = F;
  data.sv_xp = [x;p];
else
  if isfield(data, 'sv_J')
    % compute Broyden update of Jacobian
    s = [x;p] - data.sv_xp;
    if norm(s)>=max(eps, 0.1*prob.corr.TOL)
      y = F - data.sv_F;
      J = data.sv_J + ( (y-data.sv_J*s) * s' )/(s'*s);
      
      data.sv_J  = J;
      data.sv_F  = F;
      data.sv_xp = [x;p];
    else
      J = data.sv_J;
    end
  else
    % initialise data for Broyden updates
    data.sv_F  = F;
    J          = eye( numel(data.sv_F), numel(xp) );
    data.sv_J  = J;
    data.sv_xp = [x;p];
  end
end

F = data.scale*F;
J = data.scale*J;

cdata       = coco_get_chart_data(chart, data.tbid);
cdata.scale = data.scale;
cdata.J     = J;
chart       = coco_set_chart_data(chart, data.tbid, cdata);

end

function data = continex_update(prob, data, cseg, varargin)

xidx       = coco_get_func_data(prob, data.tbid, 'xidx');
x0         = cseg.src_chart.x(xidx);
data.scale = data.scaleF(x0(data.x_idx), x0(data.p_idx));
data.jac   = true;

end

function data = continex_set_jac(prob, data, varargin) %#ok<INUSL>
% Disable FDM recomputation of Jacobian while Newton is sampling points.
data.data.jac = false;
end

function data = continex_force_FDM(prob, data, varargin) %#ok<INUSL>
% Force FDM recomputation of Jacobian in next call to FDF.
data.data.DFCount=data.data.NJac+1;
end

function [data res] = continex_bddat(prob, data, cmd, chart)

switch cmd
  case 'init'
    if isempty(data.oid)
      res = { 'FC' 'PARS' '||C||' };
    else
      res = coco_get_id(data.tbid, { 'FC' 'PARS' '||C||' });
    end
  case 'data'
    [fdata xidx] = coco_get_func_data(prob, data.tbid, 'data', 'xidx');
    x   = chart.x(xidx);
    c   = x(fdata.x_idx);
    res = { c x(fdata.p_idx) norm(c) };
end

end

function data = continex_print(prob, data, command, prio, chart, x) %#ok<INUSL>

switch command
  case 'init'
    coco_print(prob, prio, '%10s', '||C||');
  case 'data'
    [fdata xidx] = coco_get_func_data(prob, data.tbid, 'data', 'xidx');
    coco_print(prob, prio, '%10.2e', norm(x(xidx(fdata.x_idx))));
end

end

function data = continex_plot(prob, data)

continex = data.continex;
if any(cellfun(@isempty, {continex.phan continex.xpar continex.ypar}))
  return
end

x = coco_bd_col(prob.bd, continex.xpar);
y = coco_bd_col(prob.bd, continex.ypar);
if isempty(x)
  return
end
if isfield(data, 'bdhan')
  delete(data.bdhan);
end
hold(continex.phan, 'on')
data.bdhan = plot(continex.phan, x, continex.pfunc(y), 'b.-');
hold(continex.phan, 'off')
grid(continex.phan, 'on')
drawnow

end
