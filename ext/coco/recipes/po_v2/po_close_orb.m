function prob = po_close_orb(prob, tbid, data)
%PO_CLOSE_ORB   Append an instance of 'po' to problem.
%
% Add boundary and integral phase conditions to an instance of 'coll'.
%
% Differs from po_v1 by providing support for bifurcation analysis of
% periodic orbits.
%
% PROB = PO_CLOSE_ORB(PROB, TBID, DATA)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_close_orb.m 2839 2015-03-05 17:09:01Z fschild $

data.tbid = tbid;
data = coco_func_data(data); % Convert to func_data class for shared access
prob = coco_add_slot(prob, tbid, @po_update, data, 'update'); % Phase condition update function
segtbid      = coco_get_id(tbid, 'seg.coll'); % Create 'coll' toolbox instance identifier
[fdata uidx] = coco_get_func_data(prob, segtbid, 'data', 'uidx'); % Extract 'coll' data structure and context-dependent index array
prob = coco_add_func(prob, tbid, @po_F, @po_DFDU, data, 'zero', ...
  'uidx', uidx(fdata.xbp_idx));
fid  = coco_get_id(tbid, 'period');
prob = coco_add_pars(prob, fid, uidx(fdata.T_idx), fid, 'active'); % Monitor period
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

if data.po.bifus % Optional monitor function for bifurcation detection
  segoid = coco_get_id(tbid, 'seg');
  prob = var_coll_add(prob, segoid); % Append variational equations
  data.var_id = coco_get_id(segoid, 'var');
  tfid = coco_get_id(tbid, 'test');
  data.tfid = tfid;
  tfps = coco_get_id(tfid, {'SN' 'PD' 'NS' 'stab'});
  prob = coco_add_chart_data(prob, tfid, [], []); % Allocate chart data for eigenvalue information
  prob = coco_add_func(prob, tfid, @po_TF, data, ...
    'regular', tfps, 'requires', data.var_id, 'passChart'); % Chart is passed to monitor function
  prob = coco_add_event(prob, 'SN', tfps{1}, 0); % SN - event type
  data.efid = coco_get_id(tbid, 'PD');
  prob = coco_add_chart_data(prob, data.efid, [], []); % Allocate chart data for branch switching at period doubling
  prob = coco_add_event(prob, @po_evhan_PD, data, tfps{2}, 0); % Period-doubling event handler
  prob = coco_add_event(prob, @po_evhan_NS, data, tfps{3}, 0); % Neimark-Sacker event handler
end

end

function [data y] = po_F(prob, data, u)
%PO_F   Evaluate boundary conditions and phase condition.
% 
% Integral phase condition is expressed as discretized linear operator.
%
% Identical to po_v1.

x0 = u(data.x0_idx); % Extract trajectory end point at t=0
x1 = u(data.x1_idx); % Extract trajectory end point at t=1

y = [x0-x1; data.xp0*u];

end

function [data J] = po_DFDU(prob, data, u)
%PO_DFDU   Evaluate linearization of boundary conditions and phase condition.
% 
% Integral phase condition is expressed as discretized linear operator.
%
% Identical to po_v1.

J = data.J;
end

function data = po_update(prob, data, cseg, varargin)
%PO_UPDATE   Update discretized linear operator.
%
% Use information about current solution to update discretization of linear
% operator corresponding to integral phase condition.
%
% Identical to po_v1.

fid           = coco_get_id(data.tbid, 'seg.coll');
[fdata uidx]  = coco_get_func_data(prob, fid, 'data', 'uidx'); % 'coll' toolbox data and context-dependent index set
u             = cseg.src_chart.x; % Current chart
data.xp0      = u(uidx(fdata.xbp_idx))'*data.intfac;
data.J(end,:) = data.xp0;

end

function [data chart y] = po_TF(prob, data, chart, u)
%PO_TF   Monitor function: bifurcations and stability.
%
% Return argument y: First component equals 0 when at least one eigenvalue
% equals 1. Second component equals 0 when at least one eigenvalue equals
% -1. Third component equals 0 when the product of two eigenvalues equals
% 1. Fourth component counts the number of eigenvalues outside the unit
% circle.
 
cdata = coco_get_chart_data(chart, data.tfid); % Read chart data
if ~isempty(cdata) && isfield(cdata, 'la')
  la = cdata.la;
else
  fdata    = coco_get_func_data(prob, data.var_id, 'data'); % Extract data of variational toolbox
  M        = fdata.M(fdata.M1_idx,:); % Jacobian of flow at t=1
  la       = eig(M);
  [vv idx] = sort(abs(la-1));
  la       = la(idx(2:end)); % Omit trivial eigenvalue at 1
  chart    = coco_set_chart_data(chart, data.tfid, struct('la', la)); % Write chart data
end
y(1,1) = prod(la-1);
y(2,1) = prod(la+1);
if numel(la)>1
  NS_TF  = la(data.la_idx1).*la(data.la_idx2); % Compute all products of pairs of distinct eigenvalues
  y(3,1) = prod(NS_TF(:)-1);
else
  y(3,1) = 1;
end
y(4,1) = sum(abs(la)>1);

end

function [data cseg msg] = po_evhan_NS(prob, data, cseg, cmd, msg)
%PO_EVHAN_NS   Neimark-Sacker bifurcation event handler.
%
% Distinguish between Neimark-Sacker bifurcations (two complex conjugate
% eigenvalues on the unit circle) and neutral saddle points (two reciprocal
% eigenvalues off the unit circle).

switch cmd
  case 'init'
    if isfield(msg, 'finish') || strcmp(msg.action, 'warn')
      msg.action = 'finish';
    elseif strcmp(msg.action, 'locate') % After unsuccessful location
      msg.action = 'warn';
    else % After detection
      cdata = coco_get_chart_data(cseg.ptlist{1}, data.tfid);
      la0 = cdata.la;
      cdata = coco_get_chart_data(cseg.ptlist{end}, data.tfid);
      la1 = cdata.la;
      switch abs(sum(sign(abs(la0)-1))-sum(sign(abs(la1)-1)))
        case 4 % Two eigenvalues cross the unit circle
          msg.point_type = 'NS';
          msg.action     = 'locate';
        case 0 % No eigenvalue crosses the unit circle
          msg.point_type = 'NSad';
          if data.po.NSad % Optional detection
            msg.action   = 'locate';
          else
            msg.action   = 'finish';
          end
        otherwise
          msg.point_type = 'NS';
          msg.action     = 'warn';
          msg.wmsg       = 'could not determine type of event';
      end
      msg.idx = 1;
    end
  case 'check' % Add special point to curve segment
    msg.action = 'add';
    msg.finish = true;
end

end

function [data cseg msg] = po_evhan_PD(prob, data, cseg, cmd, msg)
%PO_EVHAN_PD   Period-doubling bifurcation event handler.
%
% Support branch-switching at period-doubling points.

switch cmd
  case 'init'
    if isfield(msg, 'finish') || strcmp(msg.action, 'warn')
      msg.action = 'finish';
    elseif strcmp(msg.action, 'locate') % After unsuccessful location
      msg.action = 'warn';
    else % After detection
      msg.point_type = 'PD';
      msg.action     = 'locate';
      msg.idx = 1;
    end
  case 'check' % Add special point to curve segment
    msg.action = 'add';
    msg.finish = true;
    vdata = coco_get_func_data(prob, data.var_id, 'data');
    [uidx fdata] = coco_get_func_data(prob, vdata.tbid, ...
      'uidx', 'data'); % Extract context-dependent index array and variational toolbox data
    chart = cseg.curr_chart; % Current chart
    u = chart.x(uidx);
    x = u(fdata.xbp_idx); % Extract segment basepoint values
    T = u(fdata.T_idx);   % Extract interval length
    p = u(fdata.p_idx);   % Extract problem parameters
    
    M     = vdata.M;      % Solution to variational problem
    M1    = M(vdata.M1_idx,:);
    [v d] = eig(M1);      % Floquet multipliers. Bug: get chart data!
    [m i] = min(diag(d)+1);
    xp1   = reshape(x+0.01*M*v(:,i), fdata.xbp_shp)'; % Construct initial solution guess for period-doubled orbit
    xp1   = xp1(fdata.tbp_idx,:);
    t1    = fdata.tbp(fdata.tbp_idx)*T;
    xp2   = reshape(x-0.01*M*v(:,i), fdata.xbp_shp)'; % Construct initial solution guess for period-doubled orbit
    xp2   = xp2(fdata.tbp_idx,:);
    t2    = fdata.tbp(fdata.tbp_idx)*T; % Construct initial solution guess for period-doubled orbit
    x0    = [xp1; xp2(2:end,:)];
    t0    = [t1; T+t2(2:end)];
    chart = coco_set_chart_data(chart, data.efid, ...
      struct('pd_x0', x0, 'pd_t0', t0, 'pd_p', p)); % Write chart data
    cseg.curr_chart = chart;
end

end
