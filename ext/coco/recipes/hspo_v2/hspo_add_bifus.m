function prob = hspo_add_bifus(prob, oid, tbid, data)
%HSPO_ADD_BIFUS   Append bifurcation detection to problem.
%
% Append monitor functions and events, with associated event handlers.
%
% PROB = HSPO_ADD_BIFUS(PROB, OID, TBID, DATA)
%
% PROB - Continuation problem structure.
% OID  - Object instance identifier.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: hspo_add_bifus.m 2839 2015-03-05 17:09:01Z fschild $

cids = {};
vids = {};
msid = coco_get_id(oid, 'msbvp'); % 'msbvp' toolbox instance identifier
for i=1:numel(data.modes)
  soid = coco_get_id(msid,sprintf('seg%d', i));
  prob = var_coll_add(prob, soid); % Append variational problem to each segment
  cids = [cids {coco_get_id(soid, 'coll')}];
  vids = [vids {coco_get_id(soid, 'var')}];
end
data.msid = msid;
data.cids = cids; % 'coll' object instance identifiers
data.vids = vids; % 'varcoll' object instance identifiers
[fdata uidx] = coco_get_func_data(prob, msid, 'data', 'uidx'); % Extract 'msbvp' instance data and index array
data.p_idx = numel(fdata.x1_idx)+(1:numel(fdata.p_idx));
tfid = coco_get_id(tbid, 'test');
data.tfid = tfid;
tfps = coco_get_id(tfid, {'SN' 'PD' 'NS' 'stab'});
prob = coco_add_chart_data(prob, tfid, [], []); % Allocate chart data for eigenvalue information
prob = coco_add_func(prob, tfid, @hspo_TF, data, 'regular', tfps, ...
  'uidx', [uidx(fdata.x1_idx); uidx(fdata.p_idx)], ...
  'requires', vids, 'passChart'); % Chart is passed to monitor function
prob = coco_add_event(prob, 'SN', tfps{1}, 0); % SN - event type
data.efid = coco_get_id(tbid, 'PD');
prob = coco_add_chart_data(prob, data.efid, [], []); % Allocate chart data for branch switching at period doubling
prob = coco_add_event(prob, @hspo_evhan_PD, data, tfps{2}, 0); % Period-doubling event handler
prob = coco_add_event(prob, @hspo_evhan_NS, data, tfps{3}, 0); % Neimark-Sacker event handler

end

function [data chart y] = hspo_TF(prob, data, chart, u)
%HSPO_TF   Monitor function: bifurcations and stability.
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
  P = hspo_P(prob, data, u);
  M = P{1};
  for i=2:data.nsegs
    M = P{i}*M;
  end
  la = eig(M);
  chart = coco_set_chart_data(chart, data.tfid, ...
    struct('M', M, 'la', la, 'P', {P})); % Write chart data
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

function P = hspo_P(prob, data, u)
%HSPO_P   Compute collection of transfer matrices.
%
% For each segment, extract Jacobian of time-T map and premultiply by
% saltation matrix (correcting for difference in time-of-flight to event
% surface, and the imposition of the reset) to obtain transfer matrix.
%
% Identical to varcoll_v1.
%
% P = HSPO_P(PROB, DATA, U)
%
% P    - Cell array of transfer matrices.
% PROB - Continuation problem structure.
% DATA - hspo_mult_eigs_bddat function data.
% U    - Array of end points of segments at t=1 and problem parameters.

P = cell(1,data.nsegs);
p = u(data.p_idx); % Problem parameters
for i=1:data.nsegs
  fdata = coco_get_func_data(prob, data.vids{i}, 'data');
  M1    = fdata.M(fdata.M1_idx,:); % Jacobian of time-T map
  dim   = fdata.dim;
  
  fdata = coco_get_func_data(prob, data.cids{i}, 'data');
  x     = u(data.x1_idx{i}); % Segment end point at t=1
  fs    = fdata.fhan(x, p);
  if ~isempty(data.dhdxhan)
    hx = data.dhdxhan(x, p, data.events{i});
  else
    hx = coco_ezDFDX('f(x,p)v', data.hhan, x, p, data.events{i});
  end
  if ~isempty(data.dgdxhan)
    gx = data.dgdxhan(x, p, data.resets{i});
  else
    gx = coco_ezDFDX('f(x,p)v', data.ghan, x, p, data.resets{i});
  end
  P{i}  = gx*(eye(dim)-(fs*hx)/(hx*fs))*M1; % Multiplication by saltation matrix
end

end

function [data cseg msg] = hspo_evhan_NS(prob, data, cseg, cmd, msg)
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
      switch abs(sum(sign(abs(la0)-1)) - sum(sign(abs(la1)-1)))
        case 4 % Two eigenvalues cross the unit circle
          msg.point_type = 'NS';
          msg.action     = 'locate';
        case 0 % No eigenvalue crosses the unit circle
          msg.point_type = 'NSad';
          if data.hspo.NSad % Optional detection
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

function [data cseg msg] = hspo_evhan_PD(prob, data, cseg, cmd, msg)
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
    
    chart   = cseg.curr_chart; % Current chart
    cdata   = coco_get_chart_data(chart, data.tfid);
    M       = cdata.M;
    [v d]   = eig(M);
    [m idx] = min(diag(d)+1); % Find eigenvector corresponding to -1 eigenvalue
    v       = 0.01*v(:,idx); % Initial perturbation
    
    t0   = {};
    x0   = {};
    for i=1:data.nsegs % First loop
      vdata = coco_get_func_data(prob, data.vids{i}, 'data');
      [uidx fdata] = coco_get_func_data(prob, vdata.tbid, ...
        'uidx', 'data'); % Extract context-dependent index array and variational toolbox data
      u  = chart.x(uidx);
      x  = u(fdata.xbp_idx); % Extract segment basepoint values
      T  = u(fdata.T_idx);   % Extract interval length
      p  = u(fdata.p_idx);   % Extract problem parameters
      
      x  = reshape(x+vdata.M*v, fdata.xbp_shp)'; % Perturb segment
      t0 = [t0 {fdata.tbp(fdata.tbp_idx)*T}];
      x0 = [x0 {x(fdata.tbp_idx,:)}];
      v  = cdata.P{i}*v; % Map to next initial perturbation
    end
    for i=1:data.nsegs % Second loop
      vdata = coco_get_func_data(prob, data.vids{i}, 'data');
      [uidx fdata] = coco_get_func_data(prob, vdata.tbid, ...
        'uidx', 'data'); % Extract context-dependent index array and variational toolbox data
      u  = chart.x(uidx);
      x  = u(fdata.xbp_idx); % Extract segment basepoint values
      T  = u(fdata.T_idx);   % Extract interval length
      p  = u(fdata.p_idx);   % Extract problem parameters
      
      x  = reshape(x+vdata.M*v, fdata.xbp_shp)'; % Perturb segment
      t0 = [t0 {fdata.tbp(fdata.tbp_idx)*T}];
      x0 = [x0 {x(fdata.tbp_idx,:)}];
      v  = cdata.P{i}*v; % Map to next initial perturbation
    end
    
    pd        = struct('t0', {t0}, 'x0', {x0}, 'p', {p}); % Store new initial solution guess
    pd.modes  = [data.modes  data.modes];
    pd.events = [data.events data.events];
    pd.resets = [data.resets data.resets];
    chart = coco_set_chart_data(chart, data.efid, pd); % Write chart data
    cseg.curr_chart = chart;
end

end
