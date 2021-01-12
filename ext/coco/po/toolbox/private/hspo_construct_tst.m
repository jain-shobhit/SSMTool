function [prob, data] = hspo_construct_tst(prob, data)
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

data = init_data(prob, data);
hspo = data.hspo;
tst  = data.hspo_tst;
fid  = tst.fid;

if any([hspo.SN hspo.PD hspo.TR hspo.USTAB])
  pids = coco_get_id(fid, {'SN' 'PD' 'TR' 'USTAB'});
  
  data.sh.hspo_M = [];
  prob = coco_add_chart_data(prob, fid, [], []);
  
  [fdata, uidx] = coco_get_func_data(prob, data.bvid, 'data', 'uidx');
  prob = coco_add_func(prob, fid, @test_SN_PD_TR, data, ...
    'regular', pids, 'uidx', ...
    [uidx(fdata.bvp_bc.x1_idx); uidx(fdata.bvp_bc.p_idx)], ...
    'requires', tst.vids, 'passChart');
  
  data.protect('hspo_M');
  data.no_save = [ data.no_save { 'hspo_M' } ];
  
  prob = coco_add_slot(prob, fid, @bddat, data, 'bddat');
  
  if hspo.SN
    prob = coco_add_event(prob, @evhan_SN, data, 'SP', pids{1}, 0);
  end
  if hspo.PD
    prob = coco_add_event(prob, @evhan_PD, data, 'SP', pids{2}, 0);
  end
  if hspo.TR
    prob = coco_add_event(prob, @evhan_TR, data, 'SP', pids{3}, 0);
  end
end

end

function data = init_data(prob, data)
%INIT_TST_DATA   Initialize test function data.

fdata = coco_get_func_data(prob, data.bvid, 'data');
nsegs = fdata.nsegs;
vids  = cell(1,nsegs);
for i=1:nsegs
  vids{i} = coco_get_id(fdata.cids{i}, 'test');
end
tst.nsegs = nsegs;
tst.cids  = fdata.cids;
tst.vids  = vids;
tst.p_idx = numel(fdata.bvp_bc.x1_idx)+(1:numel(fdata.bvp_bc.p_idx));

if data.hspo.TR
  gdata = coco_get_func_data(prob, fdata.cids{1}, 'data');
  xdim        = gdata.xdim-1;
  I           = triu(true(xdim),1);
  A           = repmat((1:xdim)', [1 xdim]);
  tst.la_idx1 = A(I);
  A           = A';
  tst.la_idx2 = A(I);
else
  tst.la_idx1 = [];
  tst.la_idx2 = [];
end
tst.fid = coco_get_id(data.oid, 'hspo.test');

data.hspo_tst  = tst;
data.no_save = [ data.no_save { 'hspo_tst.la_idx1' 'hspo_tst.la_idx2' } ];

end

function [data, chart, y] = test_SN_PD_TR(prob, data, chart, u)
%TEST_SN_PD_TR   Monitor functions for stability and codim-1 bifurcations.
%
% y(1) : Saddle-node points
% y(2) : Period-doubling points
% y(3) : Neimark-Sacker and neutral saddle point
% y(4) : Number of unstable eigenvalues

pr  = data.pr;
tst = pr.hspo_tst;

cdata = coco_get_chart_data(chart, tst.fid); % Read chart data
if ~isempty(cdata) && isfield(cdata, 'la')
  la = cdata.la;
else
  P = hspo_P(prob, pr, u);
  M = P{1};
  for i=2:tst.nsegs
    M = P{i}*M;
  end
  la = eig(M);
  [~, idx] = sort(abs(la));
  la = la(idx(2:end));
  cdata.la = la;
  cdata.P  = P;
  chart = coco_set_chart_data(chart, tst.fid, cdata);
  data.hspo_M = M;
end

% Stability indicator
y(4,1) = sum(abs(la)>=1);

% Saddle-node points
y(1,1) = real(prod(la-1));

% Period-doubling points
y(2,1) = real(prod(la+1));

% Torus bifurcations and neutral saddle points
if numel(la)>1
  la     = la(tst.la_idx1).*la(tst.la_idx2); % Compute all products of pairs of distinct eigenvalues
  y(3,1) = real(prod(la-1))*1e10;
else
  y(3,1) = 1;
end

end

function P = hspo_P(prob, data, u)
%HSPO_P   Compute collection of transfer matrices.
%
% For each segment, extract Jacobian of time-T map and premultiply by
% saltation matrix (correcting for difference in time-of-flight to event
% surface, and the imposition of the reset) to obtain transfer matrix.
%
% P = HSPO_P(PROB, DATA, U)
%
% P    - Cell array of transfer matrices.
% PROB - Continuation problem structure.
% DATA - hspo_mult_eigs_bddat function data.
% U    - Array of end points of segments at t=1 and problem parameters.

bc  = data.hspo_orb;
tst = data.hspo_tst;

P = cell(1,bc.nsegs);
p = u(tst.p_idx); % Problem parameters
for i=1:bc.nsegs
  fdata = coco_get_func_data(prob, tst.cids{i}, 'data');
  ctst  = fdata.coll_tst;
  M0    = ctst.M(ctst.M0_idx,:);
  M1    = ctst.M(ctst.M1_idx,:);
  
  x     = u(bc.x1_idx{i}); % Segment end point at t=1
  dim   = numel(x);
  fs    = fdata.ode_F(fdata, 0, x, p);
  dhdx  = data.event_DFDX(data, x, p, data.events{i});
  dgdx  = data.reset_DFDX(data, x, p, data.resets{i});
  P{i}  = dgdx * (eye(dim) - (fs*dhdx)/(dhdx*fs)) * M1/M0; % Multiplication by saltation matrix
end

end

function [data, cseg, msg] = evhan_SN(prob, data, cseg, cmd, msg)
%EVHAN_SN   Saddle-node bifurcation event handler.
%
% Save restart data.

tst = data.hspo_tst;
nsegs = tst.nsegs;

switch cmd
  case 'init'
    if isfield(msg, 'finish') || strcmp(msg.action, 'warn')
      msg.action = 'finish';
    elseif strcmp(msg.action, 'locate') % After unsuccessful location
      msg.action = 'warn';
    else % After detection
      msg.point_type = 'SN';
      msg.action     = 'locate';
      msg.idx = 1;
    end
  case 'check'
    cdata    = coco_get_chart_data(cseg.curr_chart, tst.fid);
    [V, D]   = eig(data.hspo_M);
    [~, idx] = min(abs(diag(D)-1));
    
    v    = V(:,idx); % Initial perturbation
    sn.v = cell(1,nsegs);
    for i=1:nsegs
      sn.v{i} = v;
      v       = cdata.P{i}*v; % Map to next initial perturbation
    end
    cdata.sn   = sn;
    cseg.curr_chart = coco_set_chart_data(cseg.curr_chart, tst.fid, cdata);
    msg.action = 'add';
    msg.finish = true;
end

end

function [data, cseg, msg] = evhan_TR(prob, data, cseg, cmd, msg)
%EVHAN_TR   Torus bifurcation event handler.
%
% Distinguish between torus bifurcations (two complex conjugate eigenvalues
% on the unit circle) and neutral saddle points (two reciprocal eigenvalues
% off the unit circle). For torus bifurcation points, save restart data.

  function z = rotate_evec(prob, xx)
    x  = real(xx);
    y  = imag(xx);
    
    if abs(x'*y)<10*prob.corr.TOL
      z = xx;
      return
    end
    
    ga = (y'*y-x'*x)/(x'*y);
    if ga>0
      p0 = acot(1+ga);
      p1 = pi/4;
    else
      p0 = pi/4;
      p1 = atan(1-ga);
    end
    
    phi = (p0+p1)/2;
    while abs(p1-p0)>prob.corr.TOL
      if (cot(p0)-tan(p0)-ga)*(cot(phi)-tan(phi)-ga)<0
        p1 = phi;
      else
        p0 = phi;
      end
      phi = (p0+p1)/2;
    end
    z = (cos(phi)+1i*sin(phi))*xx;
  end

tst = data.hspo_tst;
nsegs = tst.nsegs;

switch cmd
  case 'init'
    if isfield(msg, 'finish') || strcmp(msg.action, 'warn')
      msg.action = 'finish';
    elseif strcmp(msg.action, 'locate') % After unsuccessful location
      msg.action = 'warn';
    else % After detection
      cdata = coco_get_chart_data(cseg.ptlist{1}, tst.fid);
      la0 = cdata.la;
      cdata = coco_get_chart_data(cseg.ptlist{end}, tst.fid);
      la1 = cdata.la;
      switch abs(sum(sign(abs(la0)-1)) - sum(sign(abs(la1)-1)))
        case 4 % Two eigenvalues cross the unit circle
          msg.point_type = 'TR';
          msg.action     = 'locate';
        case 0 % No eigenvalue crosses the unit circle
          msg.point_type = 'NSA';
          if data.hspo.NSA % Optional detection
            msg.action   = 'locate';
          else
            msg.action   = 'finish';
          end
        otherwise
          msg.point_type = 'TR';
          msg.action     = 'warn';
          msg.wmsg       = 'could not determine type of event';
      end
      msg.idx = 1;
    end
  case 'check' % Add special point to curve segment
    cdata    = coco_get_chart_data(cseg.curr_chart, tst.fid);
    [V, D]   = eig(data.hspo_M);
    D        = diag(D);
    [~, tridx] = min(abs(abs(D)-1));
    evec = rotate_evec(prob, V(:,tridx));
    th   = angle(D(tridx));
    tr.a = cos(th);
    tr.b = sin(th);
    
    v    = [real(evec) imag(evec)];
    tr.v = cell(1,nsegs);
    for i=1:nsegs
      tr.v{i} = v;
      v       = cdata.P{i}*v; % Map to next initial perturbation
    end
    cdata.tr   = tr;
    cseg.curr_chart = coco_set_chart_data(cseg.curr_chart, tst.fid, cdata);
    msg.action = 'add';
    msg.finish = true;
end

end

function [data, cseg, msg] = evhan_PD(prob, data, cseg, cmd, msg)
%EVHAN_PD   Period-doubling bifurcation event handler.
%
% Support branch-switching at period-doubling points.

pr    = data.pr;
tst   = pr.hspo_tst;
nsegs = tst.nsegs;

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
    cdata    = coco_get_chart_data(cseg.curr_chart, tst.fid);
    [V, D]   = eig(data.hspo_M);
    [~, idx] = min(abs(diag(D)+1));
    
    v    = V(:,idx); % Initial perturbation
    t0   = cell(1,nsegs);
    x10  = cell(1,nsegs);
    x20  = cell(1,nsegs);
    pd.v = cell(1,nsegs);
    for i=1:nsegs
      [uidx, fdata] = coco_get_func_data(prob, tst.cids{i}, 'uidx', 'data');
      maps = fdata.coll_seg.maps;
      mesh = fdata.coll_seg.mesh;
      ctst = fdata.coll_tst;
      u = cseg.curr_chart.x(uidx);
      x = u(maps.xbp_idx); % Extract segment basepoint values
      T = u(maps.T_idx);   % Extract interval length
      p = u(maps.p_idx);   % Extract problem parameters
      
      t0{i}   = mesh.tbp(maps.tbp_idx)*T;
      x1      = reshape(x+0.01*ctst.M*v, maps.xbp_shp)'; % First loop
      x10{i}  = x1(maps.tbp_idx,:);
      x2      = reshape(x-0.01*ctst.M*v, maps.xbp_shp)'; % Second loop
      x20{i}  = x2(maps.tbp_idx,:);
      pd.v{i} = v;
      v       = cdata.P{i}*v; % Map to next initial perturbation
    end
    pd.x0     = [x10 x20];
    pd.t0     = [t0 t0];
    pd.p0     = p;
    pd.modes  = [pr.modes  pr.modes];
    pd.events = [pr.events pr.events];
    pd.resets = [pr.resets pr.resets];
    cdata.pd  = pd;
    cseg.curr_chart = coco_set_chart_data(cseg.curr_chart, tst.fid, cdata);
    msg.action = 'add';
    msg.finish = true;
end

end

function [data, res] = bddat(prob, data, command, varargin) %#ok<INUSL>
%BDDAT   Append Floquet multipliers to BD.

res = {};
switch command
  case 'init'
    res = coco_get_id(data.oid,'eigs');
  case 'data'
    chart = varargin{1}; % Current chart
    cdata = coco_get_chart_data(chart, data.hspo_tst.fid);
    if ~isempty(cdata) && isfield(cdata, 'la')
      res = cdata.la;
    end
end

end