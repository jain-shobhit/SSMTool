function [prob, data] = po_construct_tst(prob, data)
%PO_CONSTRUCT_TST   Add codim-1 test functions.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_construct_tst.m 3053 2017-10-26 01:56:00Z hdankowicz $

data = init_data(data);
po   = data.po;
fid  = data.po_tst.fid;

if any([po.SN po.PD po.TR po.USTAB])
  pids = coco_get_id(fid, {'SN' 'PD' 'NS' 'USTAB'});
  
  data.sh.po_M = [];
  prob = coco_add_chart_data(prob, fid, [], []);
  void = coco_get_id(data.cid, 'test');
  prob = coco_add_func(prob, fid, @test_SN_PD_TR, data, ...
    'regular', pids, 'remesh', @remesh, 'requires', void, 'passChart');
  data.protect('po_M');
  data.no_save = [ data.no_save { 'po_M' } ];
  
  prob = coco_add_slot(prob, fid, @bddat, data, 'bddat');
  
  if po.SN
    prob = coco_add_event(prob, @evhan_SN, data, 'SP', pids{1}, 0);
  end
  if po.PD
    prob = coco_add_event(prob, @evhan_PD, data, 'SP', pids{2}, 0);
  end
  if po.TR
    prob = coco_add_event(prob, @evhan_TR, data, 'SP', pids{3}, 0);
  end
end

end

function data = init_data(data)
%INIT_TST_DATA   Initialize test function data.

if data.po.TR
  xdim        = data.po_orb.dim;
  I           = triu(true(xdim),1);
  A           = repmat((1:xdim)', [1 xdim]);
  tst.la_idx1 = A(I);
  A           = A';
  tst.la_idx2 = A(I);
else
  tst.la_idx1 = [];
  tst.la_idx2 = [];
end
tst.fid = coco_get_id(data.po_orb.fid, 'test');

data.po_tst  = tst;
data.no_save = [ data.no_save { 'po_tst.la_idx1' 'po_tst.la_idx2' } ];

end

function [data, chart, y] = test_SN_PD_TR(prob, data, chart, u) %#ok<INUSD>
%TEST_SN_PD_TR   Monitor functions for stability and codim-1 bifurcations.
%
% y(1) : Saddle-node points
% y(2) : Period-doubling points
% y(3) : Neimark-Sacker and neutral saddle point
% y(4) : Number of unstable eigenvalues

tst = data.po_tst;

cdata = coco_get_chart_data(chart, tst.fid);
if ~isempty(cdata) && isfield(cdata, 'la') % check this!
  la = cdata.la;
else
  fdata  = coco_get_func_data(prob, data.cid, 'data');
  ctst   = fdata.coll_tst;
  M0     = ctst.M(ctst.M0_idx,:);
  M1     = ctst.M(ctst.M1_idx,:);
  M      = M1/M0;
  [D, la] = eig(full(M)); %#ok<ASGLU>
  la     = diag(la); % check not getting D!
  if data.ode.autonomous
    [~, idx] = sort(abs(la-1));
    la = la(idx(2:end));
  end
  cdata.la = la;
  chart  = coco_set_chart_data(chart, tst.fid, cdata);
  data.po_M = M;
end

% Stability indicator
y(4,1) = sum(abs(la)>=1);

% Saddle-node points
y(1,1) = real(prod(la-1));

% Period-doubling points
y(2,1) = real(prod(la+1));

% Torus bifurcations and neutral saddle points
if numel(la)>1
  la  = la(tst.la_idx1).*la(tst.la_idx2);
  y(3,1) = real(prod(la-1));
else
  y(3,1) = 1;
end

end

function [data, cseg, msg] = evhan_SN(prob, data, cseg, cmd, msg)
%EVHAN_SN   Saddle-node bifurcation event handler.
%
% Save restart data.

tst = data.po_tst;

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
    cdata = coco_get_chart_data(cseg.curr_chart, tst.fid);
    if data.ode.autonomous
      [fdata, uidx] = coco_get_func_data(prob, data.cid, 'data', 'uidx');
      u     = cseg.curr_chart.x(uidx);
      maps  = fdata.coll_seg.maps;
      f0    = fdata.ode_F(fdata, 0, u(maps.x0_idx), u(maps.p_idx));
      M     = [ data.po_M-eye(numel(f0)) f0 ; f0' 0];
      [V, D]   = eig(M);
      [~, idx] = min(abs(diag(D)));
      V     = V(:,idx); % generalized eigenvector for double eigenvalue at 1
      V     = V/norm(V(1:end-1));
      sn.v = V(1:end-1);
      sn.b = V(end);
      sn.M = data.po_M;
    else
      [V, D]   = eig(data.po_M);
      [~, idx] = min(abs(diag(D)-1));
      sn.v = V(:,idx);
      sn.b = [];
    end
    cdata.sn   = sn;
    cseg.curr_chart = coco_set_chart_data(cseg.curr_chart, tst.fid, cdata);
    msg.action = 'add';
    msg.finish = true;
end

end

function [data, cseg, msg] = evhan_TR(prob, data, cseg, cmd, msg)
%EVHAN_TR   Neimark-Sacker bifurcation event handler.
%
% Distinguish between torus bifurcations (two complex conjugate eigenvalues
% on the unit circle) and neutral saddle points (two real eigenvalues whose
% product equals 1). For torus bifurcation points, save restart data.

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

tst = data.po_tst;

switch cmd
  
  case 'init'
    if isfield(msg, 'finish') || strcmp(msg.action, 'warn')
      msg.action = 'finish';
    elseif strcmp(msg.action, 'locate')
      msg.action = 'warn';
    else
      cdata = coco_get_chart_data(cseg.ptlist{1}, tst.fid);
      la0   = cdata.la;
      cdata = coco_get_chart_data(cseg.ptlist{end}, tst.fid);
      la1   = cdata.la;
      switch abs(sum(sign(abs(la0)-1))-sum(sign(abs(la1)-1)))
        case 4
          msg.point_type = 'TR';
          msg.action     = 'locate';
        case 0
          msg.point_type = 'NSA';
          if data.po.NSA
            msg.action   = 'locate';
          else
            msg.action   = 'finish';
          end
        otherwise
          msg.point_type = '?TR?';
          msg.action     = 'warn';
          msg.wmsg       = 'could not determine type of bifurcation point';
      end
      msg.idx = 1;
    end
    
  case 'check'
    cdata = coco_get_chart_data(cseg.curr_chart, tst.fid);
    [V,D] = eig(data.po_M);
    D = diag(D);
    if data.ode.autonomous
      [~, idx] = sort(abs(D-1));
      D = D(idx(2:end));
      V = V(:,idx(2:end));
    end
    [~, tridx] = min(abs(abs(D)-1));
    evec = rotate_evec(prob, V(:,tridx));
    tr.v = [real(evec) imag(evec)];
    th   = angle(D(tridx));
    tr.a = cos(th);
    tr.b = sin(th);
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

tst = data.po_tst;

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
    cdata = coco_get_chart_data(cseg.curr_chart, tst.fid);
    [V,D]      = eig(data.po_M);
    [~, idx]   = min(abs(diag(D)+1));
    pd.v = V(:,idx);

    [fdata, uidx] = coco_get_func_data(prob, data.cid, 'data', 'uidx');
    maps  = fdata.coll_seg.maps;
    mesh  = fdata.coll_seg.mesh;
    ctst  = fdata.coll_tst;
    u = cseg.curr_chart.x(uidx);
    x = u(maps.xbp_idx);
    T = u(maps.T_idx);
    p = u(maps.p_idx);
    
    % Construct initial solution guess for period-doubled orbit
    t0    = mesh.tbp(maps.tbp_idx)*T;
    t0    = [t0; T+t0(2:end)];
    xp1   = reshape(x+0.01*ctst.M*V(:,idx), maps.xbp_shp)';
    xp1   = xp1(maps.tbp_idx,:);
    xp2   = reshape(x-0.01*ctst.M*V(:,idx), maps.xbp_shp)';
    xp2   = xp2(maps.tbp_idx,:);
    x0    = [xp1; xp2(2:end,:)];
    pd.x0 = x0;
    pd.t0 = t0;
    pd.p0 = p;
    cdata.pd   = pd;
    cseg.curr_chart = coco_set_chart_data(cseg.curr_chart, tst.fid, cdata);
    msg.action = 'add';
    msg.finish = true;
end

end

function [prob, stat, xtr] = remesh(prob, data, chart, ub, Vb) %#ok<INUSD>

xtr  = [];
prob = coco_change_func(prob, data);
stat = 'success';

end

function [data, res] = bddat(prob, data, command, varargin) %#ok<INUSL>
%BDDAT   Append Floquet multipliers to BD.

res = {};
switch command
  case 'init'
    res = coco_get_id(data.oid,'eigs');
  case 'data'
    chart = varargin{1}; % Current chart
    cdata = coco_get_chart_data(chart, data.po_tst.fid);
    if ~isempty(cdata) && isfield(cdata, 'la')
      res = cdata.la;
    end
end

end