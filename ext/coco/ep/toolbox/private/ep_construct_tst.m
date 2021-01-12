function [prob, data] = ep_construct_tst(prob, data)
%EP_CONSTRUCT_TST   Add codim-1 test functions.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_construct_tst.m 3041 2017-10-04 14:14:18Z hdankowicz $

data = init_data(data);
ep   = data.ep;
fid  = data.ep_tst.fid;

if any([ep.USTAB ep.SN ep.HB])
  pids = coco_get_id(fid, {'SN' 'HB' 'USTAB'});

  data.sh.ep_X = [];
  prob = coco_add_chart_data(prob, fid, [], []);
  uidx = coco_get_func_data(prob, data.ep_eqn.fid, 'uidx');
  prob = coco_add_func(prob, fid, @test_SN_HB_NSA, data, ...
    'regular', pids, 'uidx', uidx, 'passChart');
  data.protect('ep_X');
  data.no_save = [ data.no_save { 'ep_X' } ];
  
  prob = coco_add_slot(prob, fid, @bddat, data, 'bddat');
  
  if ep.SN
    prob = coco_add_event(prob, 'SN', pids{1}, 0);
  end
  if ep.HB
    prob = coco_add_event(prob, @evhan_HB_NSA, data, 'SP', pids{2}, 0);
  end
end

end

function data = init_data(data)
%INIT_TST_DATA   Initialize test function data.

if data.ep.HB
  xdim        = data.xdim;
  I           = triu(true(xdim),1);
  A           = repmat((1:xdim)', [1 xdim]);
  tst.la_idx1 = A(I);
  A           = A';
  tst.la_idx2 = A(I);
else
  tst.la_idx1 = [];
  tst.la_idx2 = [];
end
tst.fid = coco_get_id(data.ep_eqn.fid, 'test');

data.ep_tst  = tst;
data.no_save = [ data.no_save { 'ep_tst.la_idx1' 'ep_tst.la_idx2' } ];

end

function [data, chart, y] = test_SN_HB_NSA(~, data, chart, u)
%TEST_SN_HB_NSA   Monitor functions for stability and codim-1 bifurcations.
%
% y(1) : Saddle-node points
% y(2) : Hopf- and neutral saddle point
% y(3) : Number of unstable eigenvalues

pr  = data.pr;
eqn = pr.ep_eqn;
tst = pr.ep_tst;

cdata = coco_get_chart_data(chart, tst.fid);
if ~isempty(cdata) && isfield(cdata, 'la')
  la = cdata.la;
else
  x      = u(eqn.x_idx);
  p      = u(eqn.p_idx);
  Jx     = pr.ode_DFDX(pr, [], x, p);
  [X, la] = eig(full(Jx));
  la     = diag(la);
  cdata.la = la;
  chart  = coco_set_chart_data(chart, tst.fid, cdata);
  data.ep_X = X;
end

% Stability indicator
y(3,1) = sum(real(la)>=0);

% Saddle-node points
rla    = real(la);
s      = sign(rla);
w      = min(abs(la));
y(1,1) = w*prod(s);

% Hopf- and neutral saddle points
la  = la(tst.la_idx1)+la(tst.la_idx2);
w   = min([1 min(abs(la))]);
rla = real(la);
rla = (rla>=0)-(rla<0);
la  = complex(rla,imag(la));
la  = la./abs(la);
s   = sign(real(prod(la)));

y(2,1) = w*s;

end

function [data, cseg, msg] = evhan_HB_NSA(~, data, cseg, cmd, msg)
%EVHAN_HB_NSA   Hopf bifurcation event handler.
%
% Distinguish between Hopf bifurcations (two complex conjugate, purely
% imaginary eigenvalues) and neutral saddle points (two equal in magnitude,
% opposite in sign, real eigenvalues). For Hopf bifurcation points, save
% restart data.

pr  = data.pr;
fid = pr.ep_tst.fid;

switch cmd
  
  case 'init'
    if isfield(msg, 'finish') || strcmp(msg.action, 'warn')
      msg.action = 'finish';
    elseif strcmp(msg.action, 'locate')
      msg.action = 'warn';
    else
      cdata = coco_get_chart_data(cseg.ptlist{1}, fid);
      la0   = cdata.la;
      cdata = coco_get_chart_data(cseg.ptlist{end}, fid);
      la1   = cdata.la;
      switch abs(sum(sign(real(la0)))-sum(sign(real(la1))))
        case 4
          msg.point_type = 'HB';
          msg.action     = 'locate';
        case 0
          msg.point_type = 'NSA';
          if pr.ep.NSA
            msg.action   = 'locate';
          else
            msg.action   = 'finish';
          end
        otherwise
          msg.point_type = '?HB?';
          msg.action     = 'warn';
          msg.wmsg       = 'could not determine type of bifurcation point';
      end
      msg.idx = 1;
    end
    
  case 'check'
    cdata       = coco_get_chart_data(cseg.curr_chart, fid);
    X           = data.ep_X;
    D           = cdata.la;
    [~, idx]    = min(abs(real(D)));
    la          = D(idx);
    X           = X(:,idx);
    if imag(la)<0
      la        = conj(la);
      X         = conj(X);
    end
    cdata.hb.la = la;
    cdata.hb.X  = X;
    cseg.curr_chart = coco_set_chart_data(cseg.curr_chart, fid, cdata);
    msg.action = 'add';
    msg.finish = true;
    
end

end

function [data, res] = bddat(prob, data, command, varargin) %#ok<INUSL>
%BDDAT   Append eigenvalues of Jacobian to BD.

res = {};
switch command
  case 'init'
    res = coco_get_id(data.oid,'eigs');
  case 'data'
    chart = varargin{1}; % Current chart
    cdata = coco_get_chart_data(chart, data.ep_tst.fid);
    if ~isempty(cdata) && isfield(cdata, 'la')
      res = cdata.la;
    end
end

end