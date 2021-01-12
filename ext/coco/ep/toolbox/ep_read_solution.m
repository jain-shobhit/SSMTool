function [sol, data] = ep_read_solution(oid, run, varargin)
%EP_READ_SOLUTION   Read solution and toolbox data from disk.
%
% [SOL DATA] = EP_READ_SOLUTION(VARARGIN)
%
% VARARGIN = { [OID] RUN LAB }
% Read solution data from solution data file of run RUN with solution label
% LAB.
%
% VARARGIN = { '' '' DATA }
% Construct solution data structure from initial solution guess for state
% vector and parameter values stored in fields with names 'x0' and 'p0' of
% the structure DATA.
%
% On input:
%
% OID : Optional object instance identifier (string, optional).
% RUN : Run identifier (string or cell-array of strings).
% LAB : Solution label (integer).
%
% DATA : Data structure with fields x0, p0 and pnames.
%
% On output:
%
% SOL  : Solution structure.
% DATA : Toolbox data structure or unchanged copy of input argument DATA.
%
% In the first calling form, EP_READ_SOLUTION reconstructs the solution and
% toolbox data structures of a saved equilibrium point and constructs
% restart information if the equilibrium point is a bifurcation point. More
% specifically, denote with
%
%   'ep'     :  a branch of equilibrium points,
%   'ep.VAR' :  a branch of equilibrium points with variational problem,
%   'ep.SN'  :  a branch of saddle-node bifurcation points,
%   'ep.HB'  :  a branch of Hopf bifurcation points,
%
% and with 'BR(SP)' a special point detected along a branch of equilibrium
% points of type BR, for example, with 'ep(HB)' a Hopf bifurcation point
% detected during equilibrium point continuation. DATA will always contain
% the fields of the ODE toolbox family. The solution structure SOL will
% have the fields
%
%   SOL.x    :  state vector
%   SOL.p    :  problem parameters
%   SOL.u    :  continuation variables u = [ x ; p ]
%   SOL.t    :  tangent vector, same size as u
%
% and additional fields encoding an initial solution point as required by
% EP_ADD. Depending on the types of the solution branch and the equilibrium
% point, the return value of SOL will have the following additional fields:
%
%   'ep(BP)' : For branch-points the field SOL.t0 will be initialized to
%      a singular vector normal to SOL.t.
%
%   'ep(SN)', 'ep(FP)' : For saddle-node bifurcations and fold points the
%      structures SOL.var and SOL.sn will be initialized with start data
%      for ep_add_var and ep_add_SN. Specifically, the field SOL.var.v is
%      set to a singular unit vector of the Jacobian of f.
%
%   'ep(HB)' : For Hopf bifurcation points the structures SOL.var and
%      SOL.hb will be initialized with start data for ep_add_var and
%      ep_add_HB. Specifically, the columns of the field SOL.var.v are set
%      to a pair of vectors v (of unit length) and J*v, in terms of the
%      Jacobian J of f, such that J*J*v=-om^2*v, where om is the Hopf
%      frequency.
%
%   'ep.var' : The structure SOL.var will be initialized with restart data
%      for ep_add_var.
%
%   'ep.SN'  : The structures SOL.var and SOL.sn will be initialized with
%      restart data for ep_add_var and ep_add_SN. In particular, the field
%      SOL.var.v is set to a singular unit vector of the Jacobian of f.
%
%   'ep.HB'  : The structures SOL.var and SOL.hb will be initialized with
%      start data for ep_add_var and ep_add_HB.  Specifically, the columns
%      of the field SOL.var.v are set to a pair of vectors v (of unit
%      length) and J*v, in terms of the Jacobian J of f, such that
%      J*J*v=-om^2*v, where om is the Hopf frequency.
%
% See also: COCO_READ_SOLUTION, ODE_ISOL2EP, ODE_EP2EP, EP_ADD, EP_ADD_SN,
% EP_ADD_HB, EP_ADD_VAR

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_read_solution.m 3037 2017-09-19 20:17:56Z hdankowicz $

if isempty(oid) && isempty(run)
  [sol, data] = read_sol_from(varargin{:});
  return
end

if nargin<3
  [oid, run, lab] = coco_deal('', oid, run);
else
  lab = varargin{1};
end

[tbid, format, branch_type] = guess_format(oid, run, lab);
[data, chart, uidx] = coco_read_solution(tbid, run, lab, 'data', ...
  'chart', 'uidx');

sol = struct('format', format, 'branch_type', branch_type, ...
  'pt_type', chart.pt_type, 'u', chart.x, 't', chart.t);

if coco_is_chart_data(chart, 'ep.test')
  sol.ep_test = coco_get_chart_data(chart, 'ep.test');
end

switch format
  
  case 'ep.v3'
    [sol, data] = read_ep_v3(sol, data, chart, uidx, tbid, branch_type, ...
      run, lab);

  case 'ep.v2' % Updated release from November 2014
    [sol, data] = read_ep_v2(sol, data, chart, uidx, tbid, branch_type, ...
      run, lab);
    
  case 'ep.v1' % Frank's original release
    [sol, data] = read_ep_v1(sol, data, chart, uidx, tbid, branch_type);
end

end

function [sol, data] = read_sol_from(data)
% Construct solution structure from data.

x0 = data.x0;
p0 = data.p0;

data.xdim = numel(x0);
data.pdim = numel(p0);

assert(data.pdim==numel(data.pnames) || isempty(data.pnames), ...
  '%s: incompatible number of elements for ''p0'' and ''pnames''', ...
  mfilename);

sol.x  = x0(:);
sol.p  = p0(:);
sol.u0 = [ sol.x ; sol.p ];
sol.t0 = [];

end

function [tbid, format, branch_type] = guess_format(oid, run, lab)
% Guess format of solution file from contents. This function will accept
% data files from other toolboxes, if all fields required by ep are
% present.

[tb, sol_info] = coco_read_tb_info(oid, run, lab, 'tb', 'ep');

if isempty(tb)
  tbid = coco_get_id(oid, 'ep');
else
  tbid = coco_get_id(oid, tb);
end

if ~isempty(sol_info)
  format = sol_info.format;
  branch_type = sol_info.branch_type;
  return
else
  format = 'ep.v1';
end

% Check for contents present in v1.
data = coco_read_solution(tbid, run, lab, 'data');
assert(~isempty(data), ...
  '%s: could not find solution data of toolbox instance ''%s''', ...
  mfilename, tbid);

if isfield(data, 'x_idx') ...
    && isfield(data, 'p_idx')
  branch_type = 'ep';
  
  if isfield(data, 'w') ...
      && isfield(data, 'v_idx') ...
      && isfield(data, 'k_idx')
    branch_type = 'ep.HB';
  elseif isfield(data, 'v_idx')
    branch_type = 'ep.SN';
  end
  
  return
end

error('%s: cannot restart from given solution data.', mfilename);

end

function [sol, data] = read_ep_v3(sol, data, chart, uidx, tbid, ...
  branch_type, run, lab)

sol.x  = sol.u(data.ep_eqn.x_idx);
sol.p  = sol.u(data.ep_eqn.p_idx);
sol.u0 = sol.u;
sol.t0 = [];

switch branch_type
  
  case 'ep'
    
    switch upper(sol.pt_type)
      
      case 'BP'
        cdata = coco_get_chart_data(chart, 'lsol');
        if isempty(cdata)
          coco_warn([], 1, 1, ...
            '%s: could not find restart data for branch-switching\n', ...
            mfilename);
        else
          sol.t0 = cdata.v(uidx);
        end
        
      case { 'SN' 'FP' }
        v          = chart.t(data.ep_eqn.x_idx);
        v          = v/norm(v);
        w          = 0*v;
        sol.var.v  = v;
        sol.var.w  = w;
        sol.var.u0 = [ v ; w ];
        sol.var.t0 = [];
        
        sol.sn.u0  = [];
        sol.sn.t0  = [];
        
      case 'HB'
        tfid  = coco_get_id(tbid, 'test');
        cdata = coco_get_chart_data(chart, tfid);
        if isempty(cdata)
          coco_warn([], 1, 1, ...
            '%s: could not find restart data for HB continuation\n', ...
            mfilename);
        else
          X  = cdata.hb.X;
          la = cdata.hb.la;
          v  = real(X);
          w  = imag(X);
          om = imag(la);
          k  = om^2;
          
          al = [-w'*v v'*v];
          al = al/norm(al);
          nv = al(1)*v + al(2)*w;
          nv = nv/norm(nv);
          
          sol.var.v  = [ v/norm(v) -sqrt(k)*w/norm(v) ];
          sol.var.w  = [ -sqrt(k)*w/norm(v) -k*v/norm(v) ];
          sol.var.u0 = [ sol.var.v(:) ; sol.var.w(:) ];
          sol.var.t0 = [];
          
          sol.hb.k  = k;
          sol.hb.nv = nv;
          sol.hb.u0 = k;
          sol.hb.t0 = [];
        end
    end
    
  case 'ep.VAR'
    fid = coco_get_id(tbid, 'var');
    chart = coco_read_solution(fid, run, lab, 'chart');
    v = chart.x(data.ep_var.v_idx);
    w = chart.x(data.ep_var.w_idx);
    sol.var.v  = v;
    sol.var.w  = w;
    sol.var.u  = [ v ; w ];
    sol.var.t  = chart.t([ data.ep_var.v_idx ; data.ep_var.w_idx ]);
    sol.var.u0 = [ v ; w ];
    sol.var.t0 = [];

  case 'ep.SN'
    fid = coco_get_id(tbid, 'var');
    chart = coco_read_solution(fid, run, lab, 'chart');
    v = chart.x(data.ep_var.v_idx);
    w = chart.x(data.ep_var.w_idx);
    sol.var.v  = v;
    sol.var.w  = w;
    sol.var.u  = [ v ; w ];
    sol.var.t  = chart.t([ data.ep_var.v_idx ; data.ep_var.w_idx ]);
    sol.var.u0 = [ v ; w ];
    sol.var.t0 = [];

    sol.sn.u  = [];
    sol.sn.t  = [];
    sol.sn.u0 = [];
    sol.sn.t0 = [];
    
  case 'ep.HB'
    fid = coco_get_id(tbid, 'var');
    chart = coco_read_solution(fid, run, lab, 'chart');
    v = chart.x(data.ep_var.v_idx);
    w = chart.x(data.ep_var.w_idx);
    t = chart.t([ data.ep_var.v_idx ; data.ep_var.w_idx ]);
    sol.var.v  = v;
    sol.var.w  = w;
    sol.var.u  = [ v(:) ; w(:) ];
    sol.var.t  = t(:);
    sol.var.u0 = [ v(:) ; w(:) ];
    sol.var.t0 = [];
    
    fid = coco_get_id(tbid, 'HB');
    chart = coco_read_solution(fid, run, lab, 'chart');
    v  = chart.x(data.ep_hb.v_idx);
    w  = chart.x(data.ep_hb.w_idx);
    k  = chart.x(data.ep_hb.k_idx);
    sol.hb.v  = v;
    sol.hb.w  = w;
    sol.hb.k  = k;
    sol.hb.u  = k;
    sol.hb.t  = chart.t(data.ep_hb.k_idx);
    sol.hb.u0 = k;
    sol.hb.t0 = [];
    sol.hb.nv = data.ep_hb.nv';
end

end

function [sol, data] = read_ep_v2(sol, data, chart, uidx, tbid, branch_type, run, lab)
sol.x  = sol.u(data.ep_eqn.x_idx);
sol.p  = sol.u(data.ep_eqn.p_idx);
sol.u0 = sol.u;
sol.t0 = [];

data.xdim = numel(sol.x);
data.pdim = numel(sol.p);

switch branch_type
  
  case 'ep'
    if data.ep.bifus
      tfid  = coco_get_id(tbid, 'test');
      tdata = coco_read_solution(tfid, run, lab, 'data');
      if ~isempty(tdata)
        data = tdata;
      end
    end
    
    switch upper(sol.pt_type)
      
      case 'BP'
        cdata = coco_get_chart_data(chart, 'lsol');
        if isempty(cdata)
          coco_warn([], 1, 1, ...
            '%s: could not find restart data for branch-switching\n', ...
            mfilename);
        else
          sol.t0 = cdata.v(uidx);
        end
        
      case { 'SN' 'FP' }
        v          = chart.t(data.ep_eqn.x_idx);
        v          = v/norm(v);
        w          = 0*v;
        sol.var.v  = v;
        sol.var.w  = w;
        sol.var.u0 = [ v ; w ];
        sol.var.t0 = [];
        sol.sn.u0  = [];
        sol.sn.t0  = [];
        
      case 'HB'
        tfid  = coco_get_id(tbid, 'test');
        cdata = coco_get_chart_data(chart, tfid);
        
        if isempty(cdata)
          coco_warn([], 1, 1, ...
            '%s: could not find restart data for HB continuation\n', ...
            mfilename);
        else
          X  = cdata.hb.X;
          la = cdata.hb.la;
          v  = real(X);
          w  = imag(X);
          om = imag(la);
          k  = om^2;
          
          al = [-w'*v v'*v];
          al = al/norm(al);
          nv = al(1)*v + al(2)*w;
          nv = nv/norm(nv);
          
          sol.var.v  = [ v/norm(v) -sqrt(k)*w/norm(v) ];
          sol.var.w  = [ -sqrt(k)*w/norm(v) -k*v/norm(v) ];
          sol.var.u0 = [ sol.var.v(:) ; sol.var.w(:) ];
          sol.var.t0 = [];
          
          sol.hb.k  = k;
          sol.hb.nv = nv;
          sol.hb.u0 = k;
          sol.hb.t0 = [];
          
        end
        
    end
    
  case 'ep.SN'
    fid = coco_get_id(tbid, 'var_sn');
    [data, chart] = coco_read_solution(fid, run, lab, 'data', 'chart');
    v = chart.x(data.ep_var.v_idx);
    w = chart.x(data.ep_var.w_idx);
    sol.var.v  = v;
    sol.var.w  = w;
    sol.var.u  = [ v ; w ];
    sol.var.t  = chart.t([ data.ep_var.v_idx ; data.ep_var.w_idx ]);
    sol.var.u0 = [ v ; w ];
    sol.var.t0 = [];
    
    fid  = coco_get_id(tbid, 'SN');
    data = coco_read_solution(fid, run, lab, 'data');
    sol.sn.u  = [];
    sol.sn.t  = [];
    sol.sn.u0 = [];
    sol.sn.t0 = [];
    
  case 'ep.HB'
    fid = coco_get_id(tbid, 'var_hb1');
    [data, chart] = coco_read_solution(fid, run, lab, 'data', 'chart');
    v = chart.x(data.ep_var.v_idx);
    w = chart.x(data.ep_var.w_idx);
    tv = chart.t(data.ep_var.v_idx);
    tw = chart.t(data.ep_var.w_idx);
    fid = coco_get_id(tbid, 'var_hb2');
    [data, chart] = coco_read_solution(fid, run, lab, 'data', 'chart');
    sol.var.v  = [ v chart.x(data.ep_var.v_idx) ];
    sol.var.w  = [ w chart.x(data.ep_var.w_idx) ];
    sol.var.u  = [ sol.var.v(:) ; sol.var.w(:) ];
    sol.var.t  = [ tv; chart.t(data.ep_var.v_idx); ...
      tw; chart.t(data.ep_var.w_idx) ];
    sol.var.u0 = [ sol.var.v(:) ; sol.var.w(:) ];
    sol.var.t0 = [];
    
    fid = coco_get_id(tbid, 'HB');
    [data, chart] = coco_read_solution(fid, run, lab, 'data', 'chart');
    v  = chart.x(data.ep_hb.v_idx);
    w  = chart.x(data.ep_hb.w_idx);
    k  = chart.x(data.ep_hb.k_idx);
    sol.hb.v  = v;
    sol.hb.w  = w;
    sol.hb.k  = k;
    sol.hb.u  = k;
    sol.hb.t  = chart.t(data.ep_hb.k_idx);
    sol.hb.u0 = k;
    sol.hb.t0 = [];
    sol.hb.nv = data.ep_hb.nv';
end
end

function [sol, data] = read_ep_v1(sol, data, chart, uidx, tbid, branch_type)
sol.x = sol.u(data.x_idx);
sol.p = sol.u(data.p_idx);

data.xdim = numel(sol.x);
data.pdim = numel(sol.p);

switch branch_type
  case 'ep'
    
    switch upper(sol.pt_type)
      
      case 'BP'
        cdata  = coco_get_chart_data(chart, 'lsol');
        if isempty(cdata)
          coco_warn([], 1, 1, ...
            '%s: could not find restart data for branch-switching\n', ...
            mfilename);
        else
          sol.t0 = cdata.v(uidx);
        end
        
      case { 'SN' 'FP' }
        v          = chart.t(data.x_idx);
        v          = v/norm(v);
        w          = 0*v;
        sol.var.v  = v;
        sol.var.w  = w;
        sol.var.u0 = [ v ; w ];
        sol.var.t0 = [];
        sol.sn.u0  = [];
        sol.sn.t0  = [];
        
      case 'HB'
        tfid  = coco_get_id(tbid, 'test');
        cdata = coco_get_chart_data(chart, tfid);
        
        if isempty(cdata)
          coco_warn([], 1, 1, ...
            '%s: could not find restart data for HB continuation\n', ...
            mfilename);
        else
          v  = cdata.v;
          w  = cdata.w;
          om = cdata.om;
          k  = om^2;
          
          al = [-w'*v v'*v];
          al = al/norm(al);
          nv = al(1)*v + al(2)*w;
          nv = nv/norm(nv);
          
          sol.var.v  = [ v/norm(v) -sqrt(k)*w/norm(v) ];
          sol.var.w  = [ -sqrt(k)*w/norm(v) -k*v/norm(v) ];
          sol.var.u0 = [ sol.var.v(:) ; sol.var.w(:) ];
          sol.var.t0 = [];
          
          sol.hb.k  = k;
          sol.hb.nv = nv;
          sol.hb.u0 = k;
          sol.hb.t0 = [];
        end
    end
    
  case 'ep.SN'
    v          = sol.u(data.v_idx);
    w          = 0*v;
    sol.var.v  = v;
    sol.var.w  = w;
    sol.var.u  = [ v ; w ];
    sol.var.t  = [ sol.t(data.v_idx) ; 0*sol.t(data.v_idx) ];
    sol.var.u0 = [ v ; w ];
    sol.var.t0 = [];
    
    sol.sn.u  = [];
    sol.sn.t  = [];
    sol.sn.u0 = [];
    sol.sn.t0 = [];
    
  case 'ep.HB'
    Jx  = ep_fhan_DFDX(data, sol.x, sol.p);
    v1  = sol.u(data.v_idx);
    v2  = Jx*v1;
    v3  = Jx*v2;
    tv1 = sol.t(data.v_idx);
    tv2 = Jx*tv1;
    tv3 = Jx*tv2;
    k   = sol.u(data.k_idx);
    
    sol.var.v  = [ v1 v2 ];
    sol.var.w  = [ v2 v3 ];
    sol.var.u  = [ sol.var.v(:) ; sol.var.w(:) ];
    sol.var.t  = [ tv1 ; tv2 ; tv2; tv3];
    sol.var.u0 = [ sol.var.v(:) ; sol.var.w(:) ];
    sol.var.t0 = [];
       
    sol.hb.v  = v1;
    sol.hb.w  = v3;
    sol.hb.k  = k;
    sol.hb.u  = k;
    sol.hb.t  = sol.t(data.k_idx);
    sol.hb.u0 = k;
    sol.hb.t0 = [];
    sol.hb.nv = data.w';
end

% remove surplus components from sol.u and sol.t
uidx   = [data.x_idx ; data.p_idx];
sol.u  = sol.u(uidx);
sol.t  = sol.t(uidx);
sol.u0 = sol.u;
sol.t0 = [];
end

function Jx = ep_fhan_DFDX(data, x, p)
%EP_FHAN_DFDX   Compute Jacobian of vector field w.r.t state vector.

if isempty(data.dfdxhan)
  Jx = coco_ezDFDX('f(x,p)', data.fhan, x, p);
else
  Jx = data.dfdxhan(x, p);
end

end
