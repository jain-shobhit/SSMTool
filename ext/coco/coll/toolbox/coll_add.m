function [prob, data] = coll_add(prob, data, sol, varargin)
%COLL_ADD   Add 'coll' instance trajectory segment zero problem.
%
% [PROB DATA] = COLL_ADD(PROB, DATA, SOL, [OPTS])
% OPTS = { '-cache-jac' | '-no-err' | '-no-var' | '-no-bddat' | '-xid' NAME}
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function adds a 'coll' instance consisting of
%
%   - a trajectory segment zero problem with identifier
%       coco_get_id(DATA.oid, 'coll'), 
%   - the corresponding error estimation test function with identifier
%       coco_get_id(DATA.oid, 'coll.err'),
%   - the corresponding variational problem test function with identifier
%       coco_get_id(DATA.oid, 'coll.test'), and
%   - a slot function with identifier coco_get_id(DATA.oid, 'coll'), adding
%       the approximate L2 norm of the state vector time history to the
%       bifurcation data returned by the coco entry-point function or
%       stored to disk.
%
% The corresponding function data is stored in the fields DATA.coll_seg
% (zero problem and error estimator), DATA.coll_tst, and DATA.coll_bddat,
% respectively. If the option '-cache-jac' is passed, the read-only field
% DATA.coll_Jx contains the vector field Jacobian corresponding to the most
% recent evaluation of the trajectory segment zero problem.
%
% On input:
%
% PROB : Continuation problem structure.
% DATA : An initialized 'coll' toolbox instance data structure (see also
%        ODE_INIT_DATA).
% SOL  : 'coll' instance solution structure (see also COLL_READ_SOLUTION).
%        The fields SOL.tbp, SOL.xbp, SOL.T, and SOL.p must contain a time
%        discretization and corresponding sampled time history along a
%        trajectory, the interval length, and the problem parameter values,
%        respectively. If the SOL.t0 field is nonempty, it must be
%        accompanied by the fields SOL.xbp_t0, SOL.T_t0, and SOL.p_t0,
%        which determine the initial continuation direction along the
%        solution manifold prior to remeshing. Typically, SOL.t0=[], which
%        means that the continuation direction is determined by the atlas
%        algorithm. A non-empty vector SOL.t0 must be given for
%        branch-switching at branch points.
% OPTS : '-cache-jac', '-no-err', '-no-var', '-no-bddat' and '-xid' NAME
%        (optional, multiple options may be given). '-cache-jac' enables
%        cacheing of the vector field Jacobian, which is required when
%        embedding the variational problem. '-no-err' disables adding of
%        error test functions and events, and '-no-var' disables adding of
%        the variational problem test function. '-no-bddat' disables the
%        addition of the approximate L2 norm of the state vector time
%        history to the bifurcation data returned by the coco entry-point
%        function or stored to disk. The option '-xid' NAME assigns
%        sprintf('||%s||', coco_get_id(DATA.oid, NAME)) as the label of the
%        bifurcation data column containing the approximate L2 norm. NAME
%        defaults to 'x'.
%
% On output:
%
% PROB : Continuation problem structure with a 'coll' instance added.
% DATA : Toolbox data structure with 'coll' trajectiry segent zero problem
%        data added. COLL_ADD will add up to three structures to DATA:
%        coll_seg containing additional function data of the 'coll'
%        trajectory segment zero problem, coll_tst (if the variational
%        problem test function is added) containing additional function
%        data of the corresponding variational problem test function, and
%        coll_bddat (if COLL bifurcation data is added) containing COLL's
%        slot function data for the slot 'bddat'. These structures contains
%        the following fields:
%
%        coll_seg.fid    : function identifier of 'coll' trajectory segment
%           zero problem; fid = coco_get_id(OID, 'coll')
%        coll_seg.int    : context-independent properties of the problem
%           formulation that are independent of the discretization order
%           and the discretization parameters
%        coll_seg.maps   : context-independent properties of the problem
%           formulation that depend on the discretization order but not on
%           the discretization parameters
%        coll_seg.mesh   : context-independent properties of the problem
%           formulation that depend on the discretization parameters
%
%        coll_tst.fid    : function identifier of 'coll' variational
%           problem test function; fid = coco_get_id(OID, 'coll.test')
%        coll_tst.M      : fundamental solution to variational problem; not
%           saved with solution data
%        coll_tst.M0_idx : context-independent indices of rows in solution
%           to variational problem associated with the initial condition;
%           not saved with solution data
%        coll_tst.M1_idx : context-independent indices of rows in solution
%           to variational problem associated with the final condition; not
%           saved with solution data
%        coll_tst.row, coll_tst.row0, coll_tst.rhs, coll_tst.rh0,
%        coll_tst.M0, coll_st.B1, coll_st.B2, coll_tst.la1, and
%        coll_tst.la2: parameters that describe the variational boundary
%           conditions and their updates; not saved with solution data
%
%        coll_bddat.xid  : name set with option '-xid NAME'; defaults to
%            'x'
%
% See also: ODE_INIT_DATA, COLL_READ_SOLUTION, ODE_ISOL2COLL,
% ODE_COLL2COLL, COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_add.m 2998 2017-01-29 16:25:42Z hdankowicz $

opts_spec = {
  '-cache-jac',  'CJ', false, 'toggle', {}
     '-no-err', 'ERR', false, 'toggle', {}
     '-no-var', 'VAR', false, 'toggle', {}
   '-no-bddat', 'DAT', false, 'toggle', {}
        '-xid', 'xid',   'x',   'read', {}
  };

opts = coco_parse_opts(opts_spec, varargin{:});
tbid = coco_get_id(data.oid, 'coll');
data.coll = coll_get_settings(prob, tbid, data.coll);
[prob, data] = coll_construct_seg(prob, data, sol, opts.CJ);
if ~opts.ERR
  [prob, data] = coll_construct_err(prob, data);
end
if data.coll.var && ~opts.VAR
  [prob, data] = coll_construct_var(prob, data);
end
if ~opts.DAT
  data.coll_bddat.xid = coco_get_id(data.oid, opts.xid);
  prob = coco_add_slot(prob, tbid, @bddat, data, 'bddat');
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

function [data, res] = bddat(prob, data, command, varargin)
%BDDAT   Append COLL bifurcation data to BD.

res = {};
switch command
  
  case 'init'
    xid   = data.coll_bddat.xid;
    nmid1 = sprintf('||%s||_{L_2[0,1]}', xid);
    nmid2 = sprintf('||%s||_{L_2[0,T]}', xid);
    ntsid = coco_get_id(data.oid, 'NTST');
    res   = { nmid1, nmid2 , ntsid };
    
  case 'data'
    seg   = data.coll_seg;
    maps  = seg.maps;
    mesh  = seg.mesh;
    chart = varargin{1}; % Current chart
    uidx  = coco_get_func_data(prob, seg.fid, 'uidx');
    u     = chart.x(uidx);
    x     = u(maps.xbp_idx);
    xcn   = reshape(maps.W*x, maps.x_shp);
    NTST  = maps.NTST;
    wts   = mesh.gwt; % Gauss weights
    kas   = mesh.gka; % Warping coefficients
    T     = u(maps.T_idx);
    nrmx  = sqrt((0.5/NTST)*sum(xcn.*xcn,1)*(wts.*kas)');
    res   = { nrmx , sqrt(T)*nrmx, NTST };
    
end

end
