function [prob, data] = po_add(prob, data, varargin)
%PO_ADD   Add 'po' instance periodic orbit zero problem.
%
% [PROB DATA] = PO_ADD(PROB, DATA, [OPTS])
% OPTS = { '-no-phase' | '-no-test'}
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function adds a 'po' instance consisting of
%
%   - a periodic orbit zero problem with identifier
%     coco_get_id(DATA.oid, 'po'),
%   - the corresponding codim-1 test functions with identifier
%     coco_get_id(DATA.oid, 'po.test'), and
%   - the orbital period monitor function with identifier
%     coco_get_id(DATA.oid, 'po.period').
%
% The corresponding function data is stored in the fields DATA.po_orb and
% DATA.po_tst, respectively.
%
% On input:
%
% PROB : Continuation problem structure. 
% DATA : An initialized 'po' toolbox instance data structure (see also
%        PO_INIT_DATA). 
%
% OPTS : '-no-phase' and'-no-test' (optional, multiple options may be
%        given). '-no-phase' disables adding an integral phase condition to
%        the periodic orbit zero problem and '-no-test' disables adding of
%        codim-1 test functions.
%
% On output:
%
% PROB : Continuation problem structure with a 'po' instance added.
% DATA : Toolbox data structure with 'po' periodic orbit zero problem data
%        added. PO_ADD will add up to two structures to DATA: po_orb
%        containing additional function data of the 'po' periodic orbit
%        zero problem and po_tst (if test functions are added) containing
%        additional function data of the corresponding codim-1 test
%        functions. These structures contains the following fields:
%
%        po_orb.dim : number of state space dimensions
%        po_orb.fid : function identifier of 'po' periodic orbit zero
%            problem; fid = coco_get_id(OID, 'po')
%        po_orb.J   : coefficient matrix of periodic orbit zero problem;
%            not saved with solution data
%        po_orb.intfac : integration factor in integral phase condition for
%            autonomous vector field; not saved with solution data
%
%        po_tst.la_idx1 and po_tst.la_idx2 : eigenvalue expansion vectors
%            required by the Hopf test function; only non-empty if Hopf
%            detection is enabled; require n^2 double values storage space;
%            not saved with solution data
%        po_tst.fid : function identifier of 'po' codim-1 test functions;
%            fid = coco_get_id(OID, 'po.test')
%
% See also: PO_INIT_DATA, PO_READ_SOLUTION, ODE_ISOL2PO, ODE_PO2PO,
% COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_add.m 2948 2016-04-10 10:41:39Z fschild $

opts_spec = {
  '-no-phase', 'PHS', false, 'toggle', {}
   '-no-test',  'NT', false, 'toggle', {}
  '-no-bddat', 'DAT', false, 'toggle', {}
       '-xid', 'xid',   'x',   'read', {}
  };
opts = coco_parse_opts(opts_spec, varargin{:});

[prob, data] = po_close_orb(prob, data, ~opts.PHS);
if ~opts.NT && data.po.bifus
  [prob, data] = po_construct_tst(prob, data);
end
tbid = coco_get_id(data.oid, 'po');
if ~opts.DAT
  data.po_bddat.xid = coco_get_id(data.oid, opts.xid);
  prob = coco_add_slot(prob, tbid, @bddat, data, 'bddat');
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

function [data, res] = bddat(prob, data, command, varargin)
%BDDAT   Append PO bifurcation data to BD.

res = {};
switch command
  
  case 'init'
    xid   = data.po_bddat.xid;
    nrmid = sprintf('||%s||_{2,MPD}', xid); % mean + deviation plotting measure
    maxid = sprintf('MAX(%s)', xid);
    minid = sprintf('MIN(%s)', xid);
    res   = { nrmid, maxid, minid };
    
  case 'data'
    chart = varargin{1};
    [fdata uidx] = coco_get_func_data(prob, data.cid, 'data', 'uidx');
    maps = fdata.coll_seg.maps;
    mesh = fdata.coll_seg.mesh;
    xbp  = reshape(chart.x(uidx(maps.xbp_idx)), maps.xbp_shp);
    iw   = (0.5/maps.NTST)*mesh.kas2*mesh.wts2*maps.W;
    xav  = sum(reshape(iw*xbp(:), maps.x_shp),2);
    dev  = bsxfun(@minus, xbp, xav);
    ndv  = reshape(iw*(dev(:).^2), maps.x_shp);
    mpd  = norm(xav)+sqrt(sum(sum(ndv,2),1));
    res  = { mpd , max(xbp,[],2) , min(xbp,[],2) };
    
end

end
