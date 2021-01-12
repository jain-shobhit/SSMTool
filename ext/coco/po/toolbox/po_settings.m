function po_settings(varargin)
%PO_SETTINGS   Show and explain 'po' instance settings of toolbox PO.
%
% PO_SETTINGS(VARARGIN)
%
% VARARGIN = { }
% Show and explain default settings.
%
% VARARGIN = { RUN [SOID] LAB }
% Show settings used to compute a saved solution.
%
% VARARGIN = { PROB }
% Show active settings.
%
% VARARGIN = { PROB OID }
% Show active settings of 'po' instance appended to OID in PROB.
%
% VARARGIN = { PROB OID RUN [SOID] LAB }
% Show settings that will be active when restarting from a saved
% solution.
%
% The arguments and their meaning are identical to ODE_PO2PO.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% See also: ODE_PO2PO, ODE_SETTINGS

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: PO_settings.m 2839 2015-03-05 17:09:01Z fschild $

str = coco_stream(varargin{:});
if isstruct(str.peek())
  prob = str.get();
else
  prob = coco_prob();
end

switch numel(str)
  
  case 0
    po   = struct();
    tbid = 'po';
%     [~,spec] = po_get_settings(prob, tbid, po);
    [sets,spec] = po_get_settings(prob, tbid, po);
    if nargin==0
      help po_settings>main_settings
%       coco_explain_settings(prob, spec, po, tbid);
      coco_explain_settings(spec, sets);
      fprintf('\n');
      ode_settings();
    else
%       coco_explain_settings(prob, spec, po, tbid);
      coco_explain_settings(spec, sets);
      ode_settings(prob);
    end
    
  case 1
    po   = struct();
    oid  = str.get();
    tbid = coco_get_id(oid, 'po');
%     [~,spec] = po_get_settings(prob, tbid, po);
%     coco_explain_settings(prob, spec, po, tbid);
    [sets,spec] = po_get_settings(prob, tbid, po);
    coco_explain_settings(spec, sets);
    ode_settings(prob, oid);
    
  case {2,3,4}
    if nargin<4
      grammar   = 'RUN [SOID] LAB';
      args_spec = {
        'RUN', 'cell', '{str}',  'run', {}, 'read', {}
        'SOID',     '',   'str', 'soid', '', 'read', {}
        'LAB',     '',   'num',  'lab', [], 'read', {}
        };
      opts_spec = {};
      args = coco_parse(grammar, args_spec, opts_spec, str);
      try % Bug: no rigorous distinction between 'RUN SOID LAB' and 'OID RUN LAB'. Priority given to the former.
        [~,data] = po_read_solution(args.soid, args.run, args.lab);
%         [~,spec] = po_get_settings(prob, 'po', data.po);
%         coco_explain_settings(prob, spec, data.po, 'po');
        [sets,spec] = po_get_settings(prob, 'po', data.po);
        coco_explain_settings(spec, sets);
        ode_settings(prob, '', data);
      catch read_error % test grammar = 'OID RUN LAB'
        try
          po_settings(prob, args.run{:}, args.soid, args.run{:}, args.lab);
        catch %#ok<CTCH>
          rethrow(read_error)
        end
      end
    else
      grammar   = 'OID RUN SOID LAB';
      args_spec = {
        'OID',     '',   'str',  'oid',     '', 'read', {}
        'RUN', 'cell', '{str}',  'run',     {}, 'read', {}
        'SOID',     '',   'str', 'soid',    '', 'read', {}
        'LAB',     '',   'num',  'lab',     [], 'read', {}
        };
      opts_spec = {};
      args = coco_parse(grammar, args_spec, opts_spec, str);
      [~,data] = po_read_solution(args.soid, args.run, args.lab);
      tbid = coco_get_id(args.oid, 'po');
%       [~,spec] = po_get_settings(prob, tbid, data.po);
%       coco_explain_settings(prob, spec, data.po, tbid);
      [sets,spec] = po_get_settings(prob, tbid, data.po);
      coco_explain_settings(spec, sets);
      ode_settings(prob, args.oid, data);
    end
    
  otherwise
    help po_settings
    
end
end

function main_settings %#ok<DEFNU>
%Main settings of 'po' instance of 'po' toolbox:
%============================
%
%bifus : Enable/disable detection of bifurcations. By default, the PO
%  toolbox will attempt to locate all bifurcation points for which a test
%  function is implemented. Set to off/false to disable detection of all
%  bifurcation points. It is possible to disable the detection of
%  particular bifurcation points; see PO_SETTINGS for details.
%
%  Note that this setting for the PO toolbox will not affect bifurcation
%  detection implemented by other toolboxes, for example, the atlas
%  classes.
%
%Full list of all 'po' instance settings with default values:
%-------------------------------------------------
end
