function ep_settings(varargin)
%EP_SETTINGS   Show and explain settings of 'ep' toolbox.
%
% EP_SETTINGS(VARARGIN)
%
% VARARGIN = { }
% Show and explain default settings of EP.
%
% VARARGIN = { RUN [SOID] LAB }
% Show settings of EP used to compute a saved solution.
%
% VARARGIN = { PROB }
% Show active settings of EP.
%
% VARARGIN = { PROB OID }
% Show active settings of instance of EP appended to OID in PROB.
%
% VARARGIN = { PROB OID RUN [SOID] LAB }
% Show settings of EP that will be active when restarting from a saved
% solution.
%
% The arguments and their meaning are identical to ODE_EP2EP.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% See also: ODE_EP2EP, ODE_SETTINGS

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_settings.m 3157 2019-12-16 23:36:58Z hdankowicz $

str = coco_stream(varargin{:});
if isstruct(str.peek)
  prob = str.get;
else
  prob = coco_prob;
end

switch numel(str)
  
  case 0
    ep   = struct;
    tbid = 'ep';
%     [~,spec] = ep_get_settings(prob, tbid, ep);
    [sets,spec] = ep_get_settings(prob, tbid, ep);
    if nargin==0
      help ep_settings>main_settings
%       coco_explain_settings(prob, spec, ep, tbid);
      coco_explain_settings(spec, sets);
      fprintf('\n');
      ode_settings;
    else
%       coco_explain_settings(prob, spec, ep, tbid);
      coco_explain_settings(spec, sets);
      ode_settings(prob);
    end
    
  case 1
    ep   = struct;
    oid  = str.get;
    tbid = coco_get_id(oid, 'ep');
%     [~,spec] = ep_get_settings(prob, tbid, ep);
    [sets,spec] = ep_get_settings(prob, tbid, ep);
%     coco_explain_settings(prob, spec, ep, tbid);
    coco_explain_settings(spec, sets);
    ode_settings(prob, oid);
    
  case {2,3,4}
    if nargin<4
      grammar   = 'RUN [SOID] LAB';
      args_spec = {
        'RUN', 'cell', '{str}',  'run', {}, 'read', {}
        'SOID',    '',   'str', 'soid', '', 'read', {}
        'LAB',     '',   'num',  'lab', [], 'read', {}
        };
      opts_spec = {};
      args = coco_parse(grammar, args_spec, opts_spec, str);
      try % Bug: no rigorous distinction between 'RUN SOID LAB' and 'OID RUN LAB'. Priority given to the former.
        [~,data] = ep_read_solution(args.soid, args.run, args.lab);
%         [~,spec] = ep_get_settings(prob, 'ep', data.ep);
%         coco_explain_settings(prob, spec, data.ep, 'ep');
        [sets,spec] = ep_get_settings(prob, 'ep', data.ep);
        coco_explain_settings(spec, sets);
        ode_settings(prob, '', data);
      catch read_error % test grammar = 'OID RUN LAB'
        try
          ep_settings(prob, args.run{:}, args.soid, args.run{:}, args.lab);
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
      [~,data] = ep_read_solution(args.soid, args.run, args.lab);
      tbid = coco_get_id(args.oid, 'ep');
%       [~,spec] = ep_get_settings(prob, tbid, data.ep);
%       coco_explain_settings(prob, spec, data.ep, tbid);
      [sets,spec] = ep_get_settings(prob, tbid, data.ep);
      coco_explain_settings(spec, sets);
      ode_settings(prob, args.oid, data);
    end
    
  otherwise
    help ep_settings
    
end
end

function main_settings %#ok<DEFNU>
%Main settings of 'ep' toolbox:
%============================
%
%bifus : Enable/disable detection of bifurcations. By default, the EP
%  toolbox will attempt to locate all bifurcation points for which a test
%  function is implemented. Set to off/false to disable detection of all
%  bifurcation points. It is possible to disable the detection of
%  particular bifurcation points; see EP_SETTINGS for details.
%
%  Note that this setting for the EP toolbox will not affect bifurcation
%  detection implemented by other toolboxes, for example, the atlas
%  classes.
%
%Full list of all 'ep' toolbox settings with default values:
%-------------------------------------------------
end
