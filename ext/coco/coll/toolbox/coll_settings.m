function coll_settings(varargin)
%COLL_SETTINGS   Show and explain 'coll' instance settings of toolbox COLL.
%
% COLL_SETTINGS(VARARGIN)
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
% Show active settings of 'coll' instance appended to OID in PROB.
%
% VARARGIN = { PROB OID RUN [SOID] LAB }
% Show settings that will be active when restarting from a saved
% solution.
%
% The arguments and their meaning are identical to ODE_COLL2COLL.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% See also: ODE_COLL2COLL, ODE_SETTINGS

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_settings.m 3157 2019-12-16 23:36:58Z hdankowicz $

str = coco_stream(varargin{:});
if isstruct(str.peek)
  prob = str.get;
else
  prob = coco_prob;
end

switch numel(str)
  
  case 0
    coll = struct;
    tbid = 'coll';
%     [~,spec] = coll_get_settings(prob, tbid, coll);
    [sets,spec] = coll_get_settings(prob, tbid, coll);
    if nargin==0
      help coll_settings>main_settings
%       coco_explain_settings(prob, spec, coll, tbid);
      coco_explain_settings(spec, sets);
      fprintf('\n');
      ode_settings;
    else
%       coco_explain_settings(prob, spec, coll, tbid);
      coco_explain_settings(spec, sets);
      ode_settings(prob);
    end
    
  case 1
    coll = struct;
    oid  = str.get;
    tbid = coco_get_id(oid, 'coll');
%     [~,spec] = coll_get_settings(prob, tbid, coll);
%     coco_explain_settings(prob, spec, coll, tbid);
    [sets,spec] = coll_get_settings(prob, tbid, coll);
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
        [~,data] = coll_read_solution(args.soid, args.run, args.lab);
%         [~,spec] = coll_get_settings(prob, 'coll', data.coll);
%         coco_explain_settings(prob, spec, data.coll, 'coll');
        [sets,spec] = coll_get_settings(prob, 'coll', data.coll);
        coco_explain_settings(spec, sets);
        ode_settings(prob, '', data);
      catch read_error % test grammar = 'OID RUN LAB'
        try
          coll_settings(prob, args.run{:}, args.soid, args.run{:}, args.lab);
        catch %#ok<CTCH>
          rethrow(read_error)
        end
      end
    else
      grammar   = 'OID RUN [SOID] LAB';
      args_spec = {
        'OID',     '',   'str',  'oid',     '', 'read', {}
        'RUN', 'cell', '{str}',  'run',     {}, 'read', {}
        'SOID',     '',   'str', 'soid', '$OID', 'read', {}
        'LAB',     '',   'num',  'lab',     [], 'read', {}
        };
      opts_spec = {};
      args = coco_parse(grammar, args_spec, opts_spec, str);
      [~,data] = coll_read_solution(args.soid, args.run, args.lab);
      tbid = coco_get_id(args.oid, 'coll');
%       [~,spec] = coll_get_settings(prob, tbid, data.coll);
%       coco_explain_settings(prob, spec, data.coll, tbid);
      [sets,spec] = coll_get_settings(prob, tbid, data.coll);
      coco_explain_settings(spec, sets);
      ode_settings(prob, args.oid, data);
    end
    
  otherwise
    help coll_settings
    
end
end

function main_settings %#ok<DEFNU>
%Main settings of 'coll' instance of COLL toolbox:
%============================
%
%NTST : Number of mesh intervals used in trajectory discretization
%NCOL : Number of collocation nodes used in each mesh interval
%TOL  : Discretization error tolerance, defaults to 2/3rds power of the global
%       tolerance, which can be set to # by the command
%       prob = coco_set(prob, 'all', 'TOL', #)
%var  : Enable/disable temporary storage of fundamental solution to variational problem
%
%Full list of all 'coll' settings with default values:
%-------------------------------------------------
end
