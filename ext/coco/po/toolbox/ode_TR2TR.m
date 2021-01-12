function prob = ode_TR2TR(prob, oid, varargin)
%ode_TR2TR   Start continuation of Neimark-Sacker bifurcation points.
%
% PROB = ode_TR2TR(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB ... }
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% This function determines the type of solution stored in the solution file
% with label LAB of run RUN and forwards the call to the appropriate
% toolbox function to start a continuation of Neimark-Sacker bifurcation
% points. See the documentation of toolboxes of the toolbox family ODE for
% information about which toolboxes support Neimark-Sacker continuation.
%
% See also: ODE_PO2TR

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_SN2SN.m 2839 2015-03-05 17:09:01Z fschild $

grammar   = 'RUN [SOID] LAB';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {};
str  = coco_stream(varargin{:});
args = coco_parse(grammar, args_spec, opts_spec, str);

sol_type = coco_read_tb_info(args.soid, args.run, args.lab, 'sol_type');
if isempty(sol_type)
  sol_type = guess_sol_type(args.soid, args.run, args.lab);
end

assert(~isempty(sol_type), ...
  '%s: could not determine type of saved solution', mfilename);

ctor_nm = sprintf('ode_%s2TR', sol_type);
if any(exist(ctor_nm, 'file') == [2 3 6])
  ctor = str2func(ctor_nm);
  prob = ctor(prob, oid, str.put(args.run, args.soid, args.lab));
else
  error('%s: could not find TR constructor for solution type ''%s''', ...
  	mfilename, sol_type);
end

end

function sol_type = guess_sol_type(oid, run, lab)
  tbid = coco_get_id(oid, 'po');
  data = coco_read_solution(tbid, run, lab, 'data');
  if ~isempty(data)
    sol_type = 'po';
    return
  end
  run_data = coco_read_solution(run, lab, 'run');
  switch run_data.sol_type
    case 'TR'
      sol_type = '';
    otherwise
      sol_type = run_data.sol_type;
  end
end
