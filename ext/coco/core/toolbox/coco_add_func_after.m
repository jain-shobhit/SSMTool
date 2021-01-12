function prob = coco_add_func_after(varargin)
%COCO_ADD_FUNC_AFTER   Add zero problem after adding other(s) first.
%
% PROB = COCO_ADD_FUNC_AFTER([PROB], DEPLIST, @CTOR, ARG1, ..., ARGN)
%
% --------------------------
% PROB = CTOR(PROB, ARGS{:})

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_add_func_after.m 2839 2015-03-05 17:09:01Z fschild $

%% parse input arguments
%  varargin = { [PROB], DEPLIST, @CTOR, ... }

argidx = 1;

if isempty(varargin{argidx}) || isstruct(varargin{argidx})
	prob   = varargin{argidx};
	argidx = argidx + 1;
else
	prob = [];
end

if ~isfield(prob, 'efunc')
  prob.efunc = efunc_new([]);
end
if ~isfield(prob.efunc, 'identifyers')
  prob.efunc = efunc_new(prob.efunc);
end

if iscell(varargin{argidx})
  deplist = varargin{argidx};
elseif ischar(varargin{argidx})
  deplist = { varargin{argidx} };
else
  deplist = varargin{argidx};
end
argidx = argidx + 1;

ctor = varargin{argidx};
if ~isa(ctor, 'function_handle')
  error('%s: argument %d must be a function handle', mfilename, argidx);
end
argidx = argidx + 1;

args = { varargin{argidx:end} };

% add function to pending list
prob.efunc.pending = [prob.efunc.pending ; { deplist ctor args } ];

% try to add pending functions
if prob.efunc.add_pending
  prob = coco_add_pending(prob);
end

end
