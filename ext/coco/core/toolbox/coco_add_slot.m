function prob = coco_add_slot(varargin)
%COCO_ADD_SLOT   Add slot function to a signal.
%
% PROB = COCO_ADD_SLOT([PROB], FID, @FUNC, DATA, SIGNAL, [@COPY])
%
% ---------------------------------------------------------
% [DATA OUT1 ... OUTN] = FUNC(PROB, DATA, ARG1, ... , ARGM)
%
% -----------------------
% DATA = COPY(OPTS, DATA)

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_add_slot.m 2839 2015-03-05 17:09:01Z fschild $

%% parse input arguments
%  varargin = { [OPTS], FID, @FUNC, DATA, SIGNAL, [@COPY] }

argidx = 1;

if isempty(varargin{argidx}) || isstruct(varargin{argidx})
	prob   = varargin{argidx};
	argidx = argidx + 1;
else
	prob = [];
end

fid     = varargin{argidx  };
coco_opts_tree.check_path(fid);
fhan    = varargin{argidx+1};
data    = varargin{argidx+2};
list    = varargin{argidx+3};
signame = lower(list);
if strcmp('funcs', signame)
  error('%s: illegal slot name ''%s''', mfilename, list);
end
if nargin>argidx+3
  copy = varargin{argidx+4};
else
  copy = [];
end

%% create slot structure or check for duplicate identifyer
if ~isfield(prob, 'slots')
  prob.slots.funcs = struct([]);
  if ~isfield(prob, 'signals')
    prob.signals = struct();
  end
end

if isfield(prob.slots, signame)
  fidx = prob.slots.(signame);
  fids = { prob.slots.funcs(fidx).identifyer };
  if any( strcmp(fid, fids) )
    error('%s: slot identifyer ''%s'' already in use at signal ''%s''', ...
      mfilename, fid, list);
  end
else
  prob.slots.(signame) = [];
end

%% check max number of slots
if isfield(prob, 'signals') && isfield(prob.signals, signame)
  if prob.signals.(signame).max_slots <= numel(prob.slots.(signame))
    error('%s: attempt to add more than max_slots slot functions to signal ''%s''', ...
      mfilename, signame);
  end
end

%% add function to list

func.identifyer      = fid;
func.F               = fhan;
func.data            = data;
func.copy            = copy;
prob.slots.funcs     = [ prob.slots.funcs func ];
idx                  = numel(prob.slots.funcs);
prob.slots.(signame) = [ prob.slots.(signame) idx ];
end
