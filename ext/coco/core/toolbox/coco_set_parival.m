function opts = coco_set_parival(varargin)
%COCO_SET_PARIVAL
%
%   OPTS = COCO_SET_PARIVAL([OPTS], PARS, VALS) defines initial values for
%   parameters specified in PARS. PARS is a parameter name or a cell array
%   of parameter names and VALS is an array or cell array of corresponding
%   values. This function is useful to assign a specific initial value to
%   an inactive monitor function. If the value assigned to a parameter is
%   empty, the initial value for the corresponding parameter is removed.
%
%   OPTS = COCO_SET_PARIVAL(OPTS, []) deletes all initial values.
%
%   NOTE: Initial values should always be cleared after a continuation run.

%% affected fields in opts
%
%    opts.efunc.parivals - cell array containing a list of pairs of
%                          parameter names and values to be assigned as
%                          initial values in efunc_init

%% check for input argument opts
%  varargin = { [opts], PARS, [VALS] }

argidx = 1;
if isempty(varargin{argidx}) || isstruct(varargin{argidx})
	opts   = varargin{argidx};
	argidx = argidx + 1;
else
	opts = [];
end

%% parse input argument pars

pars = varargin{argidx};

if isempty(pars)
	opts.efunc.parivals = {};
	return;
elseif ischar(pars)
	pars = { pars };
end

%% parse input argument vals

vals = varargin{argidx+1};

%% initialise opts.efunc.parivals if necessary

if ~isfield(opts, 'efunc')
	opts.efunc.parivals = {};
end

if ~isfield(opts.efunc, 'parivals')
	opts.efunc.parivals = {};
end

%% add entries to initial values list

for i=1:numel(pars)
  row = [];
  if ~isempty(opts.efunc.parivals)
    row = find(strcmp(pars{i}, { opts.efunc.parivals{:,1} }));
  end
	if isempty(row)
		row = size(opts.efunc.parivals,1) + 1;
		opts.efunc.parivals{row,1} = pars{i};
	end
	if isempty(vals)
		opts.efunc.parivals(row,:) = [];
	else
		if iscell(vals)
			val = vals{i};
		else
			val = vals(i);
		end
		if isempty(val)
			opts.efunc.parivals{row,2} = []; % this syntax works only this way
		else
			opts.efunc.parivals{row,2} = val;
		end
	end
end
