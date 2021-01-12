function prob = coco_add_signal(varargin)
%COCO_ADD_SIGNAL   Add signal to problem.
%
% PROB = COCO_ADD_SIGNAL(PROB, NAME, OWNER, OPTS...)
%
% OPTS = { '-max-slots' N }

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_add_signal.m 2839 2015-03-05 17:09:01Z fschild $

% varargin = { [PROB], NAME, OWNER, OPTIONS... }
s = coco_stream(varargin{:});

if isempty(s.peek) || isstruct(s.peek)
	prob = s.get;
else
	prob = [];
end

[name owner] = s.get;

sig     = struct('name', name, 'owner', owner, 'block', false);
sig     = parse_opts(sig, s);
signame = lower(name);

% check existing number of slots
if isfield(prob, 'slots') && isfield(prob.slots, signame)
  if sig.max_slots < numel(prob.slots.(signame))
    error('%s: attempt to add more than max_slots slot functions to signal ''%s''', ...
      mfilename, signame);
  end
end

if ~isfield(prob, 'signals')
  prob.signals = struct();
end

if isfield(prob.signals, signame)
  error('%s: signal ''%s'' already defined', mfilename, name);
else
  prob.signals.(signame) = sig;
end

end

function sig = parse_opts(sig, s)

sig.max_slots = inf;

o = lower(s.peek);
while ischar(o) && o(1)=='-'
  switch o
    case '-max-slots'
      s.skip;
      sig.max_slots = s.get;
    otherwise
      break
  end
  o = lower(s.peek);
end

end
