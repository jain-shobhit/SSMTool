function [opts varargout] = coco_emit(opts, list, varargin)
%coco_emit   Emit signal to slot list list.
%
%   [OPTS ...] = COCO_emit(OPTS, LIST, ...)
%

signame = lower(list);

if ~isfield(opts, 'signals') || ~isfield(opts.signals, signame)
  error('%s: attempt to emit undefined signal ''%s''', mfilename, list);
end

if opts.signals.(signame).block
  error('%s: attempt to emit signal ''%s'' recursively', mfilename, list);
end

opts.signals.(signame).block = true;

if isfield(opts, 'slots') && isfield(opts.slots, signame)
  slist = opts.slots.(signame);
  funcs = opts.slots.funcs;
  out   = {};
  lout  = cell(1,nargout-1);
  for i=slist
    lout{1} = funcs(i).identifyer;
    data    = funcs(i).data;
    [data lout{2:end}] = funcs(i).F(opts, data, varargin{:});
    opts.slots.funcs(i).data = data;
    out = [ out ; lout ]; %#ok<AGROW>
  end
  for i=1:nargout-1
    varargout{i} = out(:,i);
  end
else
  [varargout{1:nargout-1}] = deal({});
end

opts.signals.(signame).block = false;

end
