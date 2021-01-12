function opts = default_print(opts, varargin)
%COCO_DEFAULT_PRINT  Default implementation for printing additional data.
%
%   OPTS = COCO_DEFAULT_PRINT(OPTS) does nothing.
%
if nargout<1
	error('%s: too few output arguments', mfilename);
end
