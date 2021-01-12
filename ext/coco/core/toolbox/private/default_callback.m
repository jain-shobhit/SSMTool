function opts = default_callback(opts)
%COCO_DEFAULT_CALLBACK  Default handler for coco's call-back messages.
%
%   OPTS = COCO_DEFAULT_CALLBACK(OPTS) does nothing.
%
if nargout<1
	error('%s: too few output arguments', mfilename);
end
