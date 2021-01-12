function opts = default_update(opts)
%COCO_DEFAULT_UPDATE  Default handler for coco's update message.
%
%   OPTS = COCO_DEFAULT_UPDATE(OPTS) does nothing.
%

if nargout<1
	error('%s: too few output arguments', mfilename);
end
