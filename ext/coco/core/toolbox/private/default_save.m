function [opts data] = default_save(opts, varargin)
%COCO_DEFAULT_SAVE  Default implementation listing additional save data.
%
%   [OPTS CLASS_LIST] = COCO_DEFAULT_SAVE(OPTS) returns an empty
%   CLASS_LIST.
%

if nargout<2
	error('%s: too few output arguments', mfilename);
end

data = [];
