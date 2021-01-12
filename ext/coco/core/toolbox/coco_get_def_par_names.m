function [names] = coco_get_def_par_names(name, idx)
%COCO_GET_DEF_PAR_NAMES   Create list of default parameter names.
%
%   NAMES = COCO_GET_DEF_PAR_NAMES(NAME, IDX) creates a cell array of
%   indexed parameter names. The string NAME defines the base name and the
%   integer array idx the set of indices for which to create parameter
%   names. For example, the call
%
%   COCO_GET_DEF_PAR_NAMES('PAR', 1:3)
%
%   will return the 1-by-3 cell array { 'PAR(1)' 'PAR(2)' 'PAR(3)' }.

names = {};

for i=1:numel(idx)
	n = sprintf('%s(%d)', name, idx(i));
	names = { names{:} n };
end
