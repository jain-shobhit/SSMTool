function [pars] = coco_idx2par(opts, idx)
%COCO_IDX2PAR   Return parameter names at given indices.
%
%   PARS = COCO_IDX2PAR(OPTS, IDX) returns the names of the parameters
%   stored at the positions in array IDX. OPTS is coco's options structure
%   and IDX is an array with indices. PARS is a cell array containing the
%   names of the parameters. This function is useful for output purposes.
%
%   See also: COCO_ADD_FUNC, COCO_PAR2IDX

pars = {};

if ~isfield(opts, 'efunc') && ~isfield(opts.efunc, 'idx2par')
	return
end

pars = opts.efunc.idx2par(idx);
pars = reshape(pars, size(idx));
