function [idx] = coco_par2idx(opts, pars, sloppy)
%COCO_PAR2IDX   Return indices of parameters.
%
%   IDX = COCO_PAR2IDX(OPTS, PARS) returns the indices of the parameters as
%   specified by PARS. OPTS is coco's options structure and PARS is either
%   a string of a cell-array of strings. IDX is an array containing the
%   positions of the parameters. This function is used by toolboxes to
%   compute the index that was assigned to internal or external parameters
%   by COCO_ADD_FUNC.
%
%   See also: COCO_ADD_FUNC, COCO_IDX2PAR

if nargin<3
  sloppy = false;
else
  sloppy = strcmp('sloppy', sloppy);
end

idx = [];

if ~isfield(opts, 'efunc') && ~isfield(opts.efunc, 'idx2par')
	return
end

if ischar(pars)
	pars = { pars };
end

for par = pars(:)'
	iidx = find( strcmp(par{1}, opts.efunc.idx2par) );
	if isempty(iidx)
    if sloppy
      iidx = 0; % return illegal index
    else
      error('%s: parameter ''%s'' not found', mfilename, par{1});
    end
	elseif numel(iidx)>1
		error('%s: multiple definition of parameter ''%s''', mfilename, par{1});
	end
	idx = [ idx iidx ]; %#ok<AGROW>
end
idx = reshape(idx, size(pars));
