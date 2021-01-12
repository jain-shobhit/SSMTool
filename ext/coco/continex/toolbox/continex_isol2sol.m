function opts = continex_isol2sol(opts, oid, varargin)
% PROB = CONTINEX_ISOL2SOL(PROB, OID, VARARGIN)
%
% VARARGIN = { @F [FARGS ...] X0 [PNAMES] P0 [T0] }
% FARGS   = 'FPAR', PAR

s = coco_stream(varargin{:});

fhan = s.get;
assert(isa(fhan, 'function_handle'), '%s: invalid input for ''@f''', mfilename), 

% process options
fpars = cell(1,0);
while ischar(s.peek)
  fopt = s.get;
  switch lower(fopt)
    case 'fpar'
      fpars{end+1} = s.get; %#ok<AGROW>
    otherwise
      error('%s: unknown option ''%s''', mfilename, fopt);
  end
end

sol.x0 = s.get;
data.pnames = {};
if iscellstr(s.peek('cell'))
  data.pnames = s.get('cell');
end
sol.p0 = s.get;
if isnumeric(s.peek)
  sol.t0 = s.get;
else
  sol.t0 = [];
end

% initialise toolbox data
data.F     = fhan;
data.fpars = fpars;
data.x_idx = (1:numel(sol.x0))';
data.p_idx = numel(sol.x0)+(1:numel(sol.p0))';

opts = continex_create(opts, oid, data, sol);

end
