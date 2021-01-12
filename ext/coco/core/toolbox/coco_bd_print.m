function varargout = coco_bd_print(runid, varargin)

[bd bddat] = coco_bd_read(runid, 'bd', 'bddat');
if nargout>=1
  varargout{1} = bd;
end

if nargin<2
  if isempty(bddat)
    op_names = {};
    op_fmt   = '';
    op_sfmt  = '';
  else
    op_names = bddat.op_names;
    op_fmt   = bddat.op_fmt;
    op_sfmt  = bddat.op_sfmt;
  end
else
  op_names = varargin;
  op_len   = max(cellfun(@numel, op_names), 12);
  op_fmt   = sprintf(' %%%d.4e', op_len);
  op_sfmt  = sprintf(' %%%ds', op_len);
end

PT   = coco_bd_col(bd, 'PT', false);
NRMU = coco_bd_col(bd, '||U||', false);
LAB  = coco_bd_col(bd, 'LAB', false);
TYPE = coco_bd_col(bd, 'TYPE', false);

PARS = {};
for i=1:numel(op_names)
  PAR  = coco_bd_col(bd, op_names{i}, false);
  PARS = [ PARS PAR ]; %#ok<AGROW>
end

fprintf('%5s %12s %6s  %-5s', 'STEP', '||U||', 'LABEL', 'TYPE')
fprintf(op_sfmt, op_names{:});
fprintf('\n');

for i=1:numel(PT)
  if ~isempty(LAB{i})
    fprintf('%5d %12.4e %6d  %-5s', PT{i}, NRMU{i}, LAB{i}, pr_type(TYPE{i}));
    fprintf(op_fmt, PARS{i,:});
    fprintf('\n');
  end
end

end

function tpname = pr_type(tp)
  switch tp
    case {'RO' 'ROS'}
      tpname = '';
    otherwise
      tpname = tp;
  end
end
