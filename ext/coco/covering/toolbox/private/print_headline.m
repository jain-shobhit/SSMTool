function opts = print_headline(opts, prio)
%COCO_PRINT_HEADLINE  Print headline for point information.
%
%   OPTS = COCO_PRINT_HEADLINE(OPTS) prints a headline for the screen
%   output of the bifurcation diagram during computation. This function
%   calls OPTS.CONT.PRINT_HEADLINE to allow printing of a headline for
%   additional output data.
%
%   See also: coco_print_data
%

if nargin<2
  prio = 1;
end

cont  = opts.cont;
bddat = opts.bddat;

coco_print(opts, prio, '\n');

op_names = bddat.op_names;
op_sfmt  = bddat.op_sfmt;

coco_print(opts, prio, '%5s', 'STEP');
coco_print(opts, 2, ' %11s', 'STEP SIZE');
coco_print(opts, prio, '  %8s %12s %6s  %-5s', ...
  'TIME', '||U||', 'LABEL', 'TYPE');
coco_print(opts, prio, op_sfmt, op_names{:});

opts = coco_emit(opts, 'cont_print', 'init', prio);
coco_print(opts, prio, '\n');
