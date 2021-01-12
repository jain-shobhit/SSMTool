function opts = print_data(opts, chart)
%PRINT_DATA  Save data and print point information to screen.
%
%   OPTS = PRINT_DATA(OPTS) is called whenever a new solution has been
%   computed. If a solution label has been assigned to the current solution
%   point it prints a line with information about the current solution
%   point on screen. This function calls OPTS.CONT.PRINT_DATA to print
%   additional data.
%
%   See also: 
%

%% check if solution info should be printed

if isempty(chart.pt_type)
	prio = 2;
else
  prio = 1;
end

bddat = opts.bddat;

%% print solution info on screen

%  compute computation times in hours, minutes and seconds
ctm  = etime(clock, opts.cont.tm);
ctmh = floor(ctm/3600);
ctm  = ctm - 3600*ctmh;
ctmm = floor(ctm/60);
ctm  = ctm - 60*ctmm;
ctms = floor(ctm);

op_idx   = bddat.op_idx;
op_fmt   = bddat.op_fmt;
op_sfmt  = bddat.op_sfmt;

coco_print(opts, prio, '%5d', chart.pt);
coco_print(opts, 2, ' %11.2e', chart.R);
coco_print(opts, prio, '  %02d:%02d:%02d %12.4e %6s  %-5s', ...
  ctmh, ctmm, ctms, norm(chart.x), lab_name(chart), ...
  spt_name(chart.pt_type));
if isfield(chart, 'p')
  coco_print(opts, prio, op_fmt, chart.p(op_idx(:)));
else
  msg = {'---'};
  coco_print(opts, prio, op_sfmt, msg{ones(1,numel(op_idx))});
end

%  print user defined information
opts = coco_emit(opts, 'cont_print', 'data', prio, chart, chart.x);
coco_print(opts, prio, '\n');

end

function [y] = spt_name(x)

switch(x)
	case {'RO', 'ROS'}, y = ' ';
	otherwise,          y =  x ;
end

end

function [y] = lab_name(chart)

if isfield(chart, 'lab')
  switch(chart.pt_type)
    case 'ROS', y = sprintf('*%d', chart.lab);
    otherwise,  y = sprintf('%d' , chart.lab);
  end
else
  y = '';
end

end
