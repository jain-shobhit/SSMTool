function prob = coco_add_chart_data(prob, cdid, varargin)
%COCO_ADD_CHART_DATA   Add chart data to problem.
%
% COCO_ADD_CHART_DATA(PROB, CDID, VARARGIN)
% Add chart data with identifyer CDID to problem PROB.
%
% VARARGIN = { INIT }
% Add chart data with initial value INIT. Chart data is inherited as a copy
% from the source chart.
%
% VARARGIN = { INIT RESET }
% Add chart data with initial value INIT. Chart data created from a source
% chart is initially set to RESET. Only the chart data of the very first
% chart will be initialized to INIT. All subsequently constructed charts
% will have chart data initialized to RESET.
%
% VARARGIN = { INIT @CP_FUNC FUNC_DATA }
% Add chart data with initial value INIT. Chart data is inherited as the
% return value of the function CP_FUNC; see below.
%
% --------------------------------------------------------------------
% [FUNC_DATA CHART_DATA] = CP_FUNC(PROB, FUNC_DATA, CHART, CHART_DATA)
%
% See also: COCO_GET_CHART_DATA

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_add_chart_data.m 2839 2015-03-05 17:09:01Z fschild $

cp_han = @default_cp_func;
data   = [];

s     = coco_stream(varargin{:});
init1 = s.get;

if ~isempty(s)
  if isa(s.peek, 'function_handle')
    cp_han = s.get;
    data   = s.get;
  else
    data   = s.get;
    cp_han = @reset_cp_func;
  end
end

if ~isfield(prob, 'efunc')
  prob.efunc = efunc_new([]);
end
if ~isfield(prob.efunc, 'identifyers')
  prob.efunc = efunc_new(prob.efunc);
end

% This looks like a bug, but it isn't. See 'help struct' for more info.
if iscell(data)
  cfunc = struct('identifier', cdid, 'F', cp_han, 'data', {data});
else
  cfunc = struct('identifier', cdid, 'F', cp_han, 'data', data);
end

if ~isfield(prob, 'cfunc')
  prob.cfunc.identifiers = { cdid };
  prob.cfunc.funcs       =  cfunc ;
else
  if any(strcmpi(cdid, { prob.cfunc(:).identifiers } ))
    error('%s: copy function with name ''%s'' already defined', ...
      mfilename, cdid);
  end
  prob.cfunc.identifiers = [ prob.cfunc.identifiers { cdid } ];
  prob.cfunc.funcs       = [ prob.cfunc.funcs         cfunc  ];
end

chart              = prob.efunc.chart;
chart.private.data = [ chart.private.data ; { cdid init1 } ];
prob.efunc.chart   = chart;

end

% these functions must be constructed in this way to avoid having prob in
% the stack frame when constructing the function handles, which would lead
% to prob being saved whenever the returned function handle is saved
function [data cdata] = default_cp_func(prob, data, chart, cdata) %#ok<INUSL>
end

function [data cdata] = reset_cp_func(prob, data, chart, cdata) %#ok<INUSD,INUSL>
cdata = data; % reset chart data to RESET
end
