function [opts, data, chart, F] = efunc_call_F(opts, data, chart, func, x, t)
% EFUNC_CALL_F
%
% Evaluate zero or monitor function using the appropriate calling syntax
% (see COCO_ADD_FUNC).
%
% See also EFUNC_CALL_DF, EFUNC_CALL_FDF

switch func.call_mode
  case {1 2} % [d f]=F(o,d,x)
    [data, F] = func.F(opts, data, x);
  case {3 4} % [d c f]=F(o,d,c,x)
    [data, chart2, F] = func.F(opts, data, chart, x);
    chart.private.data = chart2.private.data;
  case {5 6} % [d f]=F(o,d,x,t)
    [data, F] = func.F(opts, data, x, t);
  case {7 8} % [d c f]=F(o,d,c,x,t)
    [data, chart2, F] = func.F(opts, data, chart, x, t);
    chart.private.data = chart2.private.data;
  case {9 10} % [o d f]=F(o,d,x)
    [opts, data, F] = func.F(opts, data, x);
  case {11 12} % [o d c f]=F(o,d,c,x)
    [opts, data, chart2, F] = func.F(opts, data, chart, x);
    chart.private.data = chart2.private.data;
  case {13 14} % [o d f]=F(o,d,x,t)
    [opts, data, F] = func.F(opts, data, x, t);
  case {15 16} % [o d c f]=F(o,d,c,x,t)
    [opts, data, chart2 F] = func.F(opts, data, chart, x, t);
    chart.private.data = chart2.private.data;
end

end
