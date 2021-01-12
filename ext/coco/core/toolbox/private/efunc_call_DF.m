function [opts, data, chart, J] = efunc_call_DF(opts, data, chart, func, x)
% EFUNC_CALL_DF
%
% Evaluate Jacobian of zero or monitor function using the appropriate
% calling syntax (see COCO_ADD_FUNC).
%
% See also EFUNC_CALL_F, EFUNC_CALL_FDF

switch func.call_mode
  case 1 % [d J]=DFDX(o,d,x)
    if isempty(func.DFDX)
      if func.vectorised
        [data, J] = coco_ezDFDX('f(o,d,x)v', opts, data, func.F, x);
      else
        [data, J] = coco_ezDFDX('f(o,d,x)',  opts, data, func.F, x);
      end
    else
      [data, J] = func.DFDX(opts, data, x);
      if func.chkdrv
        if func.vectorised
          [data, J2] = coco_ezDFDX('f(o,d,x)v', opts, data, func.F, x);
        else
          [data, J2] = coco_ezDFDX('f(o,d,x)',  opts, data, func.F, x);
        end
        efunc_chkDeriv(opts, func, J, J2);
      end
    end
  case 2 % [d f J]=FDF(o,d,x)
    [data, F, J] = func.F(opts, data, x); %#ok<ASGLU>
    if func.chkdrv
      if func.vectorised
        [data, J2] = coco_ezDFDX('f(o,d,x)v', opts, data, func.F, x);
      else
        [data, J2] = coco_ezDFDX('f(o,d,x)',  opts, data, func.F, x);
      end
      efunc_chkDeriv(opts, func, J, J2);
    end
  case 3 % [d c J]=DFDX(o,d,c,x)
    if isempty(func.DFDX)
      if func.vectorised
        [data, chart2, J] = coco_ezDFDX('f(o,d,c,x)v', opts, data, chart, func.F, x);
      else
        [data, chart2, J] = coco_ezDFDX('f(o,d,c,x)',  opts, data, chart, func.F, x);
      end
    else
      [data, chart2, J] = func.DFDX(opts, data, chart, x);
      if func.chkdrv
        if func.vectorised
          [data, chart3, J2] = coco_ezDFDX('f(o,d,c,x)v', opts, data, chart, func.F, x); %#ok<ASGLU>
        else
          [data, chart3, J2] = coco_ezDFDX('f(o,d,c,x)',  opts, data, chart, func.F, x); %#ok<ASGLU>
        end
        efunc_chkDeriv(opts, func, J, J2);
      end
    end
    chart.private.data = chart2.private.data;
  case 4 % [d c f J]=FDF(o,d,c,x)
    [data, chart2, F, J] = func.F(opts, data, chart, x); %#ok<ASGLU>
    chart.private.data = chart2.private.data;
    if func.chkdrv
      if func.vectorised
        [data, chart3, J2] = coco_ezDFDX('f(o,d,c,x)v', opts, data, chart, func.F, x); %#ok<ASGLU>
      else
        [data, chart3, J2] = coco_ezDFDX('f(o,d,c,x)',  opts, data, chart, func.F, x); %#ok<ASGLU>
      end
      efunc_chkDeriv(opts, func, J, J2);
    end
    
end

end
