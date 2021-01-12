function [data, F, J] = complementary_call_FDF(opts, data, func, u, l, v)
% OFUNC_CALL_FDF      Evaluate complementary zero/monitor function and corresponding Jacobian
%
% See also OFUNC_CALL_F, OFUNC_CALL_DF

if nargout<3
  [data, F] = complementary_call_F(opts, data, func, u, l, v);
  return
end

switch func.call_mode
  case 1 % [d f]=F(o,d,x); [d J]=DFDX(o,d,x)
    [data, F] = func.F(opts, data, u, l, v);
    if isempty(func.DFDX)
      if func.vectorised
        [data, J] = coco_ezDFDX('f(o,d,u,l,v)v', opts, data, func.F, u, l, v);
      else
        [data, J] = coco_ezDFDX('f(o,d,u,l,v)',  opts, data, func.F, u, l, v);
      end
    else
      [data, J] = func.DFDX(opts, data, u, l, v);
      if func.chkdrv
        if func.vectorised
          [data, J2] = coco_ezDFDX('f(o,d,u,l,v)v', opts, data, func.F, u, l, v);
        else
          [data, J2] = coco_ezDFDX('f(o,d,u,l,v)',  opts, data, func.F, u, l, v);
        end
        efunc_chkDeriv(opts, func, J, J2);
      end
    end
  case 2 % [d f J]=FDF(o,d,x)
    [data, F, J] = func.F(opts, data, u, l, v);
    if func.chkdrv
      if func.vectorised
        [data, J2] = coco_ezDFDX('f(o,d,u,l,v)v', opts, data, func.F, u, l, v);
      else
        [data, J2] = coco_ezDFDX('f(o,d,u,l,v)',  opts, data, func.F, u, l, v);
      end
      efunc_chkDeriv(opts, func, J, J2);
    end
end

end
