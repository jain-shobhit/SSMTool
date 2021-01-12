function [data, J] = complementary_call_DF(opts, data, func, u, l, v)
% OFUNC_CALL_DF     Evaluate Jacobian of complementary zero/monitor function
%
% See also OFUNC_CALL_F, OFUNC_CALL_FDF

switch func.call_mode
  case 1 % [d J]=DFDX(o,d,x)
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
    [data, F, J] = func.F(opts, data, u, l, v); %#ok<ASGLU>
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
