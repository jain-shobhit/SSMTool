function [data y xidx] = tanh_ramp(opts, data, xp)
% Ramp function for implicit step-size reduction.

if isempty(data)
  [fdata xidx] = coco_get_func_data(opts, 'continex');
  xidx         = xidx(fdata.x_idx);
  xp           = xp(xidx);
  data.init    = true;
else
  xidx = [];
end

x = sqrt(sum(xp.*xp, 1))-0.9;
y(1,:) = 1*(tanh(5*x)+tanh(3*x)+tanh(x));

