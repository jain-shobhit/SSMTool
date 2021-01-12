function [data y] = tanh_F(prob, data, u)
%TANH_F   COCO-compatible encoding of zero problem.

x = u(data.x_idx); % Extract problem variables
p = u(data.p_idx); % Extract problem parameters

y = x-tanh(p*data.t)/tanh(p);

end
