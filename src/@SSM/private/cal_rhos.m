function [data, y] = cal_rhos(prob, data, u)
% CAL_RHOS: This function converts cartesian coordinates to polar
% coordinats and return radius
%
% [DATA, Y] = CAL_RHOS(PROB, DATA, U)
%

re = u(1:2:end);
im = u(2:2:end);
y  = sqrt(re.^2+im.^2);
y  = data.scale(:).*y;
end