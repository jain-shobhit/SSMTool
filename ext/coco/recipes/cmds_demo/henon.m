function [data y] = henon(prob, data, u)
%HENON   COCO-compatible encoding of a zero function for the henon demo
%
% The return argument y equals zero when the point (u(5), u(6)) is the
% image under the Henon map of the point (u(3), u(4)), given the parameter
% values (u(1), u(2)).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: henon.m 2839 2015-03-05 17:09:01Z fschild $

  y  = [u(5)-u(4)-1+u(1)*u(3)^2; u(6)-u(2)*u(3)];
end
