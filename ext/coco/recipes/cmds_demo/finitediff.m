function [data y] = finitediff(prob, data, u)
%FINITEDIFF   COCO-compatible encoding of a zero function.
%
% The return argument y equals zero when the elements of u indexed by
% data.dep_idx satisfy a discretization of two coupled ordinary
% differential equations parameterized by the element of u indexed by
% data.par_idx.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: finitediff.m 2839 2015-03-05 17:09:01Z fschild $

dep = u(data.dep_idx);
par = u(data.par_idx);

ff = dep(data.f_idx);
gg = dep(data.g_idx);

f = [par(1)-(par(2)+1)*ff+ff.^2.*gg; par(2)*ff-ff.^2.*gg];

y = data.A*dep+data.B*f;

end
