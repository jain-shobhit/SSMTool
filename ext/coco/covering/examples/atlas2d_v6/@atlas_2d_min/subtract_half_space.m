function chart = subtract_half_space(atlas, chart, test, phi, flag, NB)
%SUBTRACT_HALF_SPACE   Subtract half-space from polygon.
%
% Henderson's algorithm chops a polygon by a line obtained from the
% intersection of a circle around the chart base point and a second circle
% around the projection of the neighboring chart base point onto the
% tangent space of the first chart.
%
% Identical to atlas2d_v5.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: subtract_half_space.m 2839 2015-03-05 17:09:01Z fschild $

k        = find(flag & ~circshift(flag, -1), 1); % Find first edge with intersection
flag     = circshift(flag, -k(1));
test     = circshift(test, -k(1));
chart.s  = circshift(chart.s, [0, -k(1)]);
chart.v  = circshift(chart.v, -k(1));
chart.nb = circshift(chart.nb, [0, -k(1)]);
j        = find(~flag & circshift(flag, -1), 1); % Find second edge with intersection
vx1      = chart.v(j)*chart.s(:,j);
vx2      = chart.v(j+1)*chart.s(:,j+1);
nvx1     = vx1-test(j)/((vx2-vx1)'*phi)*(vx2-vx1);
vx1      = chart.v(end)*chart.s(:,end);
vx2      = chart.v(1)*chart.s(:,1);
nvx2     = vx1-test(end)/((vx2-vx1)'*phi)*(vx2-vx1);
chart.s  = [chart.s(:,1:j), nvx1/norm(nvx1), nvx2/norm(nvx2)];  % Add new directions
chart.v  = [chart.v(1:j); norm(nvx1); norm(nvx2)];              % Add new vertex distances
chart.nb = [chart.nb(1:j+1), NB];                               % Add new edge
ep_flag  = (chart.ep_flag && (chart.pt>0));                     % Chart on computational domain boundary 
chart.bv = find(~ep_flag & (chart.v>atlas.cont.Rmarg*chart.R)); % Update index array of available vertices

end
