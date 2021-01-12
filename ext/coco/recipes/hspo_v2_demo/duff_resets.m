function y = duff_resets(x, p, reset)
%DUFF_RESETS   'hspo'-compatible encoding of duffing reset function.

y      = x;
y(3,:) = 0; % phase reset

end
