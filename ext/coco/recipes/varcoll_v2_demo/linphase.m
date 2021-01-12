function [data y] = linphase(prob, data, xp)
%LINPHASE   COCO-compatible encoding of lin phase condition.
  y = data.vphase*(xp(1:3)-xp(4:6));
end
