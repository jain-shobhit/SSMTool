function [data y] = lingap(prob, data, xp)
%LINGAP   COCO-compatible encoding of lin gap condition.
  y = data.gapvec*(xp(1:3)-xp(4:6));
end
