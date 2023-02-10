function [data, y] = amplitude(prob, data, u) %#ok<INUSL>

xbp = reshape(u, data.zdim, []);
y = xbp(data.dof,:);
y = max(abs(y),[],2);

end

