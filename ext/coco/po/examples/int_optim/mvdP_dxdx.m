function dJ = mvdP_dxdx(x, p)

x1 = x(1,:);
p1 = p(1,:);
p4 = p(4,:);

dJ = zeros(3,3,3,numel(1));
dJ(1,1,1,:) = -2*p4.*x1./p1;

end
