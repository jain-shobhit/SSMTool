function y = duff(x, p, mode)
%DUFF   'hspo'-compatible encoding of the duffing vector field with bang-bang excitation.

x1 = x(1,:);
x2 = x(2,:);
la = p(1,:);
al = p(2,:);
ep = p(3,:);
A  = p(4,:);
om = p(5,:);

switch mode
  case 'neg'
    A=-A;
end

y(1,:) = x2;
y(2,:) = A-la.*x2-al.*x1-ep.*x1.^3;
y(3,:) = ones(1,numel(om));

end
