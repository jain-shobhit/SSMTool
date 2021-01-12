function fnl = compute_fnl(obj,x,xd)
% COMPUTE_FNL We compute the nonlinear internal force in a second-order
% mechanical system. Currently, we do not treat velocity dependent
% nonlinearity.

assert(obj.order == 2, ' fnl can only be computed for second-order systems')

fnl = zeros(obj.n,1);
for j = 1:length(obj.fnl)
    fnl = fnl + expand_tensor(obj.fnl{j},x);
end
