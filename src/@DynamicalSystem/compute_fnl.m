function fnl = compute_fnl(obj,x,xd)
% COMPUTE_FNL We compute the nonlinear internal force in a second-order
% mechanical system. 

assert(obj.order == 2, ' fnl can only be computed for second-order systems')

fnl = zeros(obj.n,1);
for j = 1:length(obj.fnl)
    if size(obj.fnl(j).ind,2) == obj.N % check if the nonlinearity is velocity dependent as well
        z = [x;xd];
    else
        z = x;
    end
    fnl = fnl + expand_multiindex(obj.fnl(j),z);
end