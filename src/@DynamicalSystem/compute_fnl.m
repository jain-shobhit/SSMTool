function fnl = compute_fnl(obj,x,xd)
% COMPUTE_FNL We compute the nonlinear internal force in a second-order
% mechanical system. 

assert(obj.order == 2, ' fnl can only be computed for second-order systems')

fnl = zeros(obj.n,1);
for j = 1:length(obj.fnl)
    if ~isempty(obj.fnl{j})
        sizej = size(obj.fnl{j});
        if sizej(2)==obj.n
            fnl = fnl + expand_tensor(obj.fnl{j},x); % only displacement dependent
        else
            tmp = expand_tensor(obj.fnl{j},[x;xd]); % both disp and velocity dependent
            fnl = fnl+tmp(1:obj.n);
        end
    end
end

end
