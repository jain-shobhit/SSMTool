function DS = expand_tensor_derivative(T,p)
% EXPAND_TENSOR_DERIVATIVE This function evaluates the Jacobian of the
% nonlinearity given by a tensor/sptensor object T at the vector p.
% This function makes use of the Sandia Tensor Toolbox

if isempty(T)
    DS = sparse(1,1);
    return
else
    s1 = size(T,1); % number of rows in the jacobian
    s2 = length(p); % number of colums in the jacobian
    
    if nnz(T)==0
        DS = sparse(s1,s2);
        return
    end
    
    DS = sptensor([s1, s2]);    
    if nnz(T) % only proceed if there are nonzero entries in T
        n = ndims(T);
        % To evaluate the Jacobian, we multiply T by p along the last n-1 dims,
        % taken n-2 at a time and sum the outputs
        P = cell(1,n-2);
        P(:) = {p};
        DIMS = nchoosek(2:n, n-2);
        for j = 1:size(DIMS,1)
            dims = DIMS(j,:);
            DS = DS + ttv(T,P,dims);
        end
    end
    
    % convert tensor into a matrix
    if isa(DS,'sptensor')
        if ~isempty(DS.subs)
            DS = sparse(DS.subs(:,1),DS.subs(:,2),DS.vals, s1, s2);
        else
            DS = sparse(s1,s2);
        end
    else
        DS = DS.data;
    end
    
end

