function S = expand_tensor(T,p)
% EXPAND_TENSOR This function evaluates the nonlinearity given by a
% tensor/sptensor object T at the vector p. This function makes use of the
% Sandia Tensor Toolbox

if isempty(T)
    S = 0;
    return
else
    if nnz(T) 
        n = ndims(T);
        P = cell(1,n-1);
        P(:) = {p};
        dims = 2:n;        
        S = double( ttv(T,P,dims) );
    else
        s1 = size(T,1);
        S = zeros(s1,1);
    end    

end

end