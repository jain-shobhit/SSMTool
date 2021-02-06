function [subsGlobal, vals] = sparsify(X,globalIND, sumDIMS)
% The idea is that the tensor X represents the coefficients of a
% multivariate polynomial and would be acting on the same vector along the
% dimensions specified by sumDIMS.
%
% This function identifies the entries in tensor X along dimensions in
% vector sumDIMS and sums over them. The sum is stored in the unique
% locations, thus saving memory.
%
% globalIND is a cell array with size [ndims(X),1]. The j-th element of the
% array contains the global DOFs for the j-th dimension of X.

subsLocal = tt_ind2sub(size(X), (1:prod(size(X)))'); %#ok<PSIZE>

if isempty(globalIND)
    subsGlobal = subsLocal;
else
    subsGlobal = replace(subsLocal,globalIND);
end

if ~isempty(sumDIMS)
    % replace summable indices with their sorted counterparts
    subsLocal(:,sumDIMS) = sort(subsLocal(:,sumDIMS),2);
    subsGlobal(:,sumDIMS) = sort(subsGlobal(:,sumDIMS),2);
    ind = tt_sub2ind(size(X), subsLocal);
    
    % sum over repeated indices
    X_summed = accumarray(ind, X.data(:));
    
    % identify unique indices and corresponding subscripts
    ind = unique(ind);
    subsGlobal =  subsGlobal(ind,:);
    vals = X_summed(ind);
else
    vals = X.data(:);    
end
    
end

    function [subs_out] = replace(subs_in, globalSUBS)
        subs_out = zeros(size(subs_in));
        for iSub = 1:length(globalSUBS)
            i = 1:length(globalSUBS{iSub});
            j = globalSUBS{iSub};
            for k = 1:length(i)
                subs_out(subs_in(:,iSub)==i(k),iSub)=j(k);
            end
        end
        
    end