function multi_index = sub2multiind(subs,r)
% this function returns multi-index corresponding to the input subscripts  
% (subs) of a multidimensional array that acts on a vector of r variables.
% e.g. 
% Monomial x_1^2*x_3
% subs: [1 1 3], [1 3 1] or [3 1 1]
% multi_index: [2,0,1]
% for r = 3, we have variables [x_1,x_2,x_3]. The monomial x_1^2*x_3 
% may be represented in terms of multi-dimensional subscripts as [1 1 3],
% [1 3 1] or [3 1 1], but it has a unique representation in the multi-index 
% notation in terms of exponents of each variable as [2 0 1].   

[nSubs,d] = size(subs);
ne = numel(subs);
% 
M = sparse(1:ne,reshape(subs.',1,[]),ones(1,ne),ne,r);
% summation operator
dSum = cell(1,nSubs);
dSum(:) = {sparse(ones(1,d))};
dSum = blkdiag(dSum{:});

% perform rowwise sum taken d rows at a time
multi_index = dSum * M;
end