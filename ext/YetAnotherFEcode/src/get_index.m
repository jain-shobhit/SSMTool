function index = get_index(nodeIDs,nDOFPerNode)
% get the indices of DOFs associated to element j
% when the element j contains m nodes and each node has n DOFs,
% returns index = [ N1_1; N1_2; ..., N1_n;
%                   N2_1; N2_2; ..., N2_n;
%                   .
%                   .
%                   .
%                   Nm_1; Nm_2; ..., Nm_n]
%    size(index) = [n*m,1]

n = nDOFPerNode;
m = length(nodeIDs);
DOFs = repmat((nodeIDs-1)*n,n,1) + repmat((1:n).',1,m);
index = DOFs(:);
end