% findnode
% n = findnode(xt,yt,zt,nodes2search)
% 
% Returns the id of a node at coordinates (xt,yt,zt).
% INPUTS:   coordinates xt,yt,zt
%           nodes2search: matrix containing nodes' [id x y z]

function n = find_node(xt,yt,zt,nodes2search)

xyz = repmat([xt yt zt],size(nodes2search,1),1);
sn = abs(nodes2search - xyz);
[~,ind]=min(sum(sn,2));
node_labels = 1:size(nodes2search,1);
n = node_labels(ind);

