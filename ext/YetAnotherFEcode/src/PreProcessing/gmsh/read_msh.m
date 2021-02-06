function [nodes, elements, constr_nodes, nbc_el] = read_msh(mshfile)
[dmn_ids,ebc_ids,nbc_ids]=n_pid(mshfile);
%the different physical domain ids -dmn_ids , (these are numbered from 1-100)- domain ID corresponds to the property ID of these elements 
%the different natural BC ids - ebc_ids,(these are numbered from 101-200)
%the different essential BC ids - nbc_ids (these are numbered from 201-300)


%% generating nodes Matrix
[nodeid,node]=readnodes(mshfile);
nodes = node(nodeid,:);
all_nodes = []; 
for i=1:size(dmn_ids,1)
all_nodes = [all_nodes; readnodeset(mshfile,dmn_ids(i))]; % read node ids of triangular shell elements
end
all_nodes = unique(all_nodes);
red_nodes = setdiff(all_nodes,nodeid');

if~isempty(red_nodes)
    error('repeated nodes in the nodes matrix')
end
%% generating element connectivity Matrix
elements=[];
for i=1:size(dmn_ids,1)

[conn,ele_id,~]=readelements(mshfile,dmn_ids(i)); %return elements belonging to particular physical ID. ele_id, eletyp are row matrices
elements = [elements;[ele_id',conn]];

end

el_id = elements(:,1);  % el_id for mapping
elements = elements(:,2:end); % renumber elements

%% nodes that are constrained
constr_nodes = readnodeset(mshfile,ebc_ids);

%% Forcing elements
% F_nodes = readnodeset(mshfile,nbc_ids(1)); % Nodes containing Force / Pressure
[~,nbc_eid,~]=readelements(mshfile,nbc_ids(1)); % NBC elements
                                    
% pressure on NBC area
nbc_el = [];
for n= nbc_eid   
    i = find(el_id==n-1);% duplicate element for force boundary with index = original index + 1
    nbc_el = [nbc_el; i];    
end
