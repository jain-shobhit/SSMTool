% mesh3DHEXparallelepiped
%
% syntax:
% [nodes, elements, nset] = ...
%           mesh3DHEXparallelepiped(eltype,l,w,t,nx,ny,nz)
%
% Description: create the mesh for a parallelepiped using linear or
% quadratic brick elements.
%
% INPUTS    eltype: 'HEX8' or 'HEX20' for linear/quadratic brick
%                   'TET4' or 'TET10' for linear/quadratic tetrahedron
%                   'WED15' for wedges with quadratic shape functions
%           l,w,t: length, width, thickness
%           nx,ny,nz: number of elements for l,w,t respectively
%           
% OUTPUTS   nodes --> [x, y, z]
%           elements --> [node1_ID, ... , nodeN_ID]
%           nset --> struct containing node sets (IDs) of the 6 external
%           faces. They are ordered as: x=0, y=0, z=0, x=l, y=w, z=t.

function [nodes, elements, nset] = ...
    mesh_3Dparallelepiped(eltype,l,w,t,nx,ny,nz)

nel = nx*ny*nz;

lx = l/nx;
ly = w/ny;
lz = t/nz;

switch eltype
    case {'TET4','TET10'}
        % These are the nodes and the element ABAQUS gives when asked to
        % mesh a 1x1x1 cube. 12 elements. Yes, I know. Weak.
        nodesTET = [ ...
            1.,           1.,           1.;
            1.,           0.,           1.;
            1.,           0.,           0.;
            1.,           1.,           0.;
            0.,           0.,           1.;
            0.,           0.,           0.;
            0.,           1.,           0.;
            0.,           1.,           1.;
            0.497189224,  0.497189224,  0.502810717;
            0.748594642,  0.248594612,  0.251405358;
            0.5,           0.,          0.5;
            0.248594612,  0.248594612,  0.751405358;
            0.5,           0.,           1.;
            0.748594642,  0.248594612,  0.751405358;
            1.,           0.,          0.5;
            0.248594612,  0.748594642,  0.751405358;
            0.5,           1.,           1.;
            0.748594642,  0.748594642,  0.751405358;
            1.,          0.5,           1.;
            0.5,          0.5,           1.;
            0.5,           1.,          0.5;
            0.748594642,  0.748594642,  0.251405358;
            1.,           1.,          0.5;
            0.248594612,  0.248594612,  0.251405358;
            0.5,           0.,           0.;
            1.,          0.5,           0.;
            0.5,          0.5,           0.;
            0.5,           1.,           0.;
            0.248594612,  0.748594642,  0.251405358;
            0.,          0.5,           0.;
            0.,           1.,          0.5;
            0.,          0.5,           1.;
            0.,          0.5,          0.5;
            1.,          0.5,          0.5;
            0.,           0.,          0.5];
        elementsTET = [ ...
            9,  5,  3,  2, 12, 11, 10, 14, 13, 15;
            9,  1,  8,  2, 18, 17, 16, 14, 19, 20;
            9,  8,  1,  4, 16, 17, 18, 22, 21, 23;
            9,  3,  6,  4, 10, 25, 24, 22, 26, 27;
            9,  7,  4,  6, 29, 28, 22, 24, 30, 27;
            9,  8,  4,  7, 16, 21, 22, 29, 31, 28;
            9,  8,  7,  5, 16, 31, 29, 12, 32, 33;
            9,  4,  2,  3, 22, 34, 14, 10, 26, 15;
            9,  7,  6,  5, 29, 30, 24, 12, 33, 35;
            9,  5,  6,  3, 12, 35, 24, 10, 11, 25;
            9,  5,  2,  8, 12, 13, 14, 16, 32, 20;
            9,  2,  4,  1, 14, 34, 22, 18, 19, 23];
        NE=10;
        if strcmp(eltype,'TET4')
            elementsTET = elementsTET(:,1:4);
            NE=4;
        end
    case 'HEX8'
        elnodes = [         % "base" element coordinates
             0     0     0
             1     0     0
             1     1     0
             0     1     0
             0     0     1
             1     0     1
             1     1     1
             0     1     1];
         NE = 8;
    case 'HEX20'
        elnodes = [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1;
                    -1 -1 1;  1 -1 1;  1 1 1;  -1 1 1;
                    0 -1 -1;  1 0 -1; 0 1 -1; -1 0 -1;
                     0 -1 1;   1 0 1;  0 1 1;  -1 0 1;
                    -1 -1 0;  1 -1 0;  1 1 0;  -1 1 0]+1;
        elnodes = elnodes/2;
        NE = 20;
    case 'WED15'
        nodesWED = [0.,           0.,           1.
                    0.,           1.,           1.
                    1.,           0.,           1.
                    1.,           1.,           1.
                    0.,           0.,           0.
                    0.,           1.,           0.
                    1.,           0.,           0.
                    1.,           1.,           0.
                    0.,          0.5,           1.
                    0.5,           0.,           1.
                    0.5,          0.5,           1.
                    0.5,          0.5,           0.
                    0.5,           0.,           0.
                    0.,          0.5,           0.
                    1.,           0.,          0.5
                    0.,           1.,          0.5
                    0.,           0.,          0.5
                    0.5,           1.,           1.
                    1.,          0.5,           1.
                    1.,          0.5,           0.
                    0.5,           1.,           0.
                    1.,           1.,          0.5];
        elementsWED = [...
          	2,  3,  1,  6,  7,  5, 11, 10,  9, 12, 13, 14, 16, 15, 17;
            4,  3,  2,  8,  7,  6, 19, 11, 18, 20, 12, 21, 22, 15, 16];
        NE = 15;
    otherwise
        error('Select a valid element type (HEX8 or HEX20)')
end


% create nodes
switch NE
    case {8,20}
        % element coordinates in physical space
        elnodes(:,1) = elnodes(:,1)*lx;
        elnodes(:,2) = elnodes(:,2)*ly;
        elnodes(:,3) = elnodes(:,3)*lz;
        
        nodes = zeros(nel*NE,3);
        nn = 1;
        for ii = 1:nx
            for jj = 1:ny
                for kk = 1:nz
                    
                    elnodes_temp = elnodes;
                    elnodes_temp(:,1) = elnodes_temp(:,1)+lx*(ii-1);
                    elnodes_temp(:,2) = elnodes_temp(:,2)+ly*(jj-1);
                    elnodes_temp(:,3) = elnodes_temp(:,3)+lz*(kk-1);
                    nodes(nn:nn+NE-1,:) = elnodes_temp;
                    
                    nn = nn+NE;
                    
                end
            end
        end
    case {4,10}
        % element coordinates in physical space
        nodesTET(:,1) = nodesTET(:,1)*lx;
        nodesTET(:,2) = nodesTET(:,2)*ly;
        nodesTET(:,3) = nodesTET(:,3)*lz;
        
        nel = nel*12;
        nodes = zeros(nel*NE,3);
        nn = 1;
        for ii = 1:nx
            for jj = 1:ny
                for kk = 1:nz
                    for hh = 1:size(elementsTET,1)
                        elnodes_temp = nodesTET(elementsTET(hh,:),:);
                        elnodes_temp(:,1) = elnodes_temp(:,1)+lx*(ii-1);
                        elnodes_temp(:,2) = elnodes_temp(:,2)+ly*(jj-1);
                        elnodes_temp(:,3) = elnodes_temp(:,3)+lz*(kk-1);
                        nodes(nn:nn+NE-1,:) = elnodes_temp;
                        nn = nn+NE;
                    end
                end
            end
        end
    case 15
        nodesWED(:,1) = nodesWED(:,1)*lx;
        nodesWED(:,2) = nodesWED(:,2)*ly;
        nodesWED(:,3) = nodesWED(:,3)*lz;
        
        nel = nel*2;
        nodes = zeros(nel*NE,3);
        nn = 1;
        for ii = 1:nx
            for jj = 1:ny
                for kk = 1:nz
                    for hh = 1:size(elementsWED,1)
                        elnodes_temp = nodesWED(elementsWED(hh,:),:);
                        elnodes_temp(:,1) = elnodes_temp(:,1)+lx*(ii-1);
                        elnodes_temp(:,2) = elnodes_temp(:,2)+ly*(jj-1);
                        elnodes_temp(:,3) = elnodes_temp(:,3)+lz*(kk-1);
                        nodes(nn:nn+NE-1,:) = elnodes_temp;
                        nn = nn+NE;
                    end
                end
            end
        end
end

tol = 1e12;
nodes = round(nodes*tol)/tol;

% remove duplicate nodes from 'nodes' matrix
[nodes,~,ic] = unique(nodes,'rows','stable');
% ic = ibdex set where repeated rows appear with the indexing of the first
%      appearing row with the same components
elements = ic;
elements = reshape(elements,NE,nel)';

% idb = reshape(1:length(nodes(:)),3,size(nodes,1))';
% conn = zeros(nel,NE*3);
% for ii = 1:nel
%     conn(ii,:) = reshape(idb(elements(ii,:),:)',1,NE*3);
% end
% nodes    = [(1:size(nodes,1))'    nodes];       % [node id, x, y, z]
% elements = [(1:size(elements,1))' elements];    % [elem id, nodes(20)]
% idb      = [(1:size(nodes,1))'    idb];         % [node id, dofs(3)
% conn     = [(1:size(conn,1))'     conn];        % [elem id, dofs(60)]

idnodes = (1:size(nodes,1));

% external surfaces
tol = 1e-15;
face1_X0_nodes = idnodes(nodes(:,1)==0);
face2_Y0_nodes = idnodes(nodes(:,2)==0);
face3_Z0_nodes = idnodes(nodes(:,3)==0);
face4_XL_nodes = idnodes(abs(nodes(:,1)-l)<tol);
face5_YW_nodes = idnodes(abs(nodes(:,2)-w)<tol);
face6_ZT_nodes = idnodes(abs(nodes(:,3)-t)<tol);

nset = {face1_X0_nodes,face2_Y0_nodes,face3_Z0_nodes,face4_XL_nodes,...
    face5_YW_nodes,face6_ZT_nodes};

nn = size(nodes,1);
ndof = nn*3;
fprintf([' Mesh (' eltype '), ' num2str(nel) ' elements, ' num2str(nn) ...
    ' nodes, ' num2str(ndof) ' dofs.\n Dimensions: L=%.2f, W=%.2f,' ...
    ' T=%.2f \n\n'],l,w,t)



