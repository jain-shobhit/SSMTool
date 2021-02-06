function PlotMesh(Nodes,Elements,show)
%--------------------------------------------------------------------------
% Purpose:
%         To plot 2D and 3D Finite Element Method Mesh,currently
%         this function only works for meshes with a constant number of
%         nodes per element
% Synopsis :
%           PlotMesh(coordinates,nodes)
% Variable Description:
%           coordinates - The nodal coordinates of the mesh
%           -----> coordinates = [X Y Z]
%           nodes - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]
%           show - to dispaly nodal and element numbers
%                  0 (default) - do not display
%                  1           - display
%--------------------------------------------------------------------------

if nargin == 2
    show = 0 ;
end

dimension = size(Nodes,2) ;  % Dimension of the mesh
nel = size(Elements,1) ;                  % number of elements
nnode = length(Nodes) ;          % total number of nodes in system
nnel = size(Elements,2);                % number of nodes per element

if dimension == 3   % For 3D plots
    elementdim = rank(diff(Nodes(Elements(1,:),1:3)));
    
    
    if elementdim == 3 % solid in 3D when we simply plot the skin elements
        faces = getSkin3D(Elements);
        Elements = faces.';
        nel = size(Elements,1);      % total number of faces
        nnel = size(Elements,2);     % number of nodes per face
    end
    
    X = Nodes(Elements',1); X = reshape(X, nnel, nel);
    Y = Nodes(Elements',2); Y = reshape(Y, nnel, nel);
    Z = Nodes(Elements',3); Z = reshape(Z, nnel, nel);
    
    if elementdim ~= 3 % surface or line in 3D        
        [X,Y,Z] = tune_coordinates(X,Y,Z);
    end
    
    patch(X,Y,Z,'w','FaceAlpha',1.0,'EdgeAlpha',1,...
        'EdgeColor','k','LineStyle','-','DisplayName','Mesh');
    view(3)
    set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]) ;
    
    
elseif dimension == 2           % For 2D plots
    elementdim = rank(diff(Nodes(Elements(1,:),1:2))); % dimension of element
    hold on
    if elementdim == 2
        X = Nodes(Elements',1); X = reshape(X, nnel, nel);
        Y = Nodes(Elements',2); Y = reshape(Y, nnel, nel);
        patch(X,Y,'w','DisplayName','Mesh')
    else % line
        plot(Nodes(:,1),Nodes(:,2),'.-k', 'Markersize',10);
    end
    
end
% display Node numbers and Element numbers
if show ~= 0
    k = 1:nnode ;
    nd = k' ;
    for i = 1:nel
        text(X(:,i),Y(:,i),int2str(nd(Elements(i,:))),'fontsize',8,'color','k');
        text(mean(X(:,i)),mean(Y(:,i)),int2str(i),'fontsize',10,'color','r') ;
    end
end
rotate3d on;
axis equal;
axis off;
end

function [X,Y,Z] = tune_coordinates(X,Y,Z)
% ONLY for UNDEFORMED mesh:
% to avoid visualization problems when the surface is contained in
% a plane parallel to the main planes (xy,xz,yz):
if sum(X(:)==X(1))==numel(X)
    X(1)=X(1)+1e-15;
elseif sum(Y(:)==Y(1))==numel(Y)
    Y(1)=Y(1)+1e-15;
elseif sum(Z(:)==Z(1))==numel(Z)
    Z(1)=Z(1)+1e-15;
end
end