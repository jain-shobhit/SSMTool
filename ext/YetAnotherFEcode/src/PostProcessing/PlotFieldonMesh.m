function PlotFieldonMesh(Nodes,Elements,c,varargin)
%--------------------------------------------------------------------------
% Purpose:
%         To plot the profile of a component on mesh
% Synopsis :
%           ProfileonMesh(coordinates,nodes,component)
% Variable Description:
%           coordinates - The nodal coordinates of the mesh
%           -----> coordinates = [X Y Z]
%           nodes - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]
%           component - The components whose profile to be plotted
%           -----> components = a column vector in the order of node
%                               numbers
%
% Original version coded by 
%               Siva Srinivas Kolukula, PhD
% E-mail   :    allwayzitzme@gmail.com
% Last update: 10 April 2018 by Shobhit Jain, 
%                               ETH Zurich 
%                               E: shjain@ethz.ch
%--------------------------------------------------------------------------

dimension = size(Nodes,2) ;  % Dimension of the mesh
nel = length(Elements) ;                  % number of elements
nnode = length(Nodes) ;          % total number of nodes in system
nnel = size(Elements,2);                % number of nodes per element
%
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
Z = zeros(nnel,nel) ;
profile = zeros(nnel,nel) ;

if size(c,2)==1
    component = c;
else
    component = zeros(nnode,1);
end

if nargin == 4
    factor = varargin{1};
    c = factor*c;
end

%
if dimension == 3   % For 3D plots
    elementdim = rank(diff(Nodes(Elements(1,:),1:3)));

    if elementdim ~= 3 % surface in 3D 
        for iel=1:nel
            nd=Elements(iel,:);         % extract connected node for (iel)-th element
            X(:,iel)=Nodes(nd,1);    % extract x value of the node
            Y(:,iel)=Nodes(nd,2);    % extract y value of the node
            Z(:,iel)=Nodes(nd,3) ;   % extract z value of the node
            profile(:,iel) = component(nd') ; % extract component value of the node
        end
        % Plotting the FEM mesh and profile of the given component
        %         figure
        fill3(X,Y,Z,profile)
        rotate3d on ;
        %         title('Profile of component on Mesh') ;
        % Colorbar Setting
        if size(c,2)==1 % if the field is sclar
            SetColorbar
        else
            hold on
            quiver3(Nodes(:,1),Nodes(:,2),Nodes(:,3),c(:,1),c(:,2),c(:,3))
        end
        axis equal
        axis off ;
        
    else % solid in 3D
        fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
        XYZ = cell(1,nel) ;
        profile = XYZ ;
        for e=1:nel
            nd=Elements(e,:);
            X = Nodes(nd,1) ;
            Y = Nodes(nd,2) ;
            Z = Nodes(nd,3) ;
            XYZ{e} = [X  Y Z] ;
            profile{e} = component(nd) ;
        end
        %         figure
        cellfun(@patch,repmat({'Vertices'},1,nel),XYZ,.......
            repmat({'Faces'},1,nel),repmat({fm},1,nel),......
            repmat({'FaceVertexCdata'},1,nel),profile,......
            repmat({'FaceColor'},1,nel),repmat({'interp'},1,nel));
        view(3)
        set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]) ;
        % Colorbar Setting
        if size(c,2)==1 % if the field is sclar
            SetColorbar
        else
            hold on
            quiver3(Nodes(:,1),Nodes(:,2),Nodes(:,3),c(:,1),c(:,2),c(:,3))
        end
        axis equal
        axis off ;
        
        
    end
elseif dimension == 2           % For 2D plots
    elementdim = rank(diff(Nodes(Elements(1,:),1:2)));
    hold on;    
    if elementdim == 2 % surface in 2D
    for iel=1:nel
        nd=Elements(iel,:);         % extract connected node for (iel)-th element
        X(:,iel)=Nodes(nd,1);    % extract x value of the node
        Y(:,iel)=Nodes(nd,2);    % extract y value of the node
        profile(:,iel) = component(nd') ;         % extract component value of the node
    end
    
    % Plotting the FEM mesh and profile of the given component   
        fill(X,Y,profile)
    else
        h = plot(Nodes(:,1),Nodes(:,2),'.-k','Markersize',10,'linewidth',1);
        if size(c,2)==1
            maxcmap = max(component);
            mincmap = min(component);
            cmap = colormap;
            cmapX = [mincmap:(maxcmap - mincmap)/(size(cmap,1)-1):maxcmap]';
            cd = interp1(cmapX,cmap,component);
            interpCmap = cast([255*cd 255*ones(numel(component),1)]','uint8');

            set(h.Edge, 'ColorBinding','interpolated', 'ColorData',interpCmap)
        end
    end
    
    if size(c,2)==1 % if the field is sclar
        colorbar
        %         SetColorbar
    else
        quiver(Nodes(:,1),Nodes(:,2),c(:,1),c(:,2))
    end
    
    axis equal;
    axis off;
    

end









