function h = PlotFieldonDeformedMesh(Nodes,Elements,disp,varargin)
%--------------------------------------------------------------------------
% Purpose:
%         To plot the profile of a component on deformed mesh, currently
%         this function only works for meshes with a constant number of
%         nodes per element
% Variable Description:
%           Nodes - The nodal coordinates of the mesh
%           -----> Nodes = [X Y Z]
%           Elements - The nodal connectivity of the elements
%           -----> Elements = [node1 node2......]
%           disp -  Nodal displacements [UX UY UZ]
%           Optional parameters:
%           color - the color of the mesh (black by default)
%           factor - Amplification factor (Change accordingly, trial)
%           component -  The components whose profile to be plotted 
%           -----> components  can be given in the following form: 
%                               1. a column vector in the order of node
%                               numbers 
%                               2. empty array [] in which case it is
%                               treated as zeros
%                               3. 'U' : taken as norm of displacements in
%                               'disp' (default)
%                               4. 'U1', 'U2' or 'U3': disp components:
%                               disp(:,1), disp(:,2) or disp (:,3)
%                               repectively.
%--------------------------------------------------------------------------


%%
[meshcolor,factor,c] = parse_inputs(varargin{:});

nnodes = size(Nodes,1);      % number of nodes
dimension = size(Nodes,2) ;  % Dimension of the mesh
elementdim = rank(diff(Nodes(Elements(1,:),:))); % Dimension of elements

nel = size(Elements,1);      % total number of elements   
nnel = size(Elements,2);     % number of nodes per element


if isempty(c)
    c = zeros(nnodes,1);
end

hold on

if dimension == 3   % For 3D plots   
    ux = disp(:,1) ;
    uy = disp(:,2) ;
    uz = disp(:,3) ;
    d = sqrt(ux.^2 + uy.^2 + uz.^2);
    
    if elementdim == 3 % solid in 3D when we simply plot the skin elements    
        faces = getSkin3D(Elements);        
        Elements = faces.';
        nel = size(Elements,1);      % total number of faces
        nnel = size(Elements,2);     % number of nodes per face
    end
        
        
    switch c
        case 'U'
            c = d;
        case 'U1'
            c = ux;
        case 'U2'
            c = uy;
        case 'U3'
            c = uz;
    end
    % Preparing for plot with patch (Elements contains the faces in 2D)
    % X,Y,Z are nnel-by-nel matrices containing x,y,z coordinates for
    % each node in the mesh. E.g., X(i,j) contains the x-coordinate of
    % the i-th node of the j-th element in the mesh, i.e., the node
    % with the index Elements(j,i).

    X = Nodes(Elements',1); X = reshape(X, nnel, nel);
    Y = Nodes(Elements',2); Y = reshape(Y, nnel, nel);
    Z = Nodes(Elements',3); Z = reshape(Z, nnel, nel);
    
    UX = ux(Elements',1); UX = reshape(UX, nnel, nel);
    UY = uy(Elements',1); UY = reshape(UY, nnel, nel);
    UZ = uz(Elements',1); UZ = reshape(UZ, nnel, nel);
    profile = c(Elements',1);
    profile = reshape(profile,[nnel length(profile)/nnel]);
    
    % Plotting the profile of a property on the deformed mesh
    defoX = X+factor*UX ;
    defoY = Y+factor*UY ;
    defoZ = Z+factor*UZ ;
    
    view(3); hold on;
    h = patch(defoX,defoY,defoZ,profile,'EdgeColor',meshcolor,...
        'DisplayName','Deformed Mesh');
    rotate3d on;

elseif dimension == 2           % For 2D plots
    ux = disp(:,1) ;
    uy = disp(:,2) ;
    d = sqrt(ux.^2 + uy.^2);
    switch c
            case 'U'
                c = d;
            case 'U1'
                c = ux;
            case 'U2'
                c = uy;               
    end
    hold on;
    
    if elementdim == 2 % surface in 2D
        
        X = Nodes(Elements',1); X = reshape(X,[nnel nel]);
        Y = Nodes(Elements',2); Y = reshape(Y,[nnel nel]);
        UX = ux(Elements',1); UX = reshape(UX,[nnel nel]);
        UY = uy(Elements',1); UY = reshape(UY,[nnel nel]);
        profile = c(Elements',1); 
        profile = reshape(profile,[nnel length(profile)/nnel]);
        % Plotting the profile of a property on the deformed mesh
        defoX = X + factor*UX ;
        defoY = Y + factor*UY ;
        
        h{1} = patch(defoX,defoY,profile,'EdgeColor',meshcolor);
        h{2} = plot(defoX,defoY,'.','Color', meshcolor, 'Markersize',10);
    else
        h = plot(Nodes(:,1)+factor*ux,Nodes(:,2)+factor*uy,'.-','Color', meshcolor, 'Markersize',10);
        c = [];
    end
    
end
        axis equal;
        axis off;
    % Colorbar Setting
    if ~isempty(c)
        SetColorbar
    end

end

function [meshcolor,factor,c] = parse_inputs(varargin)

defaultColor = 'k';
defaultFactor = 1;
defaultComponent = 'U'; % plot norm of displacement
p = inputParser;

addParameter(p,'color',defaultColor)
addParameter(p,'factor',defaultFactor,@(x)validateattributes(x, ...
    {'numeric'},{'nonempty','positive'}) );
addParameter(p,'component',defaultComponent);
parse(p,varargin{:});

meshcolor = p.Results.color;
factor = p.Results.factor;
c = p.Results.component;

end
