classdef Mesh < handle
    properties
        nElements               % total number of elements in the mesh (including natural boundary elements)
        nNodes                  % total number of nodes in the mesh
        nodes = []              % Matrix containing nodal coordinates, row-wise for each nodes
        EBC = []                % EssentialBoundary class object
        nDOFPerNode = 0         % number of degrees of freedom associated to each node
        nDOFs = 0
        nDim = 0                % no. of dimensions in the physical domain    
        Elements = []           % A structure array containing pointer to 
                                % each element and an indicator which is true
                                % when element belongs to a natural
                                % boundary  and false when it belongs to
                                % the domain of the Mesh.
    end
    
    methods
        function self = Mesh(varargin)
            narginchk(0,1)            
            if nargin==1
                self.nodes = varargin{1};
            end
        end
        
        function [nNodes] = get.nNodes(self)
            nNodes = size(self.nodes,1);
        end
        
        function nDim = get.nDim(self)
            nDim = size(self.nodes,2);
        end
        
        function nElements = get.nElements(self)
            nElements = length(self.Elements);
        end
        
        function [nDOFs] = get.nDOFs(self)
            nDOFs = self.nNodes * self.nDOFPerNode;
        end
        
        function dofs = get_element_indices(self,j)
            % this function computes the indices corresponding to the j-th
            %  element DOFs in the full vector of unknowns
            dofs = self.Elements(j).iDOFs;
        end
        
        function [outdofFull] = get_DOF_from_location(self,outcoord)
            dist = vecnorm(self.nodes - repmat(outcoord,[self.nNodes,1]),2,2);
            [~,outnode] = min(dist);            
            outdofFull = (outnode-1)*self.nDOFPerNode + (1:self.nDOFPerNode);
        end
        
        function create_elements_table(self, elementsConnectivity,elementConstructor,varargin)
            % creates element_table containing element-wise data
            % Elements(j).Nodes = list of node id for j-th element
            %                       (different number of nodes per element is allowed)
            % Elements(j).Object = Element object for j-th element
            % Elements(j).isBoundary = true if j-th Element is Boundary element, false otherwise
            
            narginchk(3,5)
            %% check element compatibility with the mesh
            test = elementConstructor();
            if self.nDOFPerNode==0 % default value,                
                self.nDOFPerNode = test.nDOFPerNode;
            else    % check nDOFPerNode compatibility between elementConstructor and the Mesh
                if ~isequal(self.nDOFPerNode, test.nDOFPerNode)
                    error('nDOFPerNode error: The supplied element must have the same number of DOFs per node (nDOFPerNode) as the Mesh')
                end
            end
            %% parse inputs
            % check whether these are natural boundary elements or not           
            
            p = inputParser;
            addParameter(p,'isBoundary',false,@(x)validateattributes(x,{'logical'},{'nonempty'}))
            parse(p,varargin{:});
            isBoundary = p.Results.isBoundary;
            
            %% Add the supplied elements to given list
            nElements = size(elementsConnectivity,1); %#ok<*PROPLC>
            if ~isempty(elementsConnectivity)
                % initialize elements as a table
                thisElementsTable = struct('Object', cell(nElements,1), 'isBoundary', isBoundary);
                % node IDs for each element
                for j = 1:nElements
                    nodeIDs = elementsConnectivity(j,:);
                    thisElementsTable(j).Object = elementConstructor();
                    thisElementsTable(j).Object.nodeIDs = nodeIDs;
                    thisElementsTable(j).Object.nodes = self.nodes(nodeIDs,:);
                end
                
                % add new elements to exisiting  (possibly empty) set of elements
                self.Elements = [self.Elements; thisElementsTable];
            end
        end
        
        function set_essential_boundary_condition(self,DirichletNodes,ConstrainedDOF,value)
            if isempty(self.EBC)
                self.EBC = EssentialBoundary(self.nDOFs,self.nDOFPerNode);
            end            
            self.EBC.apply_Dirichlet_BC(DirichletNodes,ConstrainedDOF,value);
        end
        
        function reset_boundary(self)
            self.EBC = [];
            disp('Removed Dirichlet Boundary DOFs')
            
            % find boundary elements
            I = find([self.Elements.isBoundary]);            
            if ~isempty(I)
                disp('Following boundary elements will be deleted')
                disp(I)
                self.Elements(I) = [];
            else
                disp('No boundary elements present')
            end
        end
        
    end
    
end
