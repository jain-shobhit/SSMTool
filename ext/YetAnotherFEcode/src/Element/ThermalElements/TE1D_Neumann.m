classdef TE1D_Neumann < Element
    properties              % Basis properties derived from Element Class
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 1     % number of DOFs per node
        nNodes = 1          % number of nodes per element
        nDim = 0            % number of dimensions in local coordinates
    end

    properties
        dx          % length of the beam element by default zero
        Material    % Object of class Material
        g           % Neumann Constant  
        r           % Robin condition
    end
        
    methods
        
        function self = TE1D_Neumann(Material,g,r)
            self.Material = Material;
            self.g = g;
            self.r = r;
        end       


        %% GET methods: Define how properties are computed
        
        function dx = get.dx(self)
            dx = norm(diff(self.nodes));
        end
        
        %% Global methods
        function M = mass_matrix(~)
            % this function computes the element mass matrix in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            M = 0;
        end

        function [K, F] = tangent_stiffness_and_force(self,T)
            % this function computes the element stiffness matrix and
            % internal force vector in the global coordinates when the
            % nodes     : matrix containing Nodal coordinates
            % T         : element DOF values in global coordinates
            Te = self.extract_element_data(T);
            [K, F] = self.tangent_stiffness_and_force_local(Te);            
        end
               
        
        function xe = extract_element_data(self,x)
            % x is a vector of full DOFs            
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            xe = x(index,:);
        end
        

        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Element local methods %%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        
        function [K,F] = tangent_stiffness_and_force_local(self,T)
            k = self.Material.THERMAL_CONDUCTIVITY;
            dx = self.dx; %#ok<*PROPLC>
            g = self.g; 
            r = self.r; 
            K = r*k;
            F = k*g*ones(1,1);   
        end

    end
end