classdef TE2D_Neumann < Element
    
    properties              % Basis properties derived from Element Class
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes    
        nDOFPerNode = 1     % number of DOFs per node
        nNodes = 2          % number of nodes per element
        nDim = 1            % number of dimensions in local coordinates
    end
    
    properties
        h           % half Element's length -> |J| 
        Material    % Object of class Material
        g           % Neumann Constant  
        r           % Robin condition
    end
    
    methods
        
        function self = TE2D_Neumann(Material,g,r)
            self.Material = Material;
            self.g = g;
            self.r = r;
        end       

        function h = get.h(self)
            Re = self.compute_rotation_matrix();
            vR = self.nodes*Re';
            h = 0.5*sqrt((vR(1,1)-vR(2,1))^2 + (vR(1,2)-vR(2,2))^2);
        end
            
        function M = mass_matrix(~)
            % this function computes the element mass matrix in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            M = zeros(2,2);           
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
        
        function R = compute_rotation_matrix(self)
            R = rotation_matrix(self.nodes);
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Element local methods %%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        
        function [K,F] = tangent_stiffness_and_force_local(self,T)
            k = self.Material.THERMAL_CONDUCTIVITY;
            h = self.h;
            g = self.g; 
            r = self.r;  
            a = self.Material.CONVECTION_COEFFICIENT; 
            K = r*a*h*[2/3 1/3;1/3 2/3];;
            F = K*T + k*g*h*ones(1,1);   
        end

    end
end