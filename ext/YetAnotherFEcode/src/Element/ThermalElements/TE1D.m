classdef TE1D < Element
    properties              % Basis properties derived from Element Class
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes    
        nDOFPerNode = 1     % number of DOFs per node
        nNodes = 2          % number of nodes per element
        nDim = 1            % number of dimensions in local coordinates
    end
    
    properties
        dx          % length of the beam element by default zeo        
        Material    % Object of class Material
        uniformBodyForce
    end

    
    methods
        
        function self = TE1D(Material)
            self.Material = Material;
        end
        %% GET methods: Define how properties are computed
        
        function dx = get.dx(self)
            dx = norm(diff(self.nodes));
        end
        
        %% Global methods
        function M = mass_matrix(self)
            % this function computes the element mass matrix in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            M = self.mass_matrix_local();
        end

        function [K, F] = tangent_stiffness_and_force(self,T)
            % this function computes the element stiffness matrix and
            % internal force vector in the global coordinates when the
            % nodes     : matrix containing Nodal coordinates
            % T         : element DOF values in global coordinates
            Te = self.extract_element_data(T);
            [K, F] = self.tangent_stiffness_and_force_local(Te);            
        end
               
        function  f = get.uniformBodyForce(self)
            % this function computes the element stiffness matrix and
            % internal force vector in the global coordinates when the
            % nodes : matrix containing Nodal coordinates            
            f = self.dx/2 *ones(2,1);            
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
            alpha = self.Material.CONVECTION_COEFFICIENT;
            Emiss = self.Material.EMISSIVITY;
            StBol = self.Material.STEFAN_BOLTZMANN; 
            l = self.dx;
            KCond = k/l*[1 -1;-1 1];
            
            % Convection initialization
            
            if alpha ~= 0
                KConv = alpha*l*[1 1/2;1/2 1];
            else
                KConv = zeros(2);
            end
            
            % Radiation initialization
            
            if Emiss ~= 0
                [H,DH] = self.Radiation_local_vector(T);
                 H = Emiss*StBol*0.5*l*H;
                 DH = Emiss*StBol*0.5*l*DH; 
            else
                H = zeros(2,1);
                DH = zeros(2);
            end
            
            Klin = KCond + KConv;
            F = Klin*T + H;  
            K = Klin + DH;
        end
        
        function [K,F] = mass_matrix_local(self)
            c = self.Material.HEAT_CAPACITY;             
            l = self.dx;
            
            K = l*c/3*[1 1/2;1/2 1];
            F = zeros(2,1);             
        end
        
        
        
        function [H,DH] = Radiation_local_vector(self,T)
            T1 = T(1);
            T2 = T(2);
            
            H = [T1^4/3 + (4*T1^3*T2)/15 + (T1^2*T2^2)/5 + (2*T1*T2^3)/15 + T2^4/15;
                T1^4/15 + (2*T1^3*T2)/15 + (T1^2*T2^2)/5 + (4*T1*T2^3)/15 + T2^4/3];
            
            DH = [(4*T1^3)/3 + (4*T1^2*T2)/5 + (2*T1*T2^2)/5 + (2*T2^3)/15, (4*T1^3)/15 + (2*T1^2*T2)/5 + (2*T1*T2^2)/5 + (4*T2^3)/15;
                  (4*T1^3)/15 + (2*T1^2*T2)/5 + (2*T1*T2^2)/5 + (4*T2^3)/15,  (2*T1^3)/15 + (2*T1^2*T2)/5 + (4*T1*T2^2)/5 + (4*T2^3)/3];
        end
        
        
        
    end
    
end
