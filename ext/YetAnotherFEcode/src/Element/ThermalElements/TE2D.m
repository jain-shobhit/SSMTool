classdef TE2D < Element
    properties              % Basis properties derived from Element Class
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes    
        nDOFPerNode = 1     % number of DOFs per node
        nNodes = 3          % number of nodes per element
        nDim = 2            % number of dimensions in local coordinates
    end
    
    properties
        area           % Element's area
        jacobian       % Jacobian
        Material    % Object of class Material
        uniformBodyForce
    end
    
    
    methods
        
        function self = TE2D(Material)
            self.Material = Material;
        end
        
        function J = get.jacobian(self)
            Re = self.compute_rotation_matrix();
            vR = self.nodes*Re';
            x31 = vR(3,1)-vR(1,1);
            x21 = vR(2,1)-vR(1,1);
            y21 = vR(2,2)-vR(1,2);
            y31 = vR(3,2)-vR(1,2);
            J = [x21 y21; x31 y31];            
        end
        
        function area = get.area(self)
            area = 0.5*det(self.jacobian);
        end
        
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
               
        
        function xe = extract_element_data(self,x)
            % x is a vector of full DOFs            
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            xe = x(index,:);
        end
        
        function R = compute_rotation_matrix(self)
            R = rotation_matrix(self.nodes);
        end
        
        function  f = get.uniformBodyForce(self)
            A = self.area;
            f = A/3 *ones(3,1); % Heat source vector
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Element local methods %%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        
        function [K,F] = tangent_stiffness_and_force_local(self,T)
            k = self.Material.THERMAL_CONDUCTIVITY;  
            alpha = self.Material.CONVECTION_COEFFICIENT;
            Emiss = self.Material.EMISSIVITY;
            StBol = self.Material.STEFAN_BOLTZMANN; 
            A = self.area; 
            
            % Initialization of the conduction matrix 
            dN1 = [-1;-1];
            dN2 = [1;0];
            dN3 = [0;1];
            J = self.jacobian; %#ok<*PROPLC>
            j1 = J\dN1;
            j2 = J\dN2;
            j3 = J\dN3;
            K_11 = k*A*dot(j1,j1);
            K_12 = k*A*dot(j1,j2);
            K_13 = k*A*dot(j1,j3);
            K_22 = k*A*dot(j2,j2);
            K_23 = k*A*dot(j2,j3); 
            K_33 = k*A*dot(j3,j3);
            KCond = [K_11 K_12 K_13; K_12 K_22 K_23; K_13 K_23 K_33]; % Stiffness matrix for conduction
            
            %Initialization of the convection matrix
            if alpha ~= 0
               KConv = alpha*A*[1/6 1/12 1/12;1/12 1/6 1/12;1/12 1/12 1/6]; %Stiffness matrix for convection
            else
               KConv = zeros(3); 
            end
            
            %inintialization of radiation term
            if Emiss ~= 0
                [H, DH] = self.Radiation_local_vector(T);
                H = (Emiss*StBol*2*A)*H;
                DH = (Emiss*StBol*2*A)*DH;
            else
                H = zeros(3,1);
                DH = zeros(3,3);
            end
            
            % Local Stiffness Matrix
            Klin = KCond + KConv;
            
            F =  Klin*T + H;
            K = Klin + DH;
        end
        
        function [K] = mass_matrix_local(self)
            c = self.Material.HEAT_CAPACITY;             
            A = self.area;            
            K = c*A*[1/6 1/12 1/12;1/12 1/6 1/12;1/12 1/12 1/6]; %#ok<*PROP>
        end
         
        function [H,DH] = Radiation_local_vector(~,T)
            T1 = T(1);
            T2 = T(2);
            T3 = T(3);
            H = [T1^4/42 + (2*T1^3*T2)/105 + (2*T1^3*T3)/105 + (T1^2*T2^2)/70 + (T1^2*T2*T3)/70 + (T1^2*T3^2)/70 + (T1*T2^3)/105 + (T1*T2^2*T3)/105 + (T1*T2*T3^2)/105 + (T1*T3^3)/105 + T2^4/210 + (T2^3*T3)/210 + (T2^2*T3^2)/210 + (T2*T3^3)/210 + T3^4/210;
                 T1^4/210 + (T1^3*T2)/105 + (T1^3*T3)/210 + (T1^2*T2^2)/70 + (T1^2*T2*T3)/105 + (T1^2*T3^2)/210 + (2*T1*T2^3)/105 + (T1*T2^2*T3)/70 + (T1*T2*T3^2)/105 + (T1*T3^3)/210 + T2^4/42 + (2*T2^3*T3)/105 + (T2^2*T3^2)/70 + (T2*T3^3)/105 + T3^4/210;
                 T1^4/210 + (T1^3*T2)/210 + (T1^3*T3)/105 + (T1^2*T2^2)/210 + (T1^2*T2*T3)/105 + (T1^2*T3^2)/70 + (T1*T2^3)/210 + (T1*T2^2*T3)/105 + (T1*T2*T3^2)/70 + (2*T1*T3^3)/105 + T2^4/210 + (T2^3*T3)/105 + (T2^2*T3^2)/70 + (2*T2*T3^3)/105 + T3^4/42];
            
             DH = [(2*T1^3)/21 + (2*T1^2*T2)/35 + (2*T1^2*T3)/35 + (T1*T2^2)/35 + (T1*T2*T3)/35 + (T1*T3^2)/35 + T2^3/105 + (T2^2*T3)/105 + (T2*T3^2)/105 + T3^3/105  (2*T1^3)/105 + (T1^2*T2)/35 + (T1^2*T3)/70 + (T1*T2^2)/35 + (2*T1*T2*T3)/105 + (T1*T3^2)/105 + (2*T2^3)/105 + (T2^2*T3)/70 + (T2*T3^2)/105 + T3^3/210  (2*T1^3)/105 + (T1^2*T2)/70 + (T1^2*T3)/35 + (T1*T2^2)/105 + (2*T1*T2*T3)/105 + (T1*T3^2)/35 + T2^3/210 + (T2^2*T3)/105 + (T2*T3^2)/70 + (2*T3^3)/105;
                  (2*T1^3)/105 + (T1^2*T2)/35 + (T1^2*T3)/70 + (T1*T2^2)/35 + (2*T1*T2*T3)/105 + (T1*T3^2)/105 + (2*T2^3)/105 + (T2^2*T3)/70 + (T2*T3^2)/105 + T3^3/210     T1^3/105 + (T1^2*T2)/35 + (T1^2*T3)/105 + (2*T1*T2^2)/35 + (T1*T2*T3)/35 + (T1*T3^2)/105 + (2*T2^3)/21 + (2*T2^2*T3)/35 + (T2*T3^2)/35 + T3^3/105  T1^3/210 + (T1^2*T2)/105 + (T1^2*T3)/105 + (T1*T2^2)/70 + (2*T1*T2*T3)/105 + (T1*T3^2)/70 + (2*T2^3)/105 + (T2^2*T3)/35 + (T2*T3^2)/35 + (2*T3^3)/105;
                  (2*T1^3)/105 + (T1^2*T2)/70 + (T1^2*T3)/35 + (T1*T2^2)/105 + (2*T1*T2*T3)/105 + (T1*T3^2)/35 + T2^3/210 + (T2^2*T3)/105 + (T2*T3^2)/70 + (2*T3^3)/105  T1^3/210 + (T1^2*T2)/105 + (T1^2*T3)/105 + (T1*T2^2)/70 + (2*T1*T2*T3)/105 + (T1*T3^2)/70 + (2*T2^3)/105 + (T2^2*T3)/35 + (T2*T3^2)/35 + (2*T3^3)/105     T1^3/105 + (T1^2*T2)/105 + (T1^2*T3)/35 + (T1*T2^2)/105 + (T1*T2*T3)/35 + (2*T1*T3^2)/35 + T2^3/105 + (T2^2*T3)/35 + (2*T2*T3^2)/35 + (2*T3^3)/21]; 
        end

        
    end
    
end
