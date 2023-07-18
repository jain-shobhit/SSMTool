classdef Tet4Element < ContinuumElement
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 3     % number of DOFs per node
        nNodes = 4          % number of nodes per element
        nDim = 3            % number of dimensions in local coordinates
        nelDOFs
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType = 'TET'
    end
    
    properties
        thickness = 1       % element thickness  
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Tet4Element(Material)
            % _____________________________________________________________
            %
            % SELF-FUNCTION
            % self = Tet4Element(Material,Ngauss)
            % defines element's properties
            %______________________________________________________________
            Ngauss = 1; % note: shape function derivatives are constants,
                        % and M is a lumped-parameter mass matrix. No
                        % quadrature integration is actually needed.
            self.thickness = 1;
            self.nelDOFs = self.nNodes * self.nDOFPerNode;
            ContinuumElementConstructor(self, Material, Ngauss);
        end
        
        function Mel = mass_matrix(self)
            % _____________________________________________________________
            %
            % Mel = mass_matrix_global(self,nodes,~);
            % Mel: element-level mass matrix (in global coordinates)
            %______________________________________________________________
            rho = self.Material.DENSITY;
            m = self.vol * rho / 4; % lumped masses are better for TET4
            Mel = sparse(eye(12)*m);
        end
        
        function [G,detJ,dH] = shape_function_derivatives(self, X)
            %______________________________________________________________
            %
            % [G,detJ,dH] = shape_function_derivatives(self)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy uz vx vy vz wx wy wz}' (ux=du/dx...)
            % and p={u1,v1,w1,...,u4,v4,w4}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            xyz = self.nodes;
            % First order thetrahedron. Shape functions:
            %   N = [(1-g-h-r), g, h, r].';
            % Shape function derivatives in natural coordinates:
            dHn = [ -1, 1, 0, 0;
                    -1, 0, 1, 0;
                    -1, 0, 0, 1];    
            
            J = dHn*xyz;
            J1 = [0 0 0; 0 0 0; 0 0 0];
            J1(1,1) = J(2,2)*J(3,3) - J(2,3)*J(3,2);
            J1(2,1) = J(2,3)*J(3,1) - J(2,1)*J(3,3);
            J1(3,1) = J(2,1)*J(3,2) - J(2,2)*J(3,1);
            J1(1,2) = J(1,3)*J(3,2) - J(1,2)*J(3,3);
            J1(2,2) = J(1,1)*J(3,3) - J(1,3)*J(3,1);
            J1(3,2) = J(1,2)*J(3,1) - J(1,1)*J(3,2);
            J1(1,3) = J(1,2)*J(2,3) - J(1,3)*J(2,2);
            J1(2,3) = J(1,3)*J(2,1) - J(1,1)*J(2,3);
            J1(3,3) = J(1,1)*J(2,2) - J(1,2)*J(2,1);
            detJ = J(1,1)*J1(1,1) + J(1,2)*J1(2,1) + J(1,3)*J1(3,1);
            J1 = J1/detJ;
            dH = J1*dHn;   	% derivatives in physical coordinates,
                            % 3x4 matrix, [dNi_dx; dNi_dy; dNi_dz]
                            % with i = 1...4
            G = self.initialization.G;
            G(1:3,1:3:12) = dH;
            G(4:6,2:3:12) = dH;
            G(7:9,3:3:12) = dH;
        end
        
    end
    
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 4-NODED TETRAHEDRON
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.6 Triangular, tetrahedral, and wedge elements
            g = X(1);
            h = X(2);
            r = X(3);
            N = [(2*(1-g-h-r)-1)*(1-g-h-r)
                    (2*g-1)*g
                    (2*h-1)*h
                    (2*r-1)*r];
        end
        
        function X = natural_coordinates
            X = [ ...
                0 0 0
                1 0 0
                0 1 0
                0 0 1];
        end
        
    end

        
end % classdef
    
