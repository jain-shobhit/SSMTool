classdef Tet10Element < ContinuumElement
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 3     % number of DOFs per node
        nNodes = 10         % number of nodes per element
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
        
        function self = Tet10Element(Material, Ngauss)
            % _____________________________________________________________
            %
            % SELF-FUNCTION
            % self = Tet10Element(Material,Ngauss)
            % defines element's properties
            %______________________________________________________________
            if nargin == 1
                Ngauss = 2;
            end
            self.thickness = 1;
            self.nelDOFs = self.nNodes * self.nDOFPerNode;
            ContinuumElementConstructor(self, Material, Ngauss);
        end
        
        function [G,detJ,dH] = shape_function_derivatives(self, X)
            %______________________________________________________________
            %
            % [G,detJ,dH] = shape_function_derivatives(self,g,h,r)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy uz vx vy vz wx wy wz}' (ux=du/dx...)
            % and p={u1,v1,w1,...,u10,v10,w10}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            g = X(1);
            h = X(2);
            r = X(3);
            xyz = self.nodes;
            % shape function derivatives in natural coordinates
            dHn = [0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0];
            dHn(1,1)=4*g+4*h+4*r-3; dHn(1,2)=4*g-1;         dHn(1,5)=4-4*h-4*r-8*g;
            dHn(1,6)=4*h;           dHn(1,7)=-4*h;          dHn(1,8)=-4*r;
            dHn(1,9)=4*r;           dHn(2,1)=4*g+4*h+4*r-3; dHn(2,3)=4*h-1;
            dHn(2,5)=-4*g;          dHn(2,6)=4*g;           dHn(2,7)=4-8*h-4*r-4*g;
            dHn(2,8)=-4*r;          dHn(2,10)=4*r;          dHn(3,1)=4*g+4*h+4*r-3;
            dHn(3,4)=4*r-1;         dHn(3,5)=-4*g;          dHn(3,7)=-4*h;
            dHn(3,8)=4-4*h-8*r-4*g; dHn(3,9)=4*g;           dHn(3,10)=4*h;
            
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
                            % 3x10 matrix, [dNi_dx; dNi_dy; dNi_dz]
                            % with i = 1...10
            G = self.initialization.G;
            G(1:3,1:3:30) = dH;
            G(4:6,2:3:30) = dH;
            G(7:9,3:3:30) = dH;
        end
        
    end % methods
    
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 10-NODED TETRAHEDRON
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.6 Triangular, tetrahedral, and wedge elements
            g = X(1);
            h = X(2);
            r = X(3);
            N = [(2*(1-g-h-r)-1)*(1-g-h-r)
                    (2*g-1)*g
                    (2*h-1)*h
                    (2*r-1)*r
                    4*(1-g-h-r)*g
                    4*g*h
                    4*(1-g-h-r)*h
                    4*(1-g-h-r)*r
                    4*g*r
                    4*h*r];
        end
        
        function X = natural_coordinates
            X = [ ...
                0 0 0
                1 0 0
                0 1 0
                0 0 1
                .5 0 0
                .5 .5 0
                0 .5 0
                0 0 .5
                .5 0 .5
                0 .5 .5];
        end
        
    end

        
end % classdef
    
