classdef Hex20Element < ContinuumElement
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 3     % number of DOFs per node
        nNodes = 20         % number of nodes per element
        nDim = 3            % number of dimensions in local coordinates
        nelDOFs
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType = 'HEX'
    end
    
    properties
        thickness = 1       % element thickness  
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Hex20Element(Material, Ngauss)
            % _____________________________________________________________
            %
            % SELF-FUNCTION
            % self = Hex20Element(Material,Ngauss)
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
            % [G,detJ,dH] = G_HEX20(self,g,h,r)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy uz vx vy vz wx wy wz}' (ux=du/dx...)
            % and p={u1,v1,w1,...,u20,v20,w20}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            g = X(1);
            h = X(2);
            r = X(3);
            xyz = self.nodes;
            % shape functions derivatives (ABAQUS ORDER)
            dHn = [ ...
                (g/8 - 1/8)*(h - 1)*(r - 1) + ((h - 1)*(r - 1)*(g + h + r + 2))/8,   (g/8 + 1/8)*(h - 1)*(r - 1) - ((h - 1)*(r - 1)*(h - g + r + 2))/8, - (g/8 + 1/8)*(h + 1)*(r - 1) - ((h + 1)*(r - 1)*(g + h - r - 2))/8, - (g/8 - 1/8)*(h + 1)*(r - 1) - ((h + 1)*(r - 1)*(g - h + r + 2))/8, - (g/8 - 1/8)*(h - 1)*(r + 1) - ((h - 1)*(r + 1)*(g + h - r + 2))/8, - (g/8 + 1/8)*(h - 1)*(r + 1) - ((h - 1)*(r + 1)*(g - h + r - 2))/8, (g/8 + 1/8)*(h + 1)*(r + 1) + ((h + 1)*(r + 1)*(g + h + r - 2))/8, (g/8 - 1/8)*(h + 1)*(r + 1) + ((h + 1)*(r + 1)*(g - h - r + 2))/8, - (g/4 - 1/4)*(h - 1)*(r - 1) - ((g + 1)*(h - 1)*(r - 1))/4,                               (h/4 - 1/4)*(h + 1)*(r - 1), (g/4 - 1/4)*(h + 1)*(r - 1) + ((g + 1)*(h + 1)*(r - 1))/4,                                -(h/4 - 1/4)*(h + 1)*(r - 1), (g/4 - 1/4)*(h - 1)*(r + 1) + ((g + 1)*(h - 1)*(r + 1))/4,                                -(h/4 - 1/4)*(h + 1)*(r + 1), - (g/4 - 1/4)*(h + 1)*(r + 1) - ((g + 1)*(h + 1)*(r + 1))/4,                               (h/4 - 1/4)*(h + 1)*(r + 1),                                -(r/4 - 1/4)*(h - 1)*(r + 1),                               (r/4 - 1/4)*(h - 1)*(r + 1),                                -(r/4 - 1/4)*(h + 1)*(r + 1),                               (r/4 - 1/4)*(h + 1)*(r + 1);
                (g/8 - 1/8)*(h - 1)*(r - 1) + (g/8 - 1/8)*(r - 1)*(g + h + r + 2), - (g/8 + 1/8)*(h - 1)*(r - 1) - (g/8 + 1/8)*(r - 1)*(h - g + r + 2), - (g/8 + 1/8)*(h + 1)*(r - 1) - (g/8 + 1/8)*(r - 1)*(g + h - r - 2),   (g/8 - 1/8)*(h + 1)*(r - 1) - (g/8 - 1/8)*(r - 1)*(g - h + r + 2), - (g/8 - 1/8)*(h - 1)*(r + 1) - (g/8 - 1/8)*(r + 1)*(g + h - r + 2),   (g/8 + 1/8)*(h - 1)*(r + 1) - (g/8 + 1/8)*(r + 1)*(g - h + r - 2), (g/8 + 1/8)*(h + 1)*(r + 1) + (g/8 + 1/8)*(r + 1)*(g + h + r - 2), (g/8 - 1/8)*(r + 1)*(g - h - r + 2) - (g/8 - 1/8)*(h + 1)*(r + 1),                                -(g/4 - 1/4)*(g + 1)*(r - 1), (h/4 - 1/4)*(g + 1)*(r - 1) + ((g + 1)*(h + 1)*(r - 1))/4,                               (g/4 - 1/4)*(g + 1)*(r - 1), - (h/4 - 1/4)*(g - 1)*(r - 1) - ((g - 1)*(h + 1)*(r - 1))/4,                               (g/4 - 1/4)*(g + 1)*(r + 1), - (h/4 - 1/4)*(g + 1)*(r + 1) - ((g + 1)*(h + 1)*(r + 1))/4,                                -(g/4 - 1/4)*(g + 1)*(r + 1), (h/4 - 1/4)*(g - 1)*(r + 1) + ((g - 1)*(h + 1)*(r + 1))/4,                                -(r/4 - 1/4)*(g - 1)*(r + 1),                               (r/4 - 1/4)*(g + 1)*(r + 1),                                -(r/4 - 1/4)*(g + 1)*(r + 1),                               (r/4 - 1/4)*(g - 1)*(r + 1);
                (g/8 - 1/8)*(h - 1)*(r - 1) + (g/8 - 1/8)*(h - 1)*(g + h + r + 2), - (g/8 + 1/8)*(h - 1)*(r - 1) - (g/8 + 1/8)*(h - 1)*(h - g + r + 2),   (g/8 + 1/8)*(h + 1)*(r - 1) - (g/8 + 1/8)*(h + 1)*(g + h - r - 2), - (g/8 - 1/8)*(h + 1)*(r - 1) - (g/8 - 1/8)*(h + 1)*(g - h + r + 2),   (g/8 - 1/8)*(h - 1)*(r + 1) - (g/8 - 1/8)*(h - 1)*(g + h - r + 2), - (g/8 + 1/8)*(h - 1)*(r + 1) - (g/8 + 1/8)*(h - 1)*(g - h + r - 2), (g/8 + 1/8)*(h + 1)*(r + 1) + (g/8 + 1/8)*(h + 1)*(g + h + r - 2), (g/8 - 1/8)*(h + 1)*(g - h - r + 2) - (g/8 - 1/8)*(h + 1)*(r + 1),                                -(g/4 - 1/4)*(g + 1)*(h - 1),                               (h/4 - 1/4)*(g + 1)*(h + 1),                               (g/4 - 1/4)*(g + 1)*(h + 1),                                -(h/4 - 1/4)*(g - 1)*(h + 1),                               (g/4 - 1/4)*(g + 1)*(h - 1),                                -(h/4 - 1/4)*(g + 1)*(h + 1),                                -(g/4 - 1/4)*(g + 1)*(h + 1),                               (h/4 - 1/4)*(g - 1)*(h + 1), - (r/4 - 1/4)*(g - 1)*(h - 1) - ((g - 1)*(h - 1)*(r + 1))/4, (r/4 - 1/4)*(g + 1)*(h - 1) + ((g + 1)*(h - 1)*(r + 1))/4, - (r/4 - 1/4)*(g + 1)*(h + 1) - ((g + 1)*(h + 1)*(r + 1))/4, (r/4 - 1/4)*(g - 1)*(h + 1) + ((g - 1)*(h + 1)*(r + 1))/4];
            % jacobian
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
            dH = J1*dHn;        % derivatives in physical coordinates,
            % 3x20 matrix, [dNi_dx; dNi_dy; dNi_dz]
            % with i = 1...20
            G = self.initialization.G;
            G(1:3,1:3:60) = dH;
            G(4:6,2:3:60) = dH;
            G(7:9,3:3:60) = dH;
        end
        
    end
    
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 20-NODED HEXAHEDRON
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.4 Solid isoparametric quadrilaterals and hexahedra
            g = X(1);
            h = X(2);
            r = X(3);
            N = [ -1/8*(1-g)*(1-h)*(1-r)*(2+g+h+r); % n1
                -1/8*(1+g)*(1-h)*(1-r)*(2-g+h+r);   % n2
                -1/8*(1+g)*(1+h)*(1-r)*(2-g-h+r);   % n3
                -1/8*(1-g)*(1+h)*(1-r)*(2+g-h+r);   % n4
                -1/8*(1-g)*(1-h)*(1+r)*(2+g+h-r);   % n5
                -1/8*(1+g)*(1-h)*(1+r)*(2-g+h-r);   % n6
                -1/8*(1+g)*(1+h)*(1+r)*(2-g-h-r);   % n7
                -1/8*(1-g)*(1+h)*(1+r)*(2+g-h-r);   % n8
                +1/4*(1-g)*(1+g)*(1-h)*(1-r);       % n9
                +1/4*(1-h)*(1+h)*(1+g)*(1-r);       % n10
                +1/4*(1-g)*(1+g)*(1+h)*(1-r);       % n11
                +1/4*(1-h)*(1+h)*(1-g)*(1-r);       % n12
                +1/4*(1-g)*(1+g)*(1-h)*(1+r);       % n13
                +1/4*(1-h)*(1+h)*(1+g)*(1+r);       % n14
                +1/4*(1-g)*(1+g)*(1+h)*(1+r);       % n15
                +1/4*(1-h)*(1+h)*(1-g)*(1+r);       % n16
                +1/4*(1-r)*(1+r)*(1-g)*(1-h);       % n17
                +1/4*(1-r)*(1+r)*(1+g)*(1-h);       % n18
                +1/4*(1-r)*(1+r)*(1+g)*(1+h);       % n19
                +1/4*(1-r)*(1+r)*(1-g)*(1+h)];      % n20
        end
        
        function X = natural_coordinates
            X = [ ...
                -1 -1 -1
                +1 -1 -1
                +1 +1 -1
                -1 +1 -1
                -1 -1 +1
                +1 -1 +1
                +1 +1 +1
                -1 +1 +1
                 0 -1 -1
                +1  0 -1
                 0 +1 -1
                -1  0 -1
                 0 -1 +1
                +1  0 +1
                 0 +1 +1
                -1  0 +1
                -1 -1  0
                +1 -1  0
                +1 +1  0
                -1 +1  0];
        end
        
    end
    
end % classdef

