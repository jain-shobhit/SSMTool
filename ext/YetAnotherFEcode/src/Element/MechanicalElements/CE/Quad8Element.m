classdef Quad8Element < ContinuumElement
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 2     % number of DOFs per node
        nNodes = 8          % number of nodes per element
        nDim = 2            % number of dimensions in local coordinates
        nelDOFs
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType = 'QUAD'
    end
    
    properties
        thickness = 0       % element thickness, by default zero    
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Quad8Element(thickness, Material, Ngauss)
            % Self function (constructor)
            if nargin == 2
                Ngauss = 2;
            end
            self.thickness = thickness;
            self.nelDOFs = self.nNodes * self.nDOFPerNode;
            ContinuumElementConstructor(self, Material, Ngauss);
        end
        
        function [G,detJ,dH] = shape_function_derivatives(self, X)
            %______________________________________________________________
            %
            % [G,detJ,dH] = shape_function_derivatives(self, X)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy vx vy}' (ux=du/dx...)
            % and p={u1,v1,...,u8,v8}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            g = X(1);
            h = X(2);
            xy = self.nodes;
            % shape function derivatives in natural coordinates
            dHn = [ ...
                    -((h - 1)*(g + h + 1))/4 - ((g - 1)*(h - 1))/4, -((g - 1)*(g + h + 1))/4 - ((g - 1)*(h - 1))/4;
                    ((h - 1)*(h - g + 1))/4 - ((g + 1)*(h - 1))/4,   ((g + 1)*(h - g + 1))/4 + ((g + 1)*(h - 1))/4;
                    ((h + 1)*(g + h - 1))/4 + ((g + 1)*(h + 1))/4,   ((g + 1)*(g + h - 1))/4 + ((g + 1)*(h + 1))/4;
                    ((h + 1)*(g - h + 1))/4 + ((g - 1)*(h + 1))/4,   ((g - 1)*(g - h + 1))/4 - ((g - 1)*(h + 1))/4;
                                                        g*(h - 1),                                     g^2/2 - 1/2;
                                                      1/2 - h^2/2,                                      -h*(g + 1);
                                                       -g*(h + 1),                                     1/2 - g^2/2;
                                                      h^2/2 - 1/2,                                       h*(g - 1)]';
            J = dHn*xy;
            detJ = det(J);
            dH = J\dHn;	% derivatives in physical coordinates,
                      	% 2x16 matrix, [dNi_dx; dNi_dy]
                       	% with i = 1...10
            G = self.initialization.G;
            G(1:2,1:2:end) = dH;
            G(3:4,2:2:end) = dH;
        end
        
    end
    
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 8-NODED QUADRILATERAL
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.4 Solid isoparametric quadrilaterals and hexahedra
            g = X(1);
            h = X(2);
            N = 1/4*[...
                    -(1-g)*(1-h)*(1+g+h); 
                    -(1+g)*(1-h)*(1-g+h);
                    -(1+g)*(1+h)*(1-g-h); 
                    -(1-g)*(1+h)*(1+g-h);
                    2*(1-g)*(1+g)*(1-h);  
                    2*(1-h)*(1+h)*(1+g);
                    2*(1-g)*(1+g)*(1+h);  
                    2*(1-h)*(1+h)*(1-g)];
        end
        
        function X = natural_coordinates
            X = [ ...
                -1  -1  % node 1 (corner)
                1   -1  % node 2 (corner)
                1   1	% node 3 (corner)
                -1  1   % node 4 (corner)
                0   -1  % node 5
                1   0   % node 6
                0   1   % node 7
                -1  0]; % node 8
        end
        
    end
        
end % classdef
    
