classdef Quad4Element < ContinuumElement
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 2     % number of DOFs per node
        nNodes = 4          % number of nodes per element
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
        
        function self = Quad4Element(thickness, Material, Ngauss)
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
            % and p={u1,v1,...,u4,v4}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            r = X(1);
            s = X(2);
            xy = self.nodes;
            % shape function derivatives in natural coordinates
            dHn = 1/4*[ s-1,  r-1; 
                        1-s, -r-1; 
                        s+1,  r+1; 
                       -s-1,  1-r].';
            J = dHn*xy;
            detJ = det(J);
            dH = J\dHn;	% derivatives in physical coordinates,
                      	% 2x8 matrix, [dNi_dx; dNi_dy]
                       	% with i = 1...8
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
            N = 1/4*[ +(g - 1)*(h - 1); 
                      -(g + 1)*(h - 1);
                      +(g + 1)*(h + 1); 
                      -(g - 1)*(h + 1)];
        end
        
        function X = natural_coordinates
            X = [-1  -1   % node 1
                  1  -1   % node 2
                  1   1	  % node 3
                 -1   1]; % node 4
        end
        
    end
        
end % classdef
    
