classdef Wed15Element < ContinuumElement
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 3     % number of DOFs per node
        nNodes = 15         % number of nodes per element
        nDim = 3            % number of dimensions in local coordinates
        nelDOFs
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType = 'WED'
    end
    
    properties
        thickness = 1       % element thickness  
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Wed15Element(Material,Ngauss)
            % _____________________________________________________________
            %
            % SELF-FUNCTION
            % self = Wed15Element(Material,Ngauss)
            % defines element's properties
            %______________________________________________________________
            self.Material = Material;
            if nargin == 1
                Ngauss.lin = 3;  % Line integration order (1 to 5)
                Ngauss.tri = 03; % Triangle integration order (1,3,6,7,12)
            end
            self.thickness = 1;
            self.nelDOFs = self.nNodes * self.nDOFPerNode;
            ContinuumElementConstructor(self, Material, Ngauss);
        end
        
        function [G,detJ,dH] = shape_function_derivatives(self, X)
            %______________________________________________________________
            %
            % [G,detJ,dH] = G_WED15(self,g,h,r)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy uz vx vy vz wx wy wz}' (ux=du/dx...)
            % and p={u1,v1,w1,...,u15,v15,w15}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            g = X(1);
            h = X(2);
            r = X(3);
            xyz = self.nodes;
            % shape function derivatives in natural coordinates
            dHn = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            dHn(1,1)=1/2 - (r - 1)*(g + h - 1) - r^2/2 - ((r - 1)*(2*g + 2*h - 1))/2;
            dHn(1,2)=r^2/2 - g*(r - 1) - ((2*g - 1)*(r - 1))/2 - 1/2;
            dHn(1,4)=((r + 1)*(2*g + 2*h - 1))/2 + (r + 1)*(g + h - 1) - r^2/2 + 1/2;
            dHn(1,5)=((2*g - 1)*(r + 1))/2 + g*(r + 1) + r^2/2 - 1/2;
            dHn(1,7)=(r - 1)*(2*g + 2*h - 2) + 2*g*(r - 1);
            dHn(1,8)=-2*h*(r - 1);
            dHn(1,9)=2*h*(r - 1);
            dHn(1,10)=- (r + 1)*(2*g + 2*h - 2) - 2*g*(r + 1);
            dHn(1,11)=2*h*(r + 1);
            dHn(1,12)=-2*h*(r + 1);
            dHn(1,13)=r^2 - 1;
            dHn(1,14)=1 - r^2;
            dHn(2,1)=1/2 - (r - 1)*(g + h - 1) - r^2/2 - ((r - 1)*(2*g + 2*h - 1))/2;
            dHn(2,3)=r^2/2 - h*(r - 1) - ((2*h - 1)*(r - 1))/2 - 1/2;
            dHn(2,4)=((r + 1)*(2*g + 2*h - 1))/2 + (r + 1)*(g + h - 1) - r^2/2 + 1/2;
            dHn(2,6)=((2*h - 1)*(r + 1))/2 + h*(r + 1) + r^2/2 - 1/2;
            dHn(2,7)=2*g*(r - 1);
            dHn(2,8)=-2*g*(r - 1);
            dHn(2,9)=2*(r - 1)*(g + h - 1) + 2*h*(r - 1);
            dHn(2,10)=-2*g*(r + 1);
            dHn(2,11)=2*g*(r + 1);
            dHn(2,12)=- 2*(r + 1)*(g + h - 1) - 2*h*(r + 1);
            dHn(2,13)=r^2 - 1;
            dHn(2,15)=1 - r^2;
            dHn(3,1)=- ((g + h - 1)*(2*g + 2*h - 1))/2 - r*(g + h - 1);
            dHn(3,2)=g*r - (g*(2*g - 1))/2;
            dHn(3,3)=h*r - (h*(2*h - 1))/2;
            dHn(3,4)=((g + h - 1)*(2*g + 2*h - 1))/2 - r*(g + h - 1);
            dHn(3,5)=g*r + (g*(2*g - 1))/2;
            dHn(3,6)=h*r + (h*(2*h - 1))/2;
            dHn(3,7)=g*(2*g + 2*h - 2);
            dHn(3,8)=-2*g*h;
            dHn(3,9)=2*h*(g + h - 1);
            dHn(3,10)=-g*(2*g + 2*h - 2);
            dHn(3,11)=2*g*h;
            dHn(3,12)=-2*h*(g + h - 1);
            dHn(3,13)=2*r*(g + h - 1);
            dHn(3,14)=-2*g*r;
            dHn(3,15)=-2*h*r;
            
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
            % 3x10 matrix, [dNi_dx; dNi_dy; dNi_dz]
            % with i = 1...10
            G = self.initialization.G;
            G(1:3,1:3:45) = dH;
            G(4:6,2:3:45) = dH;
            G(7:9,3:3:45) = dH;
        end
        
    end
    
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 15-NODED WEDGE ELEMENT
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.6 Triangular, tetrahedral, and wedge elements
            g = X(1);
            h = X(2);
            r = X(3);
            N = [...
                    1/2*((1-g-h)*(2*(1-g-h)-1)*(1-r)-(1-g-h)*(1-r^2)); % u1
                    1/2*(g*(2*g-1)*(1-r)-g*(1-r^2)); % u2
                    1/2*(h*(2*h-1)*(1-r)-h*(1-r^2)); % u3
                    1/2*((1-g-h)*(2*(1-g-h)-1)*(1+r)-(1-g-h)*(1-r^2)); % u4
                    1/2*(g*(2*g-1)*(1+r)-g*(1-r^2)); % u5
                    1/2*(h*(2*h-1)*(1+r)-h*(1-r^2)); % u6
                    2*(1-g-h)*g*(1-r);  % u7
                    2*g*h*(1-r);        % u8
                    2*h*(1-g-h)*(1-r);  % u9
                    2*(1-g-h)*g*(1+r);  % u10
                    2*g*h*(1+r);        % u11
                    2*h*(1-g-h)*(1+r);  % u12
                    (1-g-h)*(1-r^2);    % u13
                    g*(1-r^2);          % u14
                    h*(1-r^2);          % u15
                ];
        end
        
        function X = natural_coordinates
            X = [ ...
                0   0   -1  % node 1 (corner)
                1   0   -1  % node 2 (corner)
                0   1   -1  % node 3 (corner)
                0   0   1   % node 4 (corner)
                1   0   1   % node 5 (corner)
                0   1   1   % node 6 (corner)
                .5  0   -1  % node 7
                .5  .5  -1  % node 8
                0   .5  -1  % node 9
                .5  0   1   % node 10
                .5  .5  1   % node 11
                0   .5  1   % node 12
                0   0   0   % node 13
                1   0   0   % node 14
                0   1   0]; % node 15
        end
        
    end
        
end % classdef
    
