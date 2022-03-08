classdef Hex20Element < Element
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 3     % number of DOFs per node
        nNodes = 20         % number of nodes per element
        nDim = 3            % number of dimensions in local coordinates
    end
    
    properties
        quadrature              % weights and points for gauss quadrature
        Material                % Object of class Material
        initialization          % some 0-matrices to speedup numerical integration
    end
    
    properties (Dependent)
        uniformBodyForce
        vol                     % volume of the element
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Hex20Element(Material,Ngauss)
            % _____________________________________________________________
            %
            % SELF-FUNCTION
            % self = Hex20Element(Material,Ngauss)
            % defines element's properties
            %______________________________________________________________
            if nargin == 1
                self.Material = Material;
                Ngauss = 2;
                [x,w]=lgwt(Ngauss,-1,1);
                self.quadrature.Ng = Ngauss;
                self.quadrature.X = x;	% gauss integration points
                self.quadrature.W = w;	% gauss integration weights
            elseif nargin==2
                self.Material = Material;
                [x,w]=lgwt(Ngauss,-1,1);
                self.quadrature.Ng = Ngauss;
                self.quadrature.X = x;
                self.quadrature.W = w;
            end
            % INIZIALIZATION of some matrices (this should speedup
            % numerical integration)
            v = self.Material.POISSONS_RATIO;
            E = self.Material.YOUNGS_MODULUS;
            C = [1-v v v 0 0 0; % isotropic material assumption
                v 1-v v 0 0 0;
                v v 1-v 0 0 0;
                0 0 0 .5*(1-2*v) 0 0;
                0 0 0 0 .5*(1-2*v) 0;
                0 0 0 0 0 .5*(1-2*v)]*E/((1+v)*(1-2*v));
            H = [1 0 0 0 0 0 0 0 0;
                0 0 0 0 1 0 0 0 0;
                0 0 0 0 0 0 0 0 1;
                0 1 0 1 0 0 0 0 0;
                0 0 1 0 0 0 1 0 0;
                0 0 0 0 0 1 0 1 0];
            self.initialization.A = zeros(6,9); % nonlinear strain
            self.initialization.G = zeros(9,60);% shape function derivatives
            self.initialization.Z = zeros(20);  % zero-matrix
            self.initialization.K = zeros(60);  % stiffness-element matrix
            self.initialization.F = zeros(60,1);% internal forces (element)
            self.initialization.C = C;          % constitutive law matrix
            self.initialization.H = H;          % linear strain
        end
        
        function Mel = mass_matrix(self)
            % _____________________________________________________________
            %
            % Mel = mass_matrix_global(self,nodes,~);
            % Mel: element-level mass matrix (in global coordinates)
            % nodes: [x1 y1 z1; ... ; x20 y20 z20] (coordinates)
            %______________________________________________________________
            X = self.quadrature.X;
            W = self.quadrature.W;
            rho = self.Material.DENSITY;
            Mel = zeros(60);
            for ii = 1:self.quadrature.Ng
                for jj = 1:self.quadrature.Ng
                    for kk = 1:self.quadrature.Ng
                        g = X(ii);
                        h = X(jj);
                        r = X(kk);
                        % shape functions and detJ (ABAQUS ORDER)
                        N = [ -1/8*(1-g)*(1-h)*(1-r)*(2+g+h+r); % n1
                            -1/8*(1+g)*(1-h)*(1-r)*(2-g+h+r); % n2
                            -1/8*(1+g)*(1+h)*(1-r)*(2-g-h+r); % n3
                            -1/8*(1-g)*(1+h)*(1-r)*(2+g-h+r); % n4
                            -1/8*(1-g)*(1-h)*(1+r)*(2+g+h-r); % n5
                            -1/8*(1+g)*(1-h)*(1+r)*(2-g+h-r); % n6
                            -1/8*(1+g)*(1+h)*(1+r)*(2-g-h-r); % n7
                            -1/8*(1-g)*(1+h)*(1+r)*(2+g-h-r); % n8
                            +1/4*(1-g)*(1+g)*(1-h)*(1-r);     % n9
                            +1/4*(1-h)*(1+h)*(1+g)*(1-r);     % n10
                            +1/4*(1-g)*(1+g)*(1+h)*(1-r);     % n11
                            +1/4*(1-h)*(1+h)*(1-g)*(1-r);     % n12
                            +1/4*(1-g)*(1+g)*(1-h)*(1+r);     % n13
                            +1/4*(1-h)*(1+h)*(1+g)*(1+r);     % n14
                            +1/4*(1-g)*(1+g)*(1+h)*(1+r);     % n15
                            +1/4*(1-h)*(1+h)*(1-g)*(1+r);     % n16
                            +1/4*(1-r)*(1+r)*(1-g)*(1-h);     % n17
                            +1/4*(1-r)*(1+r)*(1+g)*(1-h);     % n18
                            +1/4*(1-r)*(1+r)*(1+g)*(1+h);     % n19
                            +1/4*(1-r)*(1+r)*(1-g)*(1+h)];    % n20
                        [~,detJ] = G_HEX20(self,g,h,r); % <------- FUNCTION
                        NN(1,1:3:60) = N;
                        NN(2,2:3:60) = N;
                        NN(3,3:3:60) = N;
                        % integration of K and M through GAUSS QUADRATURE
                        Mel = Mel + W(ii)*W(jj)*W(kk)*(NN'*NN)*detJ;
                    end
                end
            end
            Mel = sparse(rho*Mel);
        end
        
        function [K,F] = tangent_stiffness_and_force(self,x)
            displ = self.extract_element_data(x);
            X = self.quadrature.X;
            W = self.quadrature.W;
            K = self.initialization.K;
            F = self.initialization.F;
            C = self.initialization.C;
            H = self.initialization.H;
            ZZ = self.initialization.Z;
            for ii = 1:self.quadrature.Ng
                for jj = 1:self.quadrature.Ng
                    for kk = 1:self.quadrature.Ng
                        g = X(ii);              % natural coordinates
                        h = X(jj);              % natural coordinates
                        r = X(kk);              % natural coordinates
                        we = W(ii)*W(jj)*W(kk); % weights
                        [G,detJ,dH] = G_HEX20(self,g,h,r); % <---- FUNCTION
                        th  = G*displ;
                        A = self.initialization.A;
                        A(1,1)=th(1); A(4,1)=th(2); A(5,1)=th(3); A(2,2)=th(2); A(4,2)=th(1);
                        A(6,2)=th(3); A(3,3)=th(3); A(5,3)=th(1); A(6,3)=th(2); A(1,4)=th(4);
                        A(4,4)=th(5); A(5,4)=th(6); A(2,5)=th(5); A(4,5)=th(4); A(6,5)=th(6);
                        A(3,6)=th(6); A(5,6)=th(4); A(6,6)=th(5); A(1,7)=th(7); A(4,7)=th(8);
                        A(5,7)=th(9); A(2,8)=th(8); A(4,8)=th(7); A(6,8)=th(9); A(3,9)=th(9);
                        A(5,9)=th(7); A(6,9)=th(8);
                        % Green Strain tensor
                        E = (H + 1/2*A)*th;
                        % second Piola-Kirchhoff stress tensor
                        s = C*E; % s = [S11 S22 S33 S12 S13 S23]
                        S = [s(1) s(4) s(5); s(4) s(2) s(6); s(5) s(6) s(3)];
                        Bnl = (H + A)*G;
                        % functions to integrate over volume
                        int_K1 = Bnl'*C*Bnl;
                        HSH = dH'*S*dH;
                        int_Ks = [HSH ZZ ZZ; ZZ HSH ZZ; ZZ ZZ HSH]; % (faster than blkdiag)
                        int_K = (int_K1 + int_Ks)*detJ;
                        int_F = (Bnl'*s)*detJ;
                        % integration of K and F through Gauss quadrature
                        K = K + we*int_K;
                        F = F + we*int_F;
                    end
                end
            end
        end
        
        function [K,F] = tangent_stiffness_and_force_local(self)
            % I have to include this to make the classdef happy...
            a=self.nDim; K=0*a; F=0;
            disp(' [tangent_stiffness_and_force_local]: I''m useless! :(')
            % ... but is useless for continuum elements (they're directly
            % computed in the global reference frame through isoparametric
            % mapping)
        end
        
        function xe = extract_element_data(self,x)
            % x is a vector of full DOFs
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            xe = x(index,:);
        end
        
        function  F = get.uniformBodyForce(self)
            % _____________________________________________________________
            %
            % F = uniform_body_force(self,direction)
            % This function computes a load along direction=3(Z) by
            % dividing the load on the 20 nodes according to the element
            % volume (V/20) [it might not be the best way, but still...]
            %______________________________________________________________
            F = sparse(60,1);
            F(3:3:end) = self.vol/20; % uniformly distributed pressure on the structure
        end
        
        function [T2, globalSubs] = T2(self)
            % this function computes the 3-tensor corresponding to the 
            % quadratic component of the nonlinear internal force in 
            % global coordinates at the element level.
                        
            % global DOFs associated to the element nodes
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = {index, index, index};
                        
            X = self.quadrature.X;
            W = self.quadrature.W;

            C = self.initialization.C;  % constitutive law matrix
            H = self.initialization.H;  % Linear strain matrix: eps_l = H*th

            % Quadratic strain matrix: A = L.th, eps_quad = A*th
            L = tenzeros([6,9,9]);
            L(1,1,1)=1; L(4,2,1)=1; L(5,3,1)=1; 
            L(4,1,2)=1; L(2,2,2)=1; L(6,3,2)=1; 
            L(5,1,3)=1; L(6,2,3)=1; L(3,3,3)=1;
            L(1,4,4)=1; L(4,5,4)=1; L(5,6,4)=1; 
            L(4,4,5)=1; L(2,5,5)=1; L(6,6,5)=1; 
            L(5,4,6)=1; L(6,5,6)=1; L(3,6,6)=1;
            L(1,7,7)=1; L(4,8,7)=1; L(5,9,7)=1; 
            L(4,7,8)=1; L(2,8,8)=1; L(6,9,8)=1; 
            L(5,7,9)=1; L(6,8,9)=1; L(3,9,9)=1;

            m = self.nNodes*self.nDOFPerNode;
            Q3h = tenzeros([m,m,m]);
            for ii = 1:self.quadrature.Ng
                for jj = 1:self.quadrature.Ng
                    for kk = 1:self.quadrature.Ng
                        g = X(ii);              % natural coordinates
                        h = X(jj);              % natural coordinates
                        r = X(kk);              % natural coordinates
                        we = W(ii)*W(jj)*W(kk); % weights
                        [G,detJ,~] = self.G_HEX20(g,h,r);

                        % G(x,y,z) and detJ from the position of the gauss points

                        %construct core part of the tensors for each gauss point
                        GHC = tensor((C*H*G)');
                        TG = tensor(G);  %create tensor object out of matrix
                        LGG = ttt(ttt(L,TG,3,1),TG,2,1);

                        Q3h_int = ttt(GHC,LGG,2,1);
                        Q3h = Q3h + Q3h_int*detJ*we;

                    end
                end
            end

            % build third order tensors using Q3h
            Q3ht = permute(Q3h,[3 2 1]);
            T2 = Q3h./2 + Q3ht;           
        end
        
        function [T3, globalSubs] = T3(self)
            % this function computes the 4-tensor corresponding to the 
            % quadratic component of the nonlinear internal force in 
            % global coordinates at the element level.
                        
            % global DOFs associated to the element nodes
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = cell(4,1);
            globalSubs(:) = {index};
                        
            X = self.quadrature.X;
            W = self.quadrature.W;

            C = self.initialization.C;  % constitutive law matrix

            % Quadratic strain matrix: A = L.th, eps_quad = A*th
            L = tenzeros([6,9,9]);
            L(1,1,1)=1; L(4,2,1)=1; L(5,3,1)=1; 
            L(4,1,2)=1; L(2,2,2)=1; L(6,3,2)=1; 
            L(5,1,3)=1; L(6,2,3)=1; L(3,3,3)=1;
            L(1,4,4)=1; L(4,5,4)=1; L(5,6,4)=1; 
            L(4,4,5)=1; L(2,5,5)=1; L(6,6,5)=1; 
            L(5,4,6)=1; L(6,5,6)=1; L(3,6,6)=1;
            L(1,7,7)=1; L(4,8,7)=1; L(5,9,7)=1; 
            L(4,7,8)=1; L(2,8,8)=1; L(6,9,8)=1; 
            L(5,7,9)=1; L(6,8,9)=1; L(3,9,9)=1;

            m = self.nNodes*self.nDOFPerNode;
            T3 = tenzeros([m,m,m,m]);

            for ii = 1:self.quadrature.Ng
                for jj = 1:self.quadrature.Ng
                    for kk = 1:self.quadrature.Ng
                        g = X(ii);              % natural coordinates
                        h = X(jj);              % natural coordinates
                        r = X(kk);              % natural coordinates
                        we = W(ii)*W(jj)*W(kk); % weights
                        [G,detJ,~] = self.G_HEX20(g,h,r);

                        % G(x,y,z) and detJ from the position of the gauss points

                        %construct core part of the tensors for each gauss point
                        TC = tensor(C);  %create tensor object, rename it to distinguish
                        TG = tensor(G);  %create tensor object out of matrix
                        LGG = ttt(ttt(L,TG,3,1),TG,2,1);

                        Q4h_int = ttt(ttt(permute(LGG,[2 1 3]),TC,2,1),LGG,3,1);
                        T3 = T3 + Q4h_int*detJ*we/2;

                    end
                end
            end
           
        end
        % ANCILLARY FUNCTIONS _____________________________________________
        
        function V = get.vol(self)
            xyz = self.nodes(1:8,:);
            x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
            dx = max(x)-min(x); dy = max(y)-min(y); dz = max(z)-min(z);
            alphaRadius = max([dx,dy,dz])+1;
            shp = alphaShape(x,y,z,alphaRadius);
            V = volume(shp);
        end
        
        function [G,detJ,dH] = G_HEX20(self,g,h,r)
            %______________________________________________________________
            %
            % [G,detJ,dH] = G_HEX20(self,g,h,r)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy uz vx vy vz wx wy wz}' (ux=du/dx...)
            % and p={u1,v1,w1,...,u20,v20,w20}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
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
        
    end % methods
    
    methods (Static)
        function [x,w]=lgwt(N,a,b)
            % _____________________________________________________________
            %
            % lgwt.m (downloaded from mathworks file exchange)
            %
            % This script is for computing definite integrals using
            % Legendre-Gauss quadrature. Computes the Legendre-Gauss nodes
            % and weights  on an interval [a,b] with truncation order N
            %
            % Suppose you have a continuous function f(x) which is defined
            % on [a,b] which you can evaluate at any x in [a,b]. Simply
            % evaluate it at all of the values contained in the x vector to
            % obtain a vector f. Then compute the definite integral using
            % sum(f.*w);
            %
            % Written by Greg von Winckel - 02/25/2004
            % _____________________________________________________________
            N=N-1;
            N1=N+1; N2=N+2;
            xu=linspace(-1,1,N1)';
            % Initial guess
            y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
            L=zeros(N1,N2);  % Legendre-Gauss Vandermonde Matrix
            Lp=zeros(N1,N2); % Derivative of LGVM
            % Compute the zeros of the N+1 Legendre Polynomial
            % using the recursion relation and the Newton-Raphson method
            y0=2;
            % Iterate until new points are uniformly within epsilon of old points
            while max(abs(y-y0))>eps
                L(:,1)=1;
                Lp(:,1)=0;
                L(:,2)=y;
                Lp(:,2)=1;
                for k=2:N1
                    L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
                end
                Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
                y0=y;
                y=y0-L(:,N2)./Lp;
            end
            x=(a*(1-y)+b*(1+y))/2;  % Linear map from[-1,1] to [a,b]
            w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;   % Compute the weights
        end
    end
    
end % classdef

