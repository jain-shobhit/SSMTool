classdef ContinuumElement < Element
    % ContinuumElement collects all the functions shared by all the
    % elements based on a continuum deformation model. It includes both
    % planar and solid elements, such as quadrilaterals, hexahedra, 
    % tetrahedra and wedges.
    
    properties (Abstract)
        nodes               % global coordinates of element nodes
        nodeIDs             % the index location of element nodes
        nDOFPerNode         % number of DOFs per node
        nNodes              % number of nodes per element
        nDim                % number of dimensions in local coordinates
        nelDOFs             % number of DOFs per element
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType              % 'HEX'/'QUAD'/'TET'/'TRI'/'WED'
                            % --> see method "quadrature_rule"!
    end
    
    properties (Dependent)
        uniformBodyForce
        area                % area of a 2D element (NaN otherwise)
        vol                	% volume of a 3D element (NaN otherwise)
    end
    
    methods (Abstract) % function that every continuum element must implement
        [G,detJ,dH] = shape_function_derivatives(self, X)
        N = shape_functions(X)
        X = natural_coordinates
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function ContinuumElementConstructor(self, Material, Ngauss)
            % _____________________________________________________________
            %
            % ContinuumElementConstructor(self, Material, Ngauss)
            % defines element's properties
            %______________________________________________________________
            self.Material = Material;
            ne = self.nelDOFs;
            nd = self.nDim;
            
            quadrature_rule(self, Ngauss);
            
            % INIZIALIZATION of some matrices (this should speedup
            % numerical integration)
            if nd == 2
                thickness = self.thickness;
                C = self.Material.get_stress_strain_matrix_2D * thickness;
                H = [1 0 0 0; 
                    0 0 0 1; 
                    0 1 1 0];
                % Afun: (th) vector -> matrix (A)
                self.initialization.Afun = @(th) ... % nonlinear strain
                	 [th(1)	0     th(3) 0;
                      0    	th(2) 0     th(4);
                      th(2)	th(1) th(4) th(3)];
                % Sfun: (s)  tensor -> voigt (S)
             	self.initialization.Sfun = @(s) [s(1), s(3); 
                                                 s(3), s(2)];
                % shape function derivatives
                self.initialization.G = zeros(4,ne);
            elseif nd == 3
                self.thickness = 1;
                C = self.Material.get_stress_strain_matrix_3D;
                H = [1 0 0 0 0 0 0 0 0;
                    0 0 0 0 1 0 0 0 0;
                    0 0 0 0 0 0 0 0 1;
                    0 1 0 1 0 0 0 0 0;
                    0 0 1 0 0 0 1 0 0;
                    0 0 0 0 0 1 0 1 0];
                % Afun: (th) vector -> matrix (A)
                self.initialization.Afun = @(th) ... % nonlinear strain
                	[th(1) 0     0     th(4) 0     0     th(7) 0     0;
                     0     th(2) 0     0     th(5) 0     0     th(8) 0;
                     0     0     th(3) 0     0     th(6) 0     0     th(9);
                     th(2) th(1) 0     th(5) th(4) 0     th(8) th(7) 0;
                     th(3)     0 th(1) th(6) 0     th(4) th(9) 0     th(7);
                     0     th(3) th(2) 0     th(6) th(5) 0     th(9) th(8)];
                % Sfun: (s)  tensor -> voigt (S)
             	self.initialization.Sfun = @(s) [s(1) s(4) s(5); 
                                                 s(4) s(2) s(6); 
                                                 s(5) s(6) s(3)];
                % shape function derivatives
                self.initialization.G = zeros(9,ne);
            end
            self.initialization.K = zeros(ne); 	 % stiffness-element matrix
            self.initialization.F = zeros(ne,1); % element internal forces
            self.initialization.C = C;        	 % constitutive law matrix
            self.initialization.H = H;        	 % linear strain
        end
        
        function quadrature_rule(self, Ngauss)
            DIM = self.nDim;
            ELTYPE = self.elType;
            switch ELTYPE
                case {'HEX','QUAD'}
                    [x,w]=lgwt(Ngauss,-1,1);
                    X = zeros(DIM, Ngauss^DIM);
                    W = zeros(Ngauss^DIM, 1);
                    cont = 1;
                    for ii = 1:Ngauss
                        for jj = 1:Ngauss
                            if DIM == 3
                                for kk = 1:Ngauss
                                    X(:,cont) = [x(ii) x(jj) x(kk)].';
                                    W(cont) = w(ii)*w(jj)*w(kk);
                                    cont = cont+1;
                                end
                            elseif DIM == 2
                                X(:,cont) = [x(ii) x(jj)].';
                                W(cont) = w(ii)*w(jj);
                                cont = cont+1;
                            end
                        end
                    end
                    self.quadrature.Ng = Ngauss;
                    self.quadrature.X = X;	% gauss integration points
                    self.quadrature.W = W;	% gauss integration weights
                case {'TET'}
                    [x, w] = inttet(Ngauss);
                    self.quadrature.Ng = Ngauss;
                    self.quadrature.X = x;	% gauss integration points
                    self.quadrature.W = w;	% gauss integration weights
                case {'WED'}
                    [w, x] = wedge_rule (Ngauss.lin, Ngauss.tri);
                    self.quadrature.Ng = Ngauss;
                    self.quadrature.X = x;	% gauss integration points
                    self.quadrature.W = w;	% gauss integration weights
                otherwise
                    error([' No quadrature rule found for the element "' ...
                        ELTYPE '". Accepted values are: HEX/QUAD, TET/TRI, WED'])
            end
        end
        
        function Mel = mass_matrix(self)
            % _____________________________________________________________
            %
            % Mel = mass_matrix_global(self,nodes,~);
            % Mel: element-level mass matrix (in global coordinates)
            %______________________________________________________________
            X = self.quadrature.X;
            W = self.quadrature.W;
            rho = self.Material.DENSITY;
            Mel = zeros( self.nelDOFs );
            for ii = 1:length(W)
                Xi = X(:, ii);  % quadrature points
                we = W(ii);     % quadrature weights
                N = self.shape_functions(Xi);
                NN = kron(N', eye(self.nDim));
                [~, detJ] = shape_function_derivatives(self, Xi);
                % integration of K and M through GAUSS QUADRATURE
                Mel = Mel + (NN'*NN)*(we*detJ);
            end
            Mel = sparse(Mel)*(rho*self.thickness);
        end
        
        function [K,F] = tangent_stiffness_and_force(self,x)
            displ = self.extract_element_data(x);
            X = self.quadrature.X;
            W = self.quadrature.W;
            K = self.initialization.K;
            F = self.initialization.F;
            C = self.initialization.C;
            H = self.initialization.H;
            Afun = self.initialization.Afun; % fun: (th) vector -> matrix (A)
            Sfun = self.initialization.Sfun; % fun: (s)  tensor -> voigt (S)
            for ii = 1:length(W)
                Xi = X(:,ii);   % quadrature points
                we = W(ii);     % quadrature weights
                [G,detJ,dH] = shape_function_derivatives(self, Xi);
                th  = G*displ;
                A = Afun(th);
                % Green Strain tensor
                E = (H + 1/2*A)*th;
                % second Piola-Kirchhoff stress tensor
                s = C*E;        % stress tensor
                S = Sfun(s);    % stress vector (voigt)
                Bnl = (H + A)*G;
                % functions to integrate over volume
                int_K1 = Bnl'*C*Bnl;
                HSH = dH'*S*dH;
                int_Ks = kron(HSH, eye(self.nDim));
                int_K = (int_K1 + int_Ks)*detJ;
                int_F = (Bnl'*s)*detJ;
                % integration of K and F through Gauss quadrature
                K = K + we*int_K;
                F = F + we*int_F;
            end
        end
        
        function xe = extract_element_data(self,x)
            % x is a vector of full DOFs
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            xe = x(index,:);
        end
        
        function F = get.uniformBodyForce(self)
            % _____________________________________________________________
            %
            % F = uniform_body_force(self,direction)
            % This function computes a load along direction=3(Z) (or =2, Y,
            % for 2D problems) by dividing the load on the 20 nodes 
            % according to the element volume (in 3D) or area (in 2D) [it 
            % might not be the best way, but still...]
            %______________________________________________________________
            if self.nDim == 2
                F = sparse(16,1);
                F(2:2:end) = self.area/self.nNodes; % uniformly distributed pressure on the structure
            elseif self.nDim == 3
                F = sparse(60,1);
                F(3:3:end) = self.vol/self.nNodes; % uniformly distributed pressure on the structure
            end
        end
        
        % ADVANCED FUNCTIONS ______________________________________________
        
        function [T2, globalSubs] = T2(self, varargin)
            % this function computes the 3-tensor corresponding to the 
            % quadratic component of the nonlinear internal force in 
            % global coordinates at the element level.
            
            if ~isempty(varargin)
                Ve = varargin{1};
                m = size(Ve, 2);
                Vflag = true;
            else
                m = self.nNodes*self.nDOFPerNode;
                Vflag = false;
            end
                        
            % global DOFs associated to the element nodes
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = {index, index, index};
                        
            X = self.quadrature.X;
            W = self.quadrature.W;

            C = self.initialization.C;  % constitutive law matrix
            H = self.initialization.H;  % Linear strain matrix: eps_l = H*th
            L = tensor(L_matrix(self));	% Quadratic strain matrix: A = L.th, eps_quad = A*th
            
            Q3h = tenzeros([m,m,m]);
            for ii = 1:length(self.quadrature.W)
                Xi = X(:,ii);   % quadrature points
                we = W(ii);     % quadrature weights
                [G,detJ,~] = shape_function_derivatives(self, Xi); %get shape function derivative
                % G(x,y,z) and detJ from the position of the gauss points
                if Vflag
                    G = G*Ve;
                end
                
                %construct core part of the tensors for each gauss point
                GHC = tensor((C*H*G)');
                TG = tensor(G);  %create tensor object out of matrix
                LGG = ttt(ttt(L,TG,3,1),TG,2,1);

                Q3h_int = ttt(GHC,LGG,2,1);                
                Q3h = Q3h + Q3h_int*detJ*we;        
            end
            
            % build third order tensors using Q3h
            Q3ht = permute(Q3h,[3 2 1]);
            T2 = Q3h./2 + Q3ht;
        end
        
        function f2 = F2(self,x,y)
            % this function computes the quadratic component of the
            % nonlinear internal force in global coordinates at the element
            % level.
            x_e = self.extract_element_data(x);
            y_e = self.extract_element_data(y);           
            
            T2e = self.T2();
            f2 = ttv(T2e,{x_e,y_e},[2,3]);
            f2 = f2.data;                
        end

        function Df2 = DF2(self,x,y)
            % this function computes the Jacobian of the quadratic component 
            % of the nonlinear internal force in global coordinates 
            % at the element level. The Jacobian is evaluated along the
            % direction (x) and acted on the vector y. Hence, the output is a vector.  
            
            x_e = self.extract_element_data(x);
            y_e = self.extract_element_data(y);            
            
            T2e = self.T2();
            T2e = T2e + permute(T2e,[1 3 2]);
            Df2 = ttv(T2e,{x_e,y_e},[2 3]);
            Df2 = Df2.data;
        end
        
        function [T3, globalSubs] = T3(self, varargin)
            % this function computes the 4-tensor corresponding to the 
            % quadratic component of the nonlinear internal force in 
            % global coordinates at the element level.
                  
            if ~isempty(varargin)
                Ve = varargin{1};
                m = size(Ve, 2);
                Vflag = true;
            else
                m = self.nNodes*self.nDOFPerNode;
                Vflag = false;
            end
            
            % global DOFs associated to the element nodes
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = cell(4,1);
            globalSubs(:) = {index};
                        
            X = self.quadrature.X;
            W = self.quadrature.W;

            C = self.initialization.C;  % constitutive law matrix
            L = tensor(L_matrix(self));	% Quadratic strain matrix: A = L.th, eps_quad = A*th
            
            T3 = tenzeros([m,m,m,m]);
            for ii = 1:length(self.quadrature.W)
                Xi = X(:,ii);   % quadrature points
                we = W(ii);     % quadrature weights
                [G,detJ,~] = shape_function_derivatives(self, Xi); %get shape function derivative
                % G(x,y,z) and detJ from the position of the gauss points
                if Vflag
                    G = G*Ve; %element-level projection
                end
                
                %construct core part of the tensors for each gauss point
                TC = tensor(C);  %create tensor object, rename it to distinguish
                TG = tensor(G);  %create tensor object out of matrix
                LGG = ttt(ttt(L,TG,3,1),TG,2,1);

                Q4h_int = ttt(ttt(permute(LGG,[2 1 3]),TC,2,1),LGG,3,1);                
                T3 = T3 + Q4h_int*detJ*we/2;
            end           
           
        end
        
        function f3 = F3(self,x,y,z)
            % this function computes the cubic component of the
            % nonlinear internal force in global coordinates at the element
            % level along the direction (x,y,z). 
            x_e = self.extract_element_data(x);
            y_e = self.extract_element_data(y);           
            z_e = self.extract_element_data(z);
            
            T3e = self.T3();
            f3 = ttv(T3e,{x_e,y_e,z_e},[2,3,4]);
            f3 = f3.data;                
        end
        
        function Df3 = DF3(self,x,y,z)
            % this function computes the Jacobian of the cubic component of the
            % nonlinear internal force in global coordinates at the element
            % level. The Jacobian is evaluated along the direction (x,y) 
            % and acted on the vector z. Hence, the output is a vector.  
            x_e = self.extract_element_data(x);
            y_e = self.extract_element_data(y);
            z_e = self.extract_element_data(z);
           
            T3e = self.T3();
            T3e = T3e + permute(T3e,[1 3 4 2]) + permute(T3e,[1 3 2 4]);
            Df3 = ttv(T3e, {x_e, y_e, z_e}, [2 3 4]);
            Df3 = Df3.data;
        end
        
         
        % ANCILLARY FUNCTIONS _____________________________________________
        
        function V = get.vol(self)
            % volume is given by the integral of detJ (jacobian from 
            % isoparametric to physical space) over the volume of the
            % isoparametric element
            if self.nDim == 3
                detJ = 0;
                W = self.quadrature.W;
                X = self.quadrature.X;
                for ii = 1 : length( W )
                    Xi = X(:,ii);
                    [~, detJ_i] = shape_function_derivatives(self, Xi);
                    detJ = detJ + detJ_i * W(ii);
                end
                V = detJ;
            else
                V = NaN;
                % Returned volume = NaN (this is not a 3D problem)
            end
        end
        
        function A = get.area(self)
            % Integrate detJ (jacobian from isoparametric to physical
            % coordinates) over the area to get A
            if self.nDim == 2
                detJ = 0;
                W = self.quadrature.W;
                X = self.quadrature.X;
                for ii = 1 : length( W )
                    Xi = X(:, ii);
                    [~, detJ_i] = shape_function_derivatives(self, Xi);
                    detJ = detJ + detJ_i * W(ii);
                end
                A = detJ;
            else
                A = NaN;
                % Returned area = NaN (this is not a 2D problem)
            end
        end
        
        function L = L_matrix(self)
            % L used to compute the quadratic strain matrix: 
            % A = L.th, eps_quad = A*th ("." is the contraction operation)
            if self.nDim == 2
                L = zeros([3,4,4]);
                L(1,1,1) = 1; L(3,2,1) = 1; L(3,1,2) = 1; L(2,2,2) = 1;
                L(1,3,3) = 1; L(3,4,3) = 1; L(3,3,4) = 1; L(2,4,4) = 1;
            elseif self.nDim == 3
                L = zeros([6,9,9]);
                L(1,1,1)=1; L(4,2,1)=1; L(5,3,1)=1; 
                L(4,1,2)=1; L(2,2,2)=1; L(6,3,2)=1; 
                L(5,1,3)=1; L(6,2,3)=1; L(3,3,3)=1;
                L(1,4,4)=1; L(4,5,4)=1; L(5,6,4)=1; 
                L(4,4,5)=1; L(2,5,5)=1; L(6,6,5)=1; 
                L(5,4,6)=1; L(6,5,6)=1; L(3,6,6)=1;
                L(1,7,7)=1; L(4,8,7)=1; L(5,9,7)=1; 
                L(4,7,8)=1; L(2,8,8)=1; L(6,9,8)=1; 
                L(5,7,9)=1; L(6,8,9)=1; L(3,9,9)=1;
            end
        end
        
    end % methods
    
    methods (Static)
        
        % (no static methods)
        
    end
    
end % classdef

