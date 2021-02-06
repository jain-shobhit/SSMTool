classdef TriShellElement < Element
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 6     % number of DOFs per node
        nNodes = 3          % number of nodes per element
        nDim = 2            % number of dimensions in local coordinates
    end
    
    properties
        thickness = 0       % element thickness, by default zero
        Material            % Object of class Material     
        IsUpdated  
        uniformBodyForce
        massMatrixLocal
        dampingMatrixLocal
        stiffnessMatrixLocal
        LocalQuantities
        transformationMatrix
        localNodes          % coodinates of the 3 nodes in local coodinates
        area                % surface area of the element        
    end
    
    methods
        
        function self = TriShellElement(t, Material)
            narginchk(2,2)
            self.thickness = t;
            self.Material = Material;            
        end
        
        function set.nodes(self,nodes)
            self.nodes = nodes;
            not_updated(self);
        end
        
        function not_updated(self)
            IsUpdated.uniformBodyForce = false;
            IsUpdated.massMatrixLocal = false;
            IsUpdated.transformationMatrix = false;
            IsUpdated.localNodes = false;
            IsUpdated.LocalQuantities = false;
            self.IsUpdated = IsUpdated;
        end
        
        function set.Material(self,Material)
            self.Material = Material;
            not_updated_LQM(self);
        end
        
        function set.thickness(self,t)
            self.thickness = t;
            not_updated_LQM(self);
        end
        
        function not_updated_LQM(self)
            self.IsUpdated.LocalQuantities = false;
            self.IsUpdated.massMatrixLocal = false;
        end
        
        %% GET methods: Define how properties are computed

        function A = get.area(self)
            D = [1, 1, 1;
                self.localNodes(:,1:2).'];
            A = det(D)/2;
        end
        
        function T = get.transformationMatrix(self)
            if ~self.IsUpdated.transformationMatrix
                update_transformation_matrix(self)
            end
            T = self.transformationMatrix;
        end        
        
        function update_transformation_matrix(self)
            T = speye(18,18);
            R = rotation_matrix(self.nodes);
            T(1:3,1:3) = R;
            T(4:6,4:6) = R;
            T(7:9,7:9) = R;
            T(10:12,10:12) = R;
            T(13:15,13:15) = R;
            T(16:18,16:18) = R;
            self.transformationMatrix = T;
            self.IsUpdated.transformationMatrix = true;
        end
        
        function localNodes = get.localNodes(self)
            if ~self.IsUpdated.localNodes
                update_local_nodes(self)
            end
            localNodes = self.localNodes;
        end
        
        function update_local_nodes(self)
            self.localNodes = self.nodes * self.transformationMatrix(1:3,1:3).';
            self.IsUpdated.localNodes = true;
        end
        
        function M = get.massMatrixLocal(self)
            if ~self.IsUpdated.massMatrixLocal
                update_mass_matrix_local(self)
            end
            M = self.massMatrixLocal;
        end
        
        function update_mass_matrix_local(self)
            % membrane contribution
            M = sparse(18,18);
            Mm = Allman_membrane_mass(self);
            M([1 2 6 7 8 12 13 14 18],[1 2 6 7 8 12 13 14 18]) = Mm;
            % bending contribution
            Mb =Allman_bending_mass(self);
            M([3 4 5 9 10 11 15 16 17],[3 4 5 9 10 11 15 16 17]) = Mb;
            self.massMatrixLocal = M;
            self.IsUpdated.massMatrixLocal = true;
        end
        
        function C = get.dampingMatrixLocal(self)
            % Assume proportional damping
            C = self.Material.DAMPING_MODULUS * self.LocalQuantities.Kt / self.Material.YOUNGS_MODULUS;
        end
        
        function K = get.stiffnessMatrixLocal(self)
            K = self.LocalQuantities.Kt;            
        end
        
        function LQ = get.LocalQuantities(self)
            if ~self.IsUpdated.LocalQuantities
                update_local_quantities(self)
            end
            LQ = self.LocalQuantities;
        end
        
        function update_local_quantities(self)
            E = self.Material.YOUNGS_MODULUS;
            h = self.thickness;
            nu = self.Material.POISSONS_RATIO;
            Am = h*E/(1-nu^2)*[ 1 nu        0;
                nu 1        0;
                0  0 (1-nu)/2];
            LQ.Am = Am;
            localnodes = self.localNodes;
            
            ab=1;
            b0=4/9;
            b = [1/12 5/12 1/2 0 1/3 -1/3 -1/12 -1/2 -5/12];
            
            % geometric quantities
            x = localnodes(:,1);
            y = localnodes(:,2);
            
            x12=x(1)-x(2); x23=x(2)-x(3); x31=x(3)-x(1); x21=-x12; x32=-x23; x13=-x31;
            y12=y(1)-y(2); y23=y(2)-y(3); y31=y(3)-y(1); y21=-y12; y32=-y23; y13=-y31;
            
            A = (y21*x13-x21*y13)/2; A2=2*A; A4=4*A;
            
            LQ.A = A;
            LQ.BL = sparse(3,18);
            LQ.BL(:,[1 2 6 7 8 12 13 14 18])= [[y23,0,x32];[0,x32,y23];[y23*(y13-y21),x32*(x31-x12),(x31*y13-x12*y21)*2]*ab/6; ...
                [y31,0,x13];[0,x13,y31];[y31*(y21-y32),x13*(x12-x23),(x12*y21-x23*y32)*2]*ab/6; ...
                [y12,0,x21];[0,x21,y12];[y12*(y32-y13),x21*(x23-x31),(x23*y32-x31*y13)*2]*ab/6]'/2/A;
            
            Kb = A*LQ.BL'*Am*LQ.BL;
            Tx=1/(2*A)*[-y32 -y13 -y21];
            Ty=1/(2*A)*[-x23 -x31 -x12];
            Bw=sparse(3,9);
            Bw(1,1)=1; Bw(2,4)=1; Bw(3,7)=1;
            
            LQ.Kxx=sparse(18,18); LQ.Kyy=sparse(18,18); LQ.Kxy=sparse(18,18);
            TBx = Tx*Bw;
            TBy = Ty*Bw;
            TByTBx = TBy'*TBx;
            LQ.Kxx([3 4 5 9 10 11 15 16 17],[3 4 5 9 10 11 15 16 17])=TBx'*TBx;
            LQ.Kyy([3 4 5 9 10 11 15 16 17],[3 4 5 9 10 11 15 16 17])=TBy'*TBy;
            LQ.Kxy([3 4 5 9 10 11 15 16 17],[3 4 5 9 10 11 15 16 17])= TByTBx + TByTBx';
            
            Bx = sparse(3,18); By = sparse(3,18);
            Bx(1,1)=1; Bx(2,7)=1; Bx(3,13)=1;
            By(1,2)=1; By(2,8)=1; By(3,14)=1;
            TxBx =  Tx*Bx;
            TyBx = Ty*Bx;
            TxBxTyBx = TxBx'*TyBx;
            TyBy = Ty*By;
            TxBy = Tx*By;
            TxByTyBy = TxBy'*TyBy;
            LQ.Kxy = LQ.Kxy + TxBxTyBx + TxBxTyBx' + TxByTyBy + TxByTyBy';
            
            % % hierarchical rotation to global displacement matrix
            Tqu=[[x32,y32,A4,x13,y13, 0,x21,y21, 0]; ...
                [x32,y32, 0,x13,y13,A4,x21,y21, 0]; ...
                [x32,y32, 0,x13,y13,0,x21,y21,A4]]/A4;
            
            % element sides length
            LL21=x21^2+y21^2; LL32=x32^2+y32^2; LL13=x13^2+y13^2;
            
            % cartesian to natural transformation matrix
            Te=[[y23*y13*LL21, y31*y21*LL32, y12*y32*LL13]; ...
                [x23*x13*LL21,x31*x21*LL32,x12*x32*LL13]; ...
                [(y23*x31+x32*y13)*LL21,(y31*x12+x13*y21)*LL32, (y12*x23+x21*y32)*LL13]]/(A*A4);
            
            % natural strain to hierarchical bending modes matrices
            Q1=[[b(1),b(2),b(3)]/LL21;[b(4),b(5),b(6)]/LL32;[b(7),b(8),b(9)]/LL13]*A2/3;
            Q2=[[b(9),b(7),b(8)]/LL21;[b(3),b(1),b(2)]/LL32;[b(6),b(4),b(5)]/LL13]*A2/3;
            Q3=[[b(5),b(6),b(4)]/LL21;[b(8),b(9),b(7)]/LL32;[b(2),b(3),b(1)]/LL13]*A2/3;
            Q4=(Q1+Q2)/2; Q5=(Q2+Q3)/2; Q6=(Q3+Q1)/2;
            
            % natural stress-strain matrix
            Enat=Te'*Am*Te;
            
            % higher order stiffnes matrix in local coordinates (hierarchical bending modes)
            Kq=(3/4)*b0*A*(Q4'*Enat*Q4+Q5'*Enat*Q5+Q6'*Enat*Q6);
            
            % high order stiffness
            Kh=Tqu'*Kq*Tqu;
            
            % tangential stiffness matrix
            
            Kt([3 4 5 9 10 11 15 16 17],[3 4 5 9 10 11 15 16 17]) =...
                Allman_bending_stiffness(self);
            % contribution for high order matrix
            Kt([1 2 6 7 8 12 13 14 18],[1 2 6 7 8 12 13 14 18])=Kh;
            LQ.Kt = Kt + Kb;
            self.LocalQuantities = LQ;
            self.IsUpdated.LocalQuantities = true;
        end

        function  f = get.uniformBodyForce(self)
            if ~self.IsUpdated.uniformBodyForce
                update_uniform_body_force(self);
            end
            f = self.uniformBodyForce;
        end
        
        function update_uniform_body_force(self)
            % this function computes the element stiffness matrix and
            % internal force vector in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            
            T_e = self.transformationMatrix;
            
            % uniformly distributed pressure on the structure 
            f_local = sparse(18,1);
            f_local(3:6:end) = self.area/3;           

            self.uniformBodyForce = T_e.' * f_local;
            self.IsUpdated.uniformBodyForce = true;
        end
        
        function [T2, globalSubs] = T2(self)
            % this function computes the 3-tensor corresponding to the 
            % quadratic component of the nonlinear internal force in 
            % global coordinates at the element level.
            
            T_e = self.transformationMatrix;
            
            % global DOFs associated to the element nodes
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = {index, index, index};
            
            % compute tensor in local coordinates
            T2e = self.T2_local();
            
            % rotate to global coordinates
            T2 = ttm(T2e, {T_e.',T_e.',T_e.'},[1,2,3]);
        end
        
        function [T3, globalSubs] = T3(self)
            % this function computes the 4-tensor corresponding to the 
            % quadratic component of the nonlinear internal force in 
            % global coordinates at the element level.
            
            T_e = self.transformationMatrix;
            
            % global DOFs associated to the element nodes
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = cell(4,1);
            globalSubs(:) = {index};
            
            % compute tensor in local coordinates           
            T3e = self.T3_local();
            
            % rotate to global coordinates
            T3 = ttm(T3e, {T_e.',T_e.',T_e.',T_e.'},[1,2,3,4]);
        end
        
        %% Global level methods
        function xe = extract_element_data(self,x)
            % x is a vector of full DOFs            
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            xe = x(index,:);
        end
        
        function q = transform_to_local(self,x,TYPE)
            T = self.transformationMatrix;
            switch TYPE
                case 'vector'
                    q = T * x;
                    
                case 'matrix'
                    q = T * x * T.';
            end
        end
        
        function x = transform_to_global(self,q,TYPE)
            T = self.transformationMatrix;
            switch TYPE
                case 'vector'
                    x = T.' * q;
                    
                case 'matrix'
                    x = T * q * T.';
            end
        end
            
        function M = mass_matrix(self)
            % this function computes the element mass matrix in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            T_e = self.transformationMatrix;
            M_local = self.massMatrixLocal; 
            M = T_e.' * M_local * T_e;
        end
        
        function C = damping_matrix(self)
            % this function computes the element mass matrix in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            T_e = self.transformationMatrix;
            C_local = self.dampingMatrixLocal; 
            C = T_e.' * C_local * T_e;
        end
        
        function K = stiffness_matrix(self)
            % this function computes the element mass matrix in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            T_e = self.transformationMatrix;
            K_local = self.stiffnessMatrixLocal; 
            K = T_e.' * K_local * T_e;
        end       
        
        
        function [K, F] = tangent_stiffness_and_force(self,x)         
            % this function computes the element stiffness matrix and
            % internal force vector in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            % x : vector of full DOFs in global coordinates    
            
            x_e = self.extract_element_data(x);
            T_e = self.transformationMatrix;
            % Displacements in local coordinates
            q = T_e*x_e; 
            % Forces and Stiffness in local coordinates            
            [K_local, F_local] = self.tangent_stiffness_and_force_local(q);
            % Forces and Stiffness in global coordinates            
            K = T_e.' * K_local * T_e;
            F = T_e.' * F_local;
        end  
        
        function [F] = internal_force(self,x)         
            % this function computes the element stiffness matrix and
            % internal force vector in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            % x : vector of full DOFs in global coordinates    
            
            x_e = self.extract_element_data(x);
            T_e = self.transformationMatrix;
            % Displacements in local coordinates
            q = T_e*x_e; 
            % Forces and Stiffness in local coordinates            
            [F_local] = self.internal_force_local(q);
            % Force in global coordinates            
            F = T_e.' * F_local;
        end
        
        function [F] = f2(self,x)         
            % this function computes the element stiffness matrix and
            % internal force vector in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            % x : vector of full DOFs in global coordinates    
            
            x_e = self.extract_element_data(x);
            T_e = self.transformationMatrix;
            % Displacements in local coordinates
            q = T_e*x_e; 
            % Forces and Stiffness in local coordinates            
            [F_local] = self.f2_local(q);
            % Force in global coordinates            
            F = T_e.' * F_local;
        end
        
        function [F] = f3(self,x)         
            % this function computes the element stiffness matrix and
            % internal force vector in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            % x : vector of full DOFs in global coordinates    
            
            x_e = self.extract_element_data(x);
            T_e = self.transformationMatrix;
            % Displacements in local coordinates
            q = T_e*x_e; 
            % Forces and Stiffness in local coordinates            
            [F_local] = self.f3_local(q);
            % Force in global coordinates            
            F = T_e.' * F_local;
        end
        
        function [K] = stiffness_derivative(self,x,v)
            % this function computes the element stiffness derivative matrix
            % in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            % x : element DOF values in global coordinates
            x_e = self.extract_element_data(x);
            v_e = self.extract_element_data(v);
            
            T_e = self.transformationMatrix;
            
            q = T_e * x_e;
            eta = T_e * v_e;
            [K_local] = self.stiffness_derivative_local(q,eta);
            K = Te.' * K_local * Te;
        end

        %% Local level methods
        function K = tangent_stiffness_matrix_local(self,x)
            LQ = self.LocalQuantities;
            % N
            BNL=[x.' * LQ.Kxx;
                x.' * LQ.Kyy;
                x.' * LQ.Kxy];
            
            N=LQ.Am * (LQ.BL + BNL/2) * x;
            K= LQ.Kt + LQ.A*(LQ.BL.' * LQ.Am * BNL + BNL.' * LQ.Am * LQ.BL + ...
                BNL.' * LQ.Am * BNL)+ LQ.A * (N(1) * LQ.Kxx + N(2) * LQ.Kyy +...
                N(3) * LQ.Kxy);
        end
        
        function F = internal_force_local(self,x)
            LQ = self.LocalQuantities;
            F = LQ.Kt*x;
            % N
            BNL=[x.' * LQ.Kxx;
                x.' * LQ.Kyy;
                x.' * LQ.Kxy];
            
            N=LQ.Am * (LQ.BL + BNL/2) * x;
            % Internal forces
            F= F + LQ.A * (LQ.BL + BNL).' * N - LQ.A*(LQ.BL.'*(LQ.Am*(LQ.BL*x)));
            
        end
        
        function [K,F] = tangent_stiffness_and_force_local(self,x)
            LQ = self.LocalQuantities;
            F = LQ.Kt*x;
            % N
            BNL=[x.' * LQ.Kxx;
                x.' * LQ.Kyy;
                x.' * LQ.Kxy];
            
            N=LQ.Am * (LQ.BL + BNL/2) * x;
            K= LQ.Kt + LQ.A*(LQ.BL.' * LQ.Am * BNL + BNL.' * LQ.Am * LQ.BL + ...
                BNL.' * LQ.Am * BNL)+ LQ.A * (N(1) * LQ.Kxx + N(2) * LQ.Kyy +...
                N(3) * LQ.Kxy);
            % Internal forces
            F= F + LQ.A * (LQ.BL + BNL).' * N - LQ.A*(LQ.BL.'*(LQ.Am*(LQ.BL*x)));
            
        end
        
        
        function [Kb]=Allman_bending_stiffness(self)
            % Ref: D.J. Allman "A Simple Cubic Displacement element for plate bending",
            % Int. J. Num. Meth. Eng., vol. 10, 263-281, 1976
            E = self.Material.YOUNGS_MODULUS;
            t = self.thickness;
            nu = self.Material.POISSONS_RATIO;
            nodes = self.localNodes; %#ok<*PROP>
            
            co1=nodes(1,:);
            co2=nodes(2,:);
            co3=nodes(3,:);
            
            x1=co1(1); y1=co1(2);
            x2=co2(1); y2=co2(2);
            x3=co3(1); y3=co3(2);
            
            A = self.area;
            
            x=[x1 x2 x3];
            y=[y1 y2 y3];
            
            x12=x1-x2;  y21=y2-y1;
            x23=x2-x3;  y32=y3-y2;
            x31=x3-x1;  y13=y1-y3;
            
            D=E*t^3/(12*(1-nu^2));
            
            L12=sqrt(x12^2+y21^2);
            L23=sqrt(x23^2+y32^2);
            L31=sqrt(x31^2+y13^2);
            
            
            sin12=x12/L12; cos12=y21/L12;
            sin23=x23/L23; cos23=y32/L23;
            sin31=x31/L31; cos31=y13/L31;
            
            si=[sin12 sin23 sin31];
            co=[cos12 cos23 cos31];
            
            H = zeros(7,7);
            
            H(1,:)=[                                     2;
                0;
                2*nu;
                2*x1+2*x2+2*x3;
                2/3*y1+2/3*y2+2/3*y3;
                2/3*nu*(x2-x3)+2/3*nu*(x1-x3)+2*nu*x3;
                2*nu*(y2-y3)+2*nu*(y1-y3)+6*nu*y3]';
            
            H(2,:)=[                                                         0;
                1-nu;
                0;
                0;
                1/6*(4-4*nu)*(x2-x3)+1/6*(4-4*nu)*(x1-x3)+1/2*(4-4*nu)*x3;
                1/6*(4-4*nu)*(y2-y3)+1/6*(4-4*nu)*(y1-y3)+1/2*(4-4*nu)*y3;
                0]';
            
            H(3,:)=[                                     0;
                0;
                2;
                2*nu*(x2-x3)+2*nu*(x1-x3)+6*nu*x3;
                2/3*nu*(y2-y3)+2/3*nu*(y1-y3)+2*nu*y3;
                2/3*x1+2/3*x2+2/3*x3;
                2*y1+2*y2+2*y3]';
            
            
            
            H(4,:)=[                                                                                                                                                                               0;
                0;
                0;
                3*(x2-x3)^2+3*(x1-x3)*(x2-x3)+3*(x1-x3)^2+12*x3*(x2-x3)+12*x3*(x1-x3)+18*x3^2;
                1/12*(12*x2-12*x3)*(y2-y3)+1/24*(12*x1-12*x3)*(y2-y3)+1/24*(12*x2-12*x3)*(y1-y3)+1/12*(12*x1-12*x3)*(y1-y3)+2*x3*(y2-y3)+1/6*(12*x2-12*x3)*y3+2*x3*(y1-y3)+1/6*(12*x1-12*x3)*y3+6*x3*y3;
                nu*(x2-x3)^2+nu*(x1-x3)*(x2-x3)+nu*(x1-x3)^2+4*nu*x3*(x2-x3)+4*nu*x3*(x1-x3)+6*nu*x3^2;
                3*nu*(x2-x3)*(y2-y3)+3/2*nu*(x1-x3)*(y2-y3)+3/2*nu*(x2-x3)*(y1-y3)+3*nu*(x1-x3)*(y1-y3)+6*nu*x3*(y2-y3)+6*nu*(x2-x3)*y3+6*nu*x3*(y1-y3)+6*nu*(x1-x3)*y3+18*nu*x3*y3]';
            
            H(5,:)=[                                                                                                                                                                                                                                  0;
                0;
                0;
                0;
                1/3*(8-8*nu)*x3*(x2-x3)+1/12*(8-8*nu)*(x1-x3)*(x2-x3)+1/3*(y2-y3)^2+1/12*(8-8*nu)*(x2-x3)^2+1/3*(8-8*nu)*x3*(x1-x3)+4/3*y3*(y2-y3)+1/3*(y1-y3)*(y2-y3)+2*y3^2+4/3*y3*(y1-y3)+1/3*(y1-y3)^2+1/2*(8-8*nu)*x3^2+1/12*(8-8*nu)*(x1-x3)^2;
                1/12*(8-4*nu)*(x2-x3)*(y2-y3)+1/24*(8-4*nu)*(x1-x3)*(y2-y3)+1/24*(8-4*nu)*(x2-x3)*(y1-y3)+1/12*(8-4*nu)*(x1-x3)*(y1-y3)+1/6*(8-4*nu)*x3*(y2-y3)+1/6*(8-4*nu)*(x2-x3)*y3+1/6*(8-4*nu)*x3*(y1-y3)+1/6*(8-4*nu)*(x1-x3)*y3+1/2*(8-4*nu)*x3*y3;
                nu*(y2-y3)^2+nu*(y1-y3)*(y2-y3)+nu*(y1-y3)^2+4*nu*y3*(y2-y3)+4*nu*y3*(y1-y3)+6*nu*y3^2]';
            
            
            
            H(6,:)=[                                                                                                                                                                                                                            0;
                0;
                0;
                0;
                0;
                1/12*(8-8*nu)*(y2-y3)^2+1/3*(8-8*nu)*y3*(y1-y3)+1/12*(8-8*nu)*(y1-y3)*(y2-y3)+1/2*(8-8*nu)*y3^2+1/12*(8-8*nu)*(y1-y3)^2+2*x3^2+4/3*x3*(x1-x3)+1/3*(8-8*nu)*y3*(y2-y3)+1/3*(x1-x3)^2+4/3*x3*(x2-x3)+1/3*(x1-x3)*(x2-x3)+1/3*(x2-x3)^2;
                1/12*(12*x2-12*x3)*(y2-y3)+1/24*(12*x1-12*x3)*(y2-y3)+1/24*(12*x2-12*x3)*(y1-y3)+1/12*(12*x1-12*x3)*(y1-y3)+2*x3*(y2-y3)+1/6*(12*x2-12*x3)*y3+2*x3*(y1-y3)+1/6*(12*x1-12*x3)*y3+6*x3*y3]';
            
            H(7,:)=[                                                                     0;
                0;
                0;
                0;
                0;
                0;
                3*(y2-y3)^2+3*(y1-y3)*(y2-y3)+3*(y1-y3)^2+12*y3*(y2-y3)+12*y3*(y1-y3)+18*y3^2]';
            
            H=2*A*D*(H'+H-diag(diag(H)));
            
            T=zeros(12,9);
            T(1,1)=1;
            T(2,4)=1;
            T(3,7)=1;
            
            T(4,1:6)=[L12/2 -L12^2/12*sin12 L12^2/12*cos12 L12/2 L12^2/12*sin12 -L12^2/12*cos12];
            T(5,4:end)=[L23/2 -L23^2/12*sin23 L23^2/12*cos23 L23/2 L23^2/12*sin23 -L23^2/12*cos23];
            T(6,[1:3 7:9])=[L31/2 L31^2/12*sin31 -L31^2/12*cos31 L31/2 -L31^2/12*sin31 L31^2/12*cos31];
            
            T(7,[2:3 5:6])=[-L12/3*cos12 -L12/3*sin12 -L12/6*cos12 -L12/6*sin12];
            T(8,[2:3 5:6])=[-L12/6*cos12 -L12/6*sin12 -L12/3*cos12 -L12/3*sin12];
            
            T(9,[5:6 8:9])=[-L23/3*cos23 -L23/3*sin23 -L23/6*cos23 -L23/6*sin23];
            T(10,[5:6 8:9])=[-L23/6*cos23 -L23/6*sin23 -L23/3*cos23 -L23/3*sin23];
            
            T(11,[2:3 8:9])=[-L31/6*cos31 -L31/6*sin31 -L31/3*cos31 -L31/3*sin31];
            T(12,[2:3 8:9])=[-L31/3*cos31 -L31/3*sin31 -L31/6*cos31 -L31/6*sin31];
            
            S1=2*sin12*cos12-2*sin31*cos31;
            S2=2*sin23*cos23-2*sin12*cos12;
            S3=2*sin31*cos31-2*sin23*cos23;
            C1=cos12^2-sin12^2-(cos31^2-sin31^2);
            C2=cos23^2-sin23^2-(cos12^2-sin12^2);
            C3=cos31^2-sin31^2-(cos23^2-sin23^2);
            
            S=[S1 S2 S3];
            C=[C1 C2 C3];
            
            B=zeros(7,12);
            
            B(1,1:3)=D*(1-nu)*S;
            B(2,1:3)=-D*(1-nu)*C;
            B(3,1:3)=-B(1,1:3);
            B(4,1:3)=3*D*(1-nu)*x.*S;
            B(5,1:3)=-D*(1-nu)*(-y.*S + 2*x.*C);
            B(6,1:3)=-D*(1-nu)*(x.*S + 2*y.*C);
            B(7,1:3)=-3*D*(1-nu)*y.*S;
           
            B(4,4:6) = -6*D*co.*(1+(1-nu)*si.^2);
            B(5,4:6) = -2*D*si.*((2-nu)*si.^2 + (2*nu-1) * co.^2);
            B(6,4:6) = -2*D*co.*((2-nu)*co.^2 + (2*nu-1) * si.^2);
            B(7,4:6) = -6*D*si.*(1 + (1-nu)*co.^2);
            
            B(1,7:2:11)= -2*D*(co.^2 + nu*si.^2);
            B(2,7:2:11)= -D*(1-nu)*2*si.*co;
            B(3,7:2:11)= -2*D*(nu*co.^2+si.^2);
            B(4,7:2:11)= -6*D*x.*(co.^2 + nu * si.^2);
            B(5,7:2:11)= -2*D*(y.*(co.^2 + nu*si.^2) + (1-nu)*2*x.*si.*co);
            B(6,7:2:11)= -2*D*(x.*(nu*co.^2 + si.^2) + (1-nu)*2*y.*si.*co);
            B(7,7:2:11)= -6*D*y.*(nu*co.^2 + si.^2);
            
            B(1,8:2:12) = -2*D*(co.^2 + nu*si.^2);
            B(2,8:2:12) = -D*(1-nu)*2*si.*co;
            B(3,8:2:12) = -2*D*(nu*co.^2+si.^2);
            B(4,8:2:12) = -6*D*x([2 3 1]).*(co.^2 + nu*si.^2);
            B(5,8:2:12) = -2*D*(y([2 3 1]).*(co.^2 + nu*si.^2) + ...
                (1-nu)*2*x([2 3 1]).*si.*co);
            B(6,8:2:12) = -2*D*(x([2 3 1]).*(nu*co.^2 + si.^2) + ...
                (1-nu)*2*y([2 3 1]).*si.*co);
            B(7,8:2:12) = -6*D*y([2 3 1]).*(nu*co.^2 + si.^2);
            
            
            Kb=(B*T)'*(H\(B*T));
            
            Kb([2 3 5 6 8 9],:)=Kb([3 2 6 5 9 8],:);
            Kb(:,[2 3 5 6 8 9])=Kb(:,[3 2 6 5 9 8]);
            
            Kb([3 6 9],:)=-Kb([3 6 9],:);
            Kb(:,[3 6 9])=-Kb(:,[3 6 9]);
            
            Kb=1/2*(Kb'+Kb);
        end
        
        function [Mmlu]=Allman_membrane_mass(self)
            
            % Allman triangle membrane consistent mass matrix
            %
            % Ref: D.J. Allman
            %"Implementation of a flat facet shell finite element for applcication in structural dynamics"
            % Computers & Structures, val 59, pp 657-663
            
            rho = self.Material.DENSITY;
            h = self.thickness;
            nodes = self.localNodes;
            
            
            co1=nodes(1,:);
            co2=nodes(2,:);
            co3=nodes(3,:);
            
            x1=co1(1); y1=co1(2);
            x2=co2(1); y2=co2(2);
            x3=co3(1); y3=co3(2);
            
            x12=x1-x2;  y21=y2-y1;
            x23=x2-x3;  y32=y3-y2;
            x31=x3-x1;  y13=y1-y3;
            
            A = self.area;
            
            X12=-(x1-x2)/(4*A);  Y12=-(y1-y2)/(4*A);
            X23=-(x2-x3)/(4*A);  Y23=-(y2-y3)/(4*A);
            X31=-(x3-x1)/(4*A);  Y31=-(y3-y1)/(4*A);
            
            
            Bu=zeros(3,9);
            Bu(1,1)=1; Bu(2,4)=1; Bu(3,7)=1;
            
            Bv=zeros(3,9);
            Bv(1,2)=1; Bv(2,5)=1; Bv(3,8)=1;
            
            Bau=zeros(6,9);
            
            Bau(1,[3 6])=1/2*[-y21 y21];
            Bau(2,[6 9])=1/2*[-y32 y32];
            Bau(3,[3 9])=1/2*[y13 -y13];
            Bau(4,1:8)=[X23*y21 Y23*y21 1/2*y21 X31*y21 Y31*y21 1/2*y21 X12*y21 Y12*y21];
            Bau(5,[1:2 4:9])=[X23*y32 Y23*y32 X31*y32 Y31*y32 1/2*y32 X12*y32 Y12*y32 1/2*y32];
            Bau(6,[1:5 7:9])=[X23*y13 Y23*y13 1/2*y13 X31*y13 Y31*y13 X12*y13 Y12*y13 1/2*y13];
            
            
            Bav=zeros(6,9);
            
            Bav(1,[3 6])=1/2*[-x12 x12];
            Bav(2,[6 9])=1/2*[-x23 x23];
            Bav(3,[3 9])=1/2*[x31 -x31];
            Bav(4,1:8)=[X23*x12 Y23*x12 1/2*x12 X31*x12 Y31*x12 1/2*x12 X12*x12 Y12*x12];
            Bav(5,[1:2 4:9])=[X23*x23 Y23*x23 X31*x23 Y31*x23 1/2*x23 X12*x23 Y12*x23 1/2*x23];
            Bav(6,[1:5 7:9])=[X23*x31 Y23*x31 1/2*x31 X31*x31 Y31*x31 X12*x31 Y12*x31 1/2*x31];
            
            
            intNN=[1/6 1/12 1/12;
                0   1/6   1/12;
                0     0   1/6];
            
            intNN=rho*h*A*1/2*(intNN'+intNN);
            
            intNNa=rho*h*A*[1/30 1/60 1/30 -1/180      0  -1/180;
                1/30 1/30 1/60  1/180 -1/180       0;
                1/60 1/30 1/30      0  1/180  -1/180];
            
            intNaN=intNNa';
            
            intNaNa=h*rho*A*[1/90 1/180 1/180       0 -1/1260  1/1260;
                0  1/90 1/180  1/1260       0 -1/1260;
                0     0  1/90 -1/1260  1/1260       0;
                0     0     0   1/840 -1/2520 -1/2520;
                0     0     0       0   1/840 -1/2520;
                0     0     0       0       0   1/840];
            
            intNaNa=1/2*(intNaNa+intNaNa');
            
            Mu=(Bu'*intNN*Bu+Bu'*intNNa*Bau+Bau'*intNaN*Bu+Bau'*intNaNa*Bau);
            
            Mv=(Bv'*intNN*Bv+Bv'*intNNa*Bav+Bav'*intNaN*Bv+Bav'*intNaNa*Bav);
            Mm=Mu+Mv;
            
            mass = A*h*rho;
            s = Mm(2,2) + Mm(5,5) + Mm(8,8);
            lumped = diag(Mm).*(mass./s);
            Mmlu = diag(lumped);
        end
        
        function [Mblu]=Allman_bending_mass(self)
            
            % Allman triangle bending consistent mass matrix
            %
            % Ref: D.J. Allman
            %"Implementation of a flat facet shell finite element for application in structural dynamics"
            % Computers & Structures, vol 59, pp 657-663
                        
            rho = self.Material.DENSITY;
            h = self.thickness;
            nodes = self.localNodes;

            co1=nodes(1,:);
            co2=nodes(2,:);
            co3=nodes(3,:);
            
            x1=co1(1); y1=co1(2);
            x2=co2(1); y2=co2(2);
            x3=co3(1); y3=co3(2);
            
            x12=x1-x2;  y21=y2-y1;
            x23=x2-x3;  y32=y3-y2;
            x31=x3-x1;  y13=y1-y3;            
            
            Bw=zeros(3,9);
            Bw(1,1)=1; Bw(2,4)=1; Bw(3,7)=1;
            
            Baw=zeros(6,9);
            Baw(1,[2 3 5 6])=1/2*[y21 x12 -y21 -x12];
            Baw(2,[5 6 8 9])=1/2*[y32 x23 -y32 -x23];
            Baw(3,[2 3 8 9])=1/2*[-y13 -x31 y13 x31];
            Baw(4,1:6)    =1/2*[-2 -y21 -x12 2 -y21 -x12];
            Baw(5,4:9)    =1/2*[-2 -y32 -x23 2 -y32 -x23];
            Baw(6,[1:3 7:9])=1/2*[2 -y13 -x31 -2 -y13 -x31];
            
           
            A = self.area;
            
            intNN=[1/6 1/12 1/12;
                0   1/6   1/12;
                0     0   1/6];
            
            intNN=rho*h*A*(intNN+tril(intNN',-1));
            
            intNNa=rho*h*A*[1/30 1/60 1/30 -1/180      0  -1/180;
                1/30 1/30 1/60  1/180 -1/180       0;
                1/60 1/30 1/30      0  1/180  -1/180];
            
            intNaN=intNNa';
            
            intNaNa=h*rho*A*[1/90 1/180 1/180       0 -1/1260  1/1260;
                0  1/90 1/180  1/1260       0 -1/1260;
                0     0  1/90 -1/1260  1/1260       0;
                0     0     0   1/840 -1/2520 -1/2520;
                0     0     0       0   1/840 -1/2520;
                0     0     0       0       0   1/840];
            
            intNaNa=(intNaNa+tril(intNaNa',-1));
            
            Mb=Bw'*intNN*Bw+Bw'*intNNa*Baw+Baw'*intNaN*Bw+Baw'*intNaNa*Baw;
            
            mass = A*h*rho;
            
            s = Mb(1,1) + Mb(4,4) + Mb(7,7);
            
            lumped = diag(Mb).*(mass./s);
            
            Mblu = diag(lumped);
        end
        
        function F2 = f2_local(self,x)
            LQ = self.LocalQuantities;            
            % N
            BNL=[x.' * LQ.Kxx;
                x.' * LQ.Kyy;
                x.' * LQ.Kxy];
            
            N=LQ.Am * (LQ.BL + BNL/2) * x;
            % Internal forces
            F2 =  LQ.A * BNL.' * LQ.Am * LQ.BL * x + LQ.A * LQ.BL.' * LQ.Am * (BNL/2) * x;
            LQ.A * (LQ.BL + BNL).' * LQ.Am * (LQ.BL + BNL/2) * x - LQ.A*(LQ.BL.'*(LQ.Am*(LQ.BL*x)));
        end
        
        function T2 = T2_local(self)
            m = 18;
            T2 = tenzeros([m,m,m]);
            LQ = self.LocalQuantities;
            
            C = self.area * LQ.BL.' * LQ.Am;
            D = C.'; %self.area * LQ.Am * LQ.BL; 
            for i = 1:m
                K2_1 = 1/2*(C(i,1)*(LQ.Kxx) + C(i,2)*(LQ.Kyy) + C(i,3)*(LQ.Kxy)) ;             % corresponds to A*BL'*Am*BNL
                K2_2 = LQ.Kxx(:,i)*D(1,:) + LQ.Kyy(:,i)*D(2,:) + LQ.Kxy(:,i)*D(3,:);     % corresponds to BNL'*Am*BL*A
                T2(i,:,:) = K2_1 + K2_2 ;                
            end
        end
        
        function F3 = f3_local(self,x)
            LQ = self.LocalQuantities;            
            % N
            BNL=[x.' * LQ.Kxx;
                x.' * LQ.Kyy;
                x.' * LQ.Kxy];

            F3 =  LQ.A * BNL.' * LQ.Am * (BNL/2) * x  ;            
        end
        
        function T3 = T3_local(self)
            m = 18;
            T3 = tenzeros([m,m,m,m]);
            Am = self.LocalQuantities.Am;
            Kxx = self.LocalQuantities.Kxx;
            Kxy = self.LocalQuantities.Kxy;
            Kyy = self.LocalQuantities.Kyy; 
            A = self.area;
            for i = 1:m
                for j = 1:m
                    K3_1 = 1/2*A* (( Kxx(i,j)*Am(1,1)+ Kyy(i,j)*Am(2,1)+Kxy(i,j)*Am(3,1) )*Kxx + ...   % corresponds to A*BNL'*Am*BNL
                        ( Kxx(i,j)*Am(1,2)+ Kyy(i,j)*Am(2,2)+Kxy(i,j)*Am(3,2) )*Kyy + ...
                        ( Kxx(i,j)*Am(1,3)+ Kyy(i,j)*Am(2,3)+Kxy(i,j)*Am(3,3) )*Kxy) ;                    
                    T3(i,j,:,:) = K3_1;
                end
            end
        end
        
    end

end
