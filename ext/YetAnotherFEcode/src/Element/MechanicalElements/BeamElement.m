classdef BeamElement < Element

    properties              % Basis properties derived from Element Class      
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes    
        nDOFPerNode = 3     % number of DOFs per node
        nNodes = 2          % number of nodes per element
        nDim = 1            % number of dimensions in local coordinates
    end
    
    properties                  % Properties specific to BeamElement        
        area = 0                % Area of cross section of the beam
        Material = Material()   % Object of class Material
        areaMoment = 0          % Second moment of area of he beam
        transformationMatrix
        IsUpdated
    end
    
    properties (Dependent) % Special properties calculated on the fly 
        dx            
        uniformBodyForce
        massMatrixLocal
        dampingMatrixLocal
        stiffnessMatrixLocal
    end
        
    methods
        function self = BeamElement(b, h, Material)
            narginchk(3,3)
            self.area = b*h;
            self.areaMoment = (b*h^3)/12;
            self.Material = Material;            
        end
        
        function set.nodes(self,nodes)
            self.nodes = nodes;
            not_updated(self);
        end
        
        function not_updated(self)
            self.IsUpdated.transformationMatrix = false;
        end
        %% GET methods: Define how properties are computed
        
        function dx = get.dx(self)
            dx = norm(diff(self.nodes));
        end
        
        function T = get.transformationMatrix(self)
            if ~self.IsUpdated.transformationMatrix
                update_transformation_matrix(self)
            end
            T = self.transformationMatrix;
        end
        
        function update_transformation_matrix(self)
            T = speye(6,6);
            R = rotation_matrix(self.nodes);
            T(1:2,1:2) = R;
            T(4:5,4:5) = R;
            self.transformationMatrix = T;
            self.IsUpdated.transformationMatrix = true; 
        end
        
        function M = get.massMatrixLocal(self)
            l = self.dx;
            M = sparse(6,6);
            M([1 4],[1 4]) = (self.Material.DENSITY * self.area * l/6) * ...
                [2 1;
                1 2];
            M([2 3 5 6],[2 3 5 6]) = ((self.Material.DENSITY * self.area * l)/420) * ...
                [156      22*l       54      -13*l;
                22*l     4*l^2      13*l    -3*l^2;
                54       13*l       156     -22*l;
                -13*l     -3*l^2     -22*l   4*l^2];
        end
        
        function C = get.dampingMatrixLocal(self)
            % Element level linear material damping (proportional to stiffness matrix)
            l = self.dx;
            C = sparse(6,6);
            C([1 4],[1 4]) = self.Material.DAMPING_MODULUS * self.area * [  1/l, -1/l;
                -1/l,  1/l];
            C([2 3 5 6],[2 3 5 6]) = self.Material.DAMPING_MODULUS * self.areaMoment * ...
                [  12/l^3,  6/l^2, -12/l^3,  6/l^2;
                6/l^2,    4/l,  -6/l^2,    2/l;
                -12/l^3, -6/l^2,  12/l^3, -6/l^2;
                6/l^2,    2/l,  -6/l^2,    4/l];
        end
        
        function K = get.stiffnessMatrixLocal(self)
            l = self.dx;
            K = sparse(6,6);
            K([1 4],[1 4]) = self.Material.YOUNGS_MODULUS * self.area * [  1/l, -1/l;
                -1/l,  1/l];
            K([2 3 5 6],[2 3 5 6]) = self.Material.YOUNGS_MODULUS * self.areaMoment * ...
                [  12/l^3,  6/l^2, -12/l^3,  6/l^2;
                6/l^2,    4/l,  -6/l^2,    2/l;
                -12/l^3, -6/l^2,  12/l^3, -6/l^2;
                6/l^2,    2/l,  -6/l^2,    4/l];
            
        end
        
        function  f = get.uniformBodyForce(self)
            % this function computes the element stiffness matrix and
            % internal force vector in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            
            T_e = self.transformationMatrix;
            l = self.dx;
            
            % local uniform force (integration of shape functions over
            % length of the element)
            f_local = [   l/2;
                l/2;
                l^2/12;
                l/2;
                l/2;
                -l^2/12];
            
            % only transverse body force: setting axial body force set to zero'
            f_local([1, 4]) = 0;
            
            f = T_e.' * f_local;
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
        
        function [S] = strain_energy(self,x)         
            % this function computes the strain energy of the element
            % nodes : matrix containing Nodal coordinates
            % x : vector of full DOFs in global coordinates    
            
            x_e = self.extract_element_data(x);
            T_e = self.transformationMatrix;
            % Displacements in local coordinates
            q = T_e*x_e; 
            % Strain energy           
            S = self.strain_energy_local(q);
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
        
        function f2 = F2(self,x)
            % this function computes the quadratic component of the
            % nonlinear internal force in global coordinates at the element
            % level.
            x_e = self.extract_element_data(x);
            T_e = self.transformationMatrix;            
            q = T_e*x_e;
            
            f2 = T_e.' * self.S_2(q);           
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
            m = length(index);
            T2e = self.T2_local([m,m,m]);
            
            % rotate to global coordinates
            T2 = ttm(T2e, {T_e.',T_e.',T_e.'},[1,2,3]);
        end
        
        function f3 = F3(self,x)
            % this function computes the quadratic component of the
            % nonlinear internal force in global coordinates at the element
            % level.
            T_e = self.transformationMatrix;
            x_e = self.extract_element_data(x);
            q = T_e*x_e;
            
            f3 = T_e.' * self.S_3(q);           
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
            m = length(index);
            T3e = self.T3_local([m,m,m,m]);
            
            % rotate to global coordinates
            T3 = ttm(T3e, {T_e.',T_e.',T_e.',T_e.'},[1,2,3,4]);
        end
        
        %% Local level methods
        
        function K = tangent_stiffness_matrix_local(self,x)
            l = self.dx;
            E = self.Material.YOUNGS_MODULUS;
            I = self.areaMoment;
            A = self.area;
            
            u1 = x(1);
            w1 = x(2);
            t1 = x(3);
            u2 = x(4);
            w2 = x(5);
            t2 = x(6);
            
            K = [                                                                 (A*E)/l,                                                                                                                                                                                                                                        -((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                                                   - (A*E*(4*t1 - t2))/30 - (A*E*(3*w1 - 3*w2))/(30*l),                                                                        -(A*E)/l,                                                                                                                                                                                                                                         ((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                                                     (A*E*(t1 - 4*t2))/30 - (A*E*(3*w1 - 3*w2))/(30*l);
                (6*A*E*w2)/(5*l^2) - (A*E*t2)/(10*l) - (6*A*E*w1)/(5*l^2) - (A*E*t1)/(10*l), (12*E*I)/l^3 - (6*A*E*u1)/(5*l^2) + (6*A*E*u2)/(5*l^2) + (9*A*E*t1^2)/(70*l) + (9*A*E*t2^2)/(70*l) + (108*A*E*w1^2)/(35*l^3) + (108*A*E*w2^2)/(35*l^3) + (27*A*E*t1*w1)/(35*l^2) - (27*A*E*t1*w2)/(35*l^2) + (27*A*E*t2*w1)/(35*l^2) - (27*A*E*t2*w2)/(35*l^2) - (216*A*E*w1*w2)/(35*l^3),                           (6*E*I)/l^2 - (3*A*E*t1^2)/280 + (3*A*E*t2^2)/280 - (A*E*u1)/(10*l) + (A*E*u2)/(10*l) + (27*A*E*w1^2)/(70*l^2) + (27*A*E*w2^2)/(70*l^2) + (3*A*E*t1*t2)/140 + (9*A*E*t1*w1)/(35*l) - (9*A*E*t1*w2)/(35*l) - (27*A*E*w1*w2)/(35*l^2), (A*E*t1)/(10*l) + (A*E*t2)/(10*l) + (6*A*E*w1)/(5*l^2) - (6*A*E*w2)/(5*l^2), (6*A*E*u1)/(5*l^2) - (12*E*I)/l^3 - (6*A*E*u2)/(5*l^2) - (9*A*E*t1^2)/(70*l) - (9*A*E*t2^2)/(70*l) - (108*A*E*w1^2)/(35*l^3) - (108*A*E*w2^2)/(35*l^3) - (27*A*E*t1*w1)/(35*l^2) + (27*A*E*t1*w2)/(35*l^2) - (27*A*E*t2*w1)/(35*l^2) + (27*A*E*t2*w2)/(35*l^2) + (216*A*E*w1*w2)/(35*l^3),                           (6*E*I)/l^2 + (3*A*E*t1^2)/280 - (3*A*E*t2^2)/280 - (A*E*u1)/(10*l) + (A*E*u2)/(10*l) + (27*A*E*w1^2)/(70*l^2) + (27*A*E*w2^2)/(70*l^2) + (3*A*E*t1*t2)/140 + (9*A*E*t2*w1)/(35*l) - (9*A*E*t2*w2)/(35*l) - (27*A*E*w1*w2)/(35*l^2);
                (A*E*t2)/30 - (2*A*E*t1)/15 - (A*E*w1)/(10*l) + (A*E*w2)/(10*l),                                                       (6*E*I)/l^2 - (3*A*E*t1^2)/280 + (3*A*E*t2^2)/280 - (A*E*u1)/(10*l) + (A*E*u2)/(10*l) + (27*A*E*w1^2)/(70*l^2) + (27*A*E*w2^2)/(70*l^2) + (3*A*E*t1*t2)/140 + (9*A*E*t1*w1)/(35*l) - (9*A*E*t1*w2)/(35*l) - (27*A*E*w1*w2)/(35*l^2), (3*A*E*l*t1^2)/35 - (3*A*E*l*t1*t2)/140 - (3*A*E*t1*w1)/140 + (3*A*E*t1*w2)/140 + (A*E*l*t2^2)/140 + (3*A*E*t2*w1)/140 - (3*A*E*t2*w2)/140 + (9*A*E*w1^2)/(70*l) - (9*A*E*w1*w2)/(35*l) + (9*A*E*w2^2)/(70*l) - (2*A*E*u1)/15 + (2*A*E*u2)/15 + (4*E*I)/l,             (2*A*E*t1)/15 - (A*E*t2)/30 + (A*E*w1)/(10*l) - (A*E*w2)/(10*l),                                                       (3*A*E*t1^2)/280 - (6*E*I)/l^2 - (3*A*E*t2^2)/280 + (A*E*u1)/(10*l) - (A*E*u2)/(10*l) - (27*A*E*w1^2)/(70*l^2) - (27*A*E*w2^2)/(70*l^2) - (3*A*E*t1*t2)/140 - (9*A*E*t1*w1)/(35*l) + (9*A*E*t1*w2)/(35*l) + (27*A*E*w1*w2)/(35*l^2),                                                                        (A*E*u1)/30 - (A*E*u2)/30 + (2*E*I)/l - (3*A*E*l*t1^2)/280 - (3*A*E*l*t2^2)/280 + (3*A*E*t1*w1)/140 - (3*A*E*t1*w2)/140 + (3*A*E*t2*w1)/140 - (3*A*E*t2*w2)/140 + (A*E*l*t1*t2)/70;
                -(A*E)/l,                                                                                                                                                                                                                                         ((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                                                     (A*E*(4*t1 - t2))/30 + (A*E*(3*w1 - 3*w2))/(30*l),                                                                         (A*E)/l,                                                                                                                                                                                                                                        -((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                                                     (A*E*(3*w1 - 3*w2))/(30*l) - (A*E*(t1 - 4*t2))/30;
                (A*E*t1)/(10*l) + (A*E*t2)/(10*l) + (6*A*E*w1)/(5*l^2) - (6*A*E*w2)/(5*l^2), (6*A*E*u1)/(5*l^2) - (12*E*I)/l^3 - (6*A*E*u2)/(5*l^2) - (9*A*E*t1^2)/(70*l) - (9*A*E*t2^2)/(70*l) - (108*A*E*w1^2)/(35*l^3) - (108*A*E*w2^2)/(35*l^3) - (27*A*E*t1*w1)/(35*l^2) + (27*A*E*t1*w2)/(35*l^2) - (27*A*E*t2*w1)/(35*l^2) + (27*A*E*t2*w2)/(35*l^2) + (216*A*E*w1*w2)/(35*l^3),                           (3*A*E*t1^2)/280 - (6*E*I)/l^2 - (3*A*E*t2^2)/280 + (A*E*u1)/(10*l) - (A*E*u2)/(10*l) - (27*A*E*w1^2)/(70*l^2) - (27*A*E*w2^2)/(70*l^2) - (3*A*E*t1*t2)/140 - (9*A*E*t1*w1)/(35*l) + (9*A*E*t1*w2)/(35*l) + (27*A*E*w1*w2)/(35*l^2), (6*A*E*w2)/(5*l^2) - (A*E*t2)/(10*l) - (6*A*E*w1)/(5*l^2) - (A*E*t1)/(10*l), (12*E*I)/l^3 - (6*A*E*u1)/(5*l^2) + (6*A*E*u2)/(5*l^2) + (9*A*E*t1^2)/(70*l) + (9*A*E*t2^2)/(70*l) + (108*A*E*w1^2)/(35*l^3) + (108*A*E*w2^2)/(35*l^3) + (27*A*E*t1*w1)/(35*l^2) - (27*A*E*t1*w2)/(35*l^2) + (27*A*E*t2*w1)/(35*l^2) - (27*A*E*t2*w2)/(35*l^2) - (216*A*E*w1*w2)/(35*l^3),                           (3*A*E*t2^2)/280 - (3*A*E*t1^2)/280 - (6*E*I)/l^2 + (A*E*u1)/(10*l) - (A*E*u2)/(10*l) - (27*A*E*w1^2)/(70*l^2) - (27*A*E*w2^2)/(70*l^2) - (3*A*E*t1*t2)/140 - (9*A*E*t2*w1)/(35*l) + (9*A*E*t2*w2)/(35*l) + (27*A*E*w1*w2)/(35*l^2);
                (A*E*t1)/30 - (2*A*E*t2)/15 - (A*E*w1)/(10*l) + (A*E*w2)/(10*l),                                                       (6*E*I)/l^2 + (3*A*E*t1^2)/280 - (3*A*E*t2^2)/280 - (A*E*u1)/(10*l) + (A*E*u2)/(10*l) + (27*A*E*w1^2)/(70*l^2) + (27*A*E*w2^2)/(70*l^2) + (3*A*E*t1*t2)/140 + (9*A*E*t2*w1)/(35*l) - (9*A*E*t2*w2)/(35*l) - (27*A*E*w1*w2)/(35*l^2),                                                                        (A*E*u1)/30 - (A*E*u2)/30 + (2*E*I)/l - (3*A*E*l*t1^2)/280 - (3*A*E*l*t2^2)/280 + (3*A*E*t1*w1)/140 - (3*A*E*t1*w2)/140 + (3*A*E*t2*w1)/140 - (3*A*E*t2*w2)/140 + (A*E*l*t1*t2)/70,             (2*A*E*t2)/15 - (A*E*t1)/30 + (A*E*w1)/(10*l) - (A*E*w2)/(10*l),                                                       (3*A*E*t2^2)/280 - (3*A*E*t1^2)/280 - (6*E*I)/l^2 + (A*E*u1)/(10*l) - (A*E*u2)/(10*l) - (27*A*E*w1^2)/(70*l^2) - (27*A*E*w2^2)/(70*l^2) - (3*A*E*t1*t2)/140 - (9*A*E*t2*w1)/(35*l) + (9*A*E*t2*w2)/(35*l) + (27*A*E*w1*w2)/(35*l^2), (A*E*l*t1^2)/140 - (3*A*E*l*t1*t2)/140 + (3*A*E*t1*w1)/140 - (3*A*E*t1*w2)/140 + (3*A*E*l*t2^2)/35 - (3*A*E*t2*w1)/140 + (3*A*E*t2*w2)/140 + (9*A*E*w1^2)/(70*l) - (9*A*E*w1*w2)/(35*l) + (9*A*E*w2^2)/(70*l) - (2*A*E*u1)/15 + (2*A*E*u2)/15 + (4*E*I)/l]; %#ok<*PROPLC>
            
        end
        
        function F = internal_force_local(self,x)
            l = self.dx;
            E = self.Material.YOUNGS_MODULUS;
            I = self.areaMoment;
            A = self.area;
            
            u1 = x(1);
            w1 = x(2);
            t1 = x(3);
            u2 = x(4);
            w2 = x(5);
            t2 = x(6);
            
            F = [                                                                                                                                                                                                                                                                                                                                                                                       - ((A*E*(18*w1^2 - 36*w1*w2 + 18*w2^2))/30 - (A*E*l*(30*u1 - 30*u2 - 3*t1*w1 + 3*t1*w2 - 3*t2*w1 + 3*t2*w2))/30)/l^2 - (A*E*(2*t1^2 - t1*t2 + 2*t2^2))/30;
                (E*(3360*I*w1 - 3360*I*w2 + 288*A*w1^3 - 288*A*w2^3 - A*l^3*t1^3 - A*l^3*t2^3 + 1680*I*l*t1 + 1680*I*l*t2 + 864*A*w1*w2^2 - 864*A*w1^2*w2 - 28*A*l^2*t1*u1 + 28*A*l^2*t1*u2 - 28*A*l^2*t2*u1 + 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 + 3*A*l^3*t1*t2^2 + 3*A*l^3*t1^2*t2 + 36*A*l^2*t1^2*w1 - 36*A*l^2*t1^2*w2 + 36*A*l^2*t2^2*w1 - 36*A*l^2*t2^2*w2 - 336*A*l*u1*w1 + 336*A*l*u1*w2 + 336*A*l*u2*w1 - 336*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 216*A*l*t2*w1*w2))/(280*l^3);
                (E*(5040*I*w1 - 5040*I*w2 + 108*A*w1^3 - 108*A*w2^3 + 24*A*l^3*t1^3 - 3*A*l^3*t2^3 + 3360*I*l*t1 + 1680*I*l*t2 + 324*A*w1*w2^2 - 324*A*w1^2*w2 - 112*A*l^2*t1*u1 + 112*A*l^2*t1*u2 + 28*A*l^2*t2*u1 - 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 6*A*l^3*t1*t2^2 - 9*A*l^3*t1^2*t2 - 9*A*l^2*t1^2*w1 + 9*A*l^2*t1^2*w2 + 9*A*l^2*t2^2*w1 - 9*A*l^2*t2^2*w2 - 84*A*l*u1*w1 + 84*A*l*u1*w2 + 84*A*l*u2*w1 - 84*A*l*u2*w2 - 216*A*l*t1*w1*w2 + 18*A*l^2*t1*t2*w1 - 18*A*l^2*t1*t2*w2))/(840*l^2);
                ((A*E*(18*w1^2 - 36*w1*w2 + 18*w2^2))/30 - (A*E*l*(30*u1 - 30*u2 - 3*t1*w1 + 3*t1*w2 - 3*t2*w1 + 3*t2*w2))/30)/l^2 + (A*E*(2*t1^2 - t1*t2 + 2*t2^2))/30;
                -(E*(3360*I*w1 - 3360*I*w2 + 288*A*w1^3 - 288*A*w2^3 - A*l^3*t1^3 - A*l^3*t2^3 + 1680*I*l*t1 + 1680*I*l*t2 + 864*A*w1*w2^2 - 864*A*w1^2*w2 - 28*A*l^2*t1*u1 + 28*A*l^2*t1*u2 - 28*A*l^2*t2*u1 + 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 + 3*A*l^3*t1*t2^2 + 3*A*l^3*t1^2*t2 + 36*A*l^2*t1^2*w1 - 36*A*l^2*t1^2*w2 + 36*A*l^2*t2^2*w1 - 36*A*l^2*t2^2*w2 - 336*A*l*u1*w1 + 336*A*l*u1*w2 + 336*A*l*u2*w1 - 336*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 216*A*l*t2*w1*w2))/(280*l^3);
                (E*(5040*I*w1 - 5040*I*w2 + 108*A*w1^3 - 108*A*w2^3 - 3*A*l^3*t1^3 + 24*A*l^3*t2^3 + 1680*I*l*t1 + 3360*I*l*t2 + 324*A*w1*w2^2 - 324*A*w1^2*w2 + 28*A*l^2*t1*u1 - 28*A*l^2*t1*u2 - 112*A*l^2*t2*u1 + 112*A*l^2*t2*u2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 - 9*A*l^3*t1*t2^2 + 6*A*l^3*t1^2*t2 + 9*A*l^2*t1^2*w1 - 9*A*l^2*t1^2*w2 - 9*A*l^2*t2^2*w1 + 9*A*l^2*t2^2*w2 - 84*A*l*u1*w1 + 84*A*l*u1*w2 + 84*A*l*u2*w1 - 84*A*l*u2*w2 - 216*A*l*t2*w1*w2 + 18*A*l^2*t1*t2*w1 - 18*A*l^2*t1*t2*w2))/(840*l^2)];
            
        end
        
        function [K,F] = tangent_stiffness_and_force_local(self,x)
            l = self.dx;
            E = self.Material.YOUNGS_MODULUS;
            I = self.areaMoment;
            A = self.area;
            
            u1 = x(1);
            w1 = x(2);
            t1 = x(3);
            u2 = x(4);
            w2 = x(5);
            t2 = x(6);
            
            K = [                                                                 (A*E)/l,                                                                                                                                                                                                                                        -((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                                                   - (A*E*(4*t1 - t2))/30 - (A*E*(3*w1 - 3*w2))/(30*l),                                                                        -(A*E)/l,                                                                                                                                                                                                                                         ((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                                                     (A*E*(t1 - 4*t2))/30 - (A*E*(3*w1 - 3*w2))/(30*l);
                (6*A*E*w2)/(5*l^2) - (A*E*t2)/(10*l) - (6*A*E*w1)/(5*l^2) - (A*E*t1)/(10*l), (12*E*I)/l^3 - (6*A*E*u1)/(5*l^2) + (6*A*E*u2)/(5*l^2) + (9*A*E*t1^2)/(70*l) + (9*A*E*t2^2)/(70*l) + (108*A*E*w1^2)/(35*l^3) + (108*A*E*w2^2)/(35*l^3) + (27*A*E*t1*w1)/(35*l^2) - (27*A*E*t1*w2)/(35*l^2) + (27*A*E*t2*w1)/(35*l^2) - (27*A*E*t2*w2)/(35*l^2) - (216*A*E*w1*w2)/(35*l^3),                           (6*E*I)/l^2 - (3*A*E*t1^2)/280 + (3*A*E*t2^2)/280 - (A*E*u1)/(10*l) + (A*E*u2)/(10*l) + (27*A*E*w1^2)/(70*l^2) + (27*A*E*w2^2)/(70*l^2) + (3*A*E*t1*t2)/140 + (9*A*E*t1*w1)/(35*l) - (9*A*E*t1*w2)/(35*l) - (27*A*E*w1*w2)/(35*l^2), (A*E*t1)/(10*l) + (A*E*t2)/(10*l) + (6*A*E*w1)/(5*l^2) - (6*A*E*w2)/(5*l^2), (6*A*E*u1)/(5*l^2) - (12*E*I)/l^3 - (6*A*E*u2)/(5*l^2) - (9*A*E*t1^2)/(70*l) - (9*A*E*t2^2)/(70*l) - (108*A*E*w1^2)/(35*l^3) - (108*A*E*w2^2)/(35*l^3) - (27*A*E*t1*w1)/(35*l^2) + (27*A*E*t1*w2)/(35*l^2) - (27*A*E*t2*w1)/(35*l^2) + (27*A*E*t2*w2)/(35*l^2) + (216*A*E*w1*w2)/(35*l^3),                           (6*E*I)/l^2 + (3*A*E*t1^2)/280 - (3*A*E*t2^2)/280 - (A*E*u1)/(10*l) + (A*E*u2)/(10*l) + (27*A*E*w1^2)/(70*l^2) + (27*A*E*w2^2)/(70*l^2) + (3*A*E*t1*t2)/140 + (9*A*E*t2*w1)/(35*l) - (9*A*E*t2*w2)/(35*l) - (27*A*E*w1*w2)/(35*l^2);
                (A*E*t2)/30 - (2*A*E*t1)/15 - (A*E*w1)/(10*l) + (A*E*w2)/(10*l),                                                       (6*E*I)/l^2 - (3*A*E*t1^2)/280 + (3*A*E*t2^2)/280 - (A*E*u1)/(10*l) + (A*E*u2)/(10*l) + (27*A*E*w1^2)/(70*l^2) + (27*A*E*w2^2)/(70*l^2) + (3*A*E*t1*t2)/140 + (9*A*E*t1*w1)/(35*l) - (9*A*E*t1*w2)/(35*l) - (27*A*E*w1*w2)/(35*l^2), (3*A*E*l*t1^2)/35 - (3*A*E*l*t1*t2)/140 - (3*A*E*t1*w1)/140 + (3*A*E*t1*w2)/140 + (A*E*l*t2^2)/140 + (3*A*E*t2*w1)/140 - (3*A*E*t2*w2)/140 + (9*A*E*w1^2)/(70*l) - (9*A*E*w1*w2)/(35*l) + (9*A*E*w2^2)/(70*l) - (2*A*E*u1)/15 + (2*A*E*u2)/15 + (4*E*I)/l,             (2*A*E*t1)/15 - (A*E*t2)/30 + (A*E*w1)/(10*l) - (A*E*w2)/(10*l),                                                       (3*A*E*t1^2)/280 - (6*E*I)/l^2 - (3*A*E*t2^2)/280 + (A*E*u1)/(10*l) - (A*E*u2)/(10*l) - (27*A*E*w1^2)/(70*l^2) - (27*A*E*w2^2)/(70*l^2) - (3*A*E*t1*t2)/140 - (9*A*E*t1*w1)/(35*l) + (9*A*E*t1*w2)/(35*l) + (27*A*E*w1*w2)/(35*l^2),                                                                        (A*E*u1)/30 - (A*E*u2)/30 + (2*E*I)/l - (3*A*E*l*t1^2)/280 - (3*A*E*l*t2^2)/280 + (3*A*E*t1*w1)/140 - (3*A*E*t1*w2)/140 + (3*A*E*t2*w1)/140 - (3*A*E*t2*w2)/140 + (A*E*l*t1*t2)/70;
                -(A*E)/l,                                                                                                                                                                                                                                         ((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                                                     (A*E*(4*t1 - t2))/30 + (A*E*(3*w1 - 3*w2))/(30*l),                                                                         (A*E)/l,                                                                                                                                                                                                                                        -((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                                                     (A*E*(3*w1 - 3*w2))/(30*l) - (A*E*(t1 - 4*t2))/30;
                (A*E*t1)/(10*l) + (A*E*t2)/(10*l) + (6*A*E*w1)/(5*l^2) - (6*A*E*w2)/(5*l^2), (6*A*E*u1)/(5*l^2) - (12*E*I)/l^3 - (6*A*E*u2)/(5*l^2) - (9*A*E*t1^2)/(70*l) - (9*A*E*t2^2)/(70*l) - (108*A*E*w1^2)/(35*l^3) - (108*A*E*w2^2)/(35*l^3) - (27*A*E*t1*w1)/(35*l^2) + (27*A*E*t1*w2)/(35*l^2) - (27*A*E*t2*w1)/(35*l^2) + (27*A*E*t2*w2)/(35*l^2) + (216*A*E*w1*w2)/(35*l^3),                           (3*A*E*t1^2)/280 - (6*E*I)/l^2 - (3*A*E*t2^2)/280 + (A*E*u1)/(10*l) - (A*E*u2)/(10*l) - (27*A*E*w1^2)/(70*l^2) - (27*A*E*w2^2)/(70*l^2) - (3*A*E*t1*t2)/140 - (9*A*E*t1*w1)/(35*l) + (9*A*E*t1*w2)/(35*l) + (27*A*E*w1*w2)/(35*l^2), (6*A*E*w2)/(5*l^2) - (A*E*t2)/(10*l) - (6*A*E*w1)/(5*l^2) - (A*E*t1)/(10*l), (12*E*I)/l^3 - (6*A*E*u1)/(5*l^2) + (6*A*E*u2)/(5*l^2) + (9*A*E*t1^2)/(70*l) + (9*A*E*t2^2)/(70*l) + (108*A*E*w1^2)/(35*l^3) + (108*A*E*w2^2)/(35*l^3) + (27*A*E*t1*w1)/(35*l^2) - (27*A*E*t1*w2)/(35*l^2) + (27*A*E*t2*w1)/(35*l^2) - (27*A*E*t2*w2)/(35*l^2) - (216*A*E*w1*w2)/(35*l^3),                           (3*A*E*t2^2)/280 - (3*A*E*t1^2)/280 - (6*E*I)/l^2 + (A*E*u1)/(10*l) - (A*E*u2)/(10*l) - (27*A*E*w1^2)/(70*l^2) - (27*A*E*w2^2)/(70*l^2) - (3*A*E*t1*t2)/140 - (9*A*E*t2*w1)/(35*l) + (9*A*E*t2*w2)/(35*l) + (27*A*E*w1*w2)/(35*l^2);
                (A*E*t1)/30 - (2*A*E*t2)/15 - (A*E*w1)/(10*l) + (A*E*w2)/(10*l),                                                       (6*E*I)/l^2 + (3*A*E*t1^2)/280 - (3*A*E*t2^2)/280 - (A*E*u1)/(10*l) + (A*E*u2)/(10*l) + (27*A*E*w1^2)/(70*l^2) + (27*A*E*w2^2)/(70*l^2) + (3*A*E*t1*t2)/140 + (9*A*E*t2*w1)/(35*l) - (9*A*E*t2*w2)/(35*l) - (27*A*E*w1*w2)/(35*l^2),                                                                        (A*E*u1)/30 - (A*E*u2)/30 + (2*E*I)/l - (3*A*E*l*t1^2)/280 - (3*A*E*l*t2^2)/280 + (3*A*E*t1*w1)/140 - (3*A*E*t1*w2)/140 + (3*A*E*t2*w1)/140 - (3*A*E*t2*w2)/140 + (A*E*l*t1*t2)/70,             (2*A*E*t2)/15 - (A*E*t1)/30 + (A*E*w1)/(10*l) - (A*E*w2)/(10*l),                                                       (3*A*E*t2^2)/280 - (3*A*E*t1^2)/280 - (6*E*I)/l^2 + (A*E*u1)/(10*l) - (A*E*u2)/(10*l) - (27*A*E*w1^2)/(70*l^2) - (27*A*E*w2^2)/(70*l^2) - (3*A*E*t1*t2)/140 - (9*A*E*t2*w1)/(35*l) + (9*A*E*t2*w2)/(35*l) + (27*A*E*w1*w2)/(35*l^2), (A*E*l*t1^2)/140 - (3*A*E*l*t1*t2)/140 + (3*A*E*t1*w1)/140 - (3*A*E*t1*w2)/140 + (3*A*E*l*t2^2)/35 - (3*A*E*t2*w1)/140 + (3*A*E*t2*w2)/140 + (9*A*E*w1^2)/(70*l) - (9*A*E*w1*w2)/(35*l) + (9*A*E*w2^2)/(70*l) - (2*A*E*u1)/15 + (2*A*E*u2)/15 + (4*E*I)/l]; %#ok<*PROPLC>
            
            F = [                                                                                                                                                                                                                                                                                                                                                                                       - ((A*E*(18*w1^2 - 36*w1*w2 + 18*w2^2))/30 - (A*E*l*(30*u1 - 30*u2 - 3*t1*w1 + 3*t1*w2 - 3*t2*w1 + 3*t2*w2))/30)/l^2 - (A*E*(2*t1^2 - t1*t2 + 2*t2^2))/30;
                (E*(3360*I*w1 - 3360*I*w2 + 288*A*w1^3 - 288*A*w2^3 - A*l^3*t1^3 - A*l^3*t2^3 + 1680*I*l*t1 + 1680*I*l*t2 + 864*A*w1*w2^2 - 864*A*w1^2*w2 - 28*A*l^2*t1*u1 + 28*A*l^2*t1*u2 - 28*A*l^2*t2*u1 + 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 + 3*A*l^3*t1*t2^2 + 3*A*l^3*t1^2*t2 + 36*A*l^2*t1^2*w1 - 36*A*l^2*t1^2*w2 + 36*A*l^2*t2^2*w1 - 36*A*l^2*t2^2*w2 - 336*A*l*u1*w1 + 336*A*l*u1*w2 + 336*A*l*u2*w1 - 336*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 216*A*l*t2*w1*w2))/(280*l^3);
                (E*(5040*I*w1 - 5040*I*w2 + 108*A*w1^3 - 108*A*w2^3 + 24*A*l^3*t1^3 - 3*A*l^3*t2^3 + 3360*I*l*t1 + 1680*I*l*t2 + 324*A*w1*w2^2 - 324*A*w1^2*w2 - 112*A*l^2*t1*u1 + 112*A*l^2*t1*u2 + 28*A*l^2*t2*u1 - 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 6*A*l^3*t1*t2^2 - 9*A*l^3*t1^2*t2 - 9*A*l^2*t1^2*w1 + 9*A*l^2*t1^2*w2 + 9*A*l^2*t2^2*w1 - 9*A*l^2*t2^2*w2 - 84*A*l*u1*w1 + 84*A*l*u1*w2 + 84*A*l*u2*w1 - 84*A*l*u2*w2 - 216*A*l*t1*w1*w2 + 18*A*l^2*t1*t2*w1 - 18*A*l^2*t1*t2*w2))/(840*l^2);
                ((A*E*(18*w1^2 - 36*w1*w2 + 18*w2^2))/30 - (A*E*l*(30*u1 - 30*u2 - 3*t1*w1 + 3*t1*w2 - 3*t2*w1 + 3*t2*w2))/30)/l^2 + (A*E*(2*t1^2 - t1*t2 + 2*t2^2))/30;
                -(E*(3360*I*w1 - 3360*I*w2 + 288*A*w1^3 - 288*A*w2^3 - A*l^3*t1^3 - A*l^3*t2^3 + 1680*I*l*t1 + 1680*I*l*t2 + 864*A*w1*w2^2 - 864*A*w1^2*w2 - 28*A*l^2*t1*u1 + 28*A*l^2*t1*u2 - 28*A*l^2*t2*u1 + 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 + 3*A*l^3*t1*t2^2 + 3*A*l^3*t1^2*t2 + 36*A*l^2*t1^2*w1 - 36*A*l^2*t1^2*w2 + 36*A*l^2*t2^2*w1 - 36*A*l^2*t2^2*w2 - 336*A*l*u1*w1 + 336*A*l*u1*w2 + 336*A*l*u2*w1 - 336*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 216*A*l*t2*w1*w2))/(280*l^3);
                (E*(5040*I*w1 - 5040*I*w2 + 108*A*w1^3 - 108*A*w2^3 - 3*A*l^3*t1^3 + 24*A*l^3*t2^3 + 1680*I*l*t1 + 3360*I*l*t2 + 324*A*w1*w2^2 - 324*A*w1^2*w2 + 28*A*l^2*t1*u1 - 28*A*l^2*t1*u2 - 112*A*l^2*t2*u1 + 112*A*l^2*t2*u2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 - 9*A*l^3*t1*t2^2 + 6*A*l^3*t1^2*t2 + 9*A*l^2*t1^2*w1 - 9*A*l^2*t1^2*w2 - 9*A*l^2*t2^2*w1 + 9*A*l^2*t2^2*w2 - 84*A*l*u1*w1 + 84*A*l*u1*w2 + 84*A*l*u2*w1 - 84*A*l*u2*w2 - 216*A*l*t2*w1*w2 + 18*A*l^2*t1*t2*w1 - 18*A*l^2*t1*t2*w2))/(840*l^2)];
            
        end
        
        function [S] = strain_energy_local(self,x)
            dx = self.dx;
            E = self.Material.YOUNGS_MODULUS;
            I = self.areaMoment;
            A = self.area;
            
            u1 = x(1);
            w1 = x(2);
            t1 = x(3);
            u2 = x(4);
            w2 = x(5);
            t2 = x(6);
            
            
            S = 0.5*dx*((A*E*((32*u2)/dx - (32*u1)/dx + (6*w1 - 6*w2 + dx*t1 + dx*t2)^2/dx^2)^2)/1024 + (E*I*(t1 - t2)^2)/dx^2) + (dx*(A*E*(2*((t1 - t2)^2/8 - (3*(2*w1 - 2*w2 + dx*t1 + dx*t2)*(6*w1 - 6*w2 + dx*t1 + dx*t2))/(16*dx^2))*(u2/dx - u1/dx + (6*w1 - 6*w2 + dx*t1 + dx*t2)^2/(32*dx^2)) + ((t1 - t2)^2*(6*w1 - 6*w2 + dx*t1 + dx*t2)^2)/(64*dx^2)) + (9*E*I*(2*w1 - 2*w2 + dx*t1 + dx*t2)^2)/dx^4))/3 + (9*A*E*(2*w1 - 2*w2 + dx*t1 + dx*t2)^4)/(1024*dx^3) + (A*E*dx*((2*(t1 - t2)^2 - (3*(2*w1 - 2*w2 + dx*t1 + dx*t2)*(6*w1 - 6*w2 + dx*t1 + dx*t2))/dx^2)^2/256 + (9*(u2/dx - u1/dx + (6*w1 - 6*w2 + dx*t1 + dx*t2)^2/(32*dx^2))*(2*w1 - 2*w2 + dx*t1 + dx*t2)^2)/(16*dx^2) - (3*(t1 - t2)^2*(2*w1 - 2*w2 + dx*t1 + dx*t2)*(6*w1 - 6*w2 + dx*t1 + dx*t2))/(32*dx^2)))/5 + (27*A*E*(2*w1 - 2*w2 + dx*t1 + dx*t2)^2*(dx^2*t1^2 - 6*dx^2*t1*t2 + dx^2*t2^2 - 8*dx*t1*w1 + 8*dx*t1*w2 - 8*dx*t2*w1 + 8*dx*t2*w2 - 12*w1^2 + 24*w1*w2 - 12*w2^2))/(1792*dx^3);
            
        end
        
        function [K] = stiffness_derivative_local(self,x,v)
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            dx = self.dx;
            
            w1 = x(2);
            t1 = x(3);
            w2 = x(5);
            t2 = x(6);
            
            vu1 = v(1);
            vw1 = v(2);
            vt1 = v(3);
            vu2 = v(4);
            vw2 = v(5);
            vt2 = v(6);
            
            K =  [                                                    0,                                                                                                                                                                           -((A*E*(36*vw1 - 36*vw2))/30 + (A*E*dx*(3*vt1 + 3*vt2))/30)/dx^2,                                                                                                                                                                                                                     - (A*E*(4*vt1 - vt2))/30 - (A*E*(3*vw1 - 3*vw2))/(30*dx),                                                    0,                                                                                                                                                                            ((A*E*(36*vw1 - 36*vw2))/30 + (A*E*dx*(3*vt1 + 3*vt2))/30)/dx^2,                                                                                                                                                                                                                       (A*E*(vt1 - 4*vt2))/30 - (A*E*(3*vw1 - 3*vw2))/(30*dx);
                -(A*E*(12*vw1 - 12*vw2 + dx*vt1 + dx*vt2))/(10*dx^2),  (3*A*E*(14*dx*vu2 - 14*dx*vu1 + 72*vw1*w1 - 72*vw1*w2 - 72*vw2*w1 + 72*vw2*w2 + 3*dx^2*t1*vt1 + 3*dx^2*t2*vt2 + 9*dx*t1*vw1 - 9*dx*t1*vw2 + 9*dx*t2*vw1 - 9*dx*t2*vw2 + 9*dx*vt1*w1 - 9*dx*vt1*w2 + 9*dx*vt2*w1 - 9*dx*vt2*w2))/(35*dx^3),                                                     (A*E*(14*dx*vu2 - 14*dx*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 - 3*dx^2*t1*vt1 + 3*dx^2*t1*vt2 + 3*dx^2*t2*vt1 + 3*dx^2*t2*vt2 + 36*dx*t1*vw1 - 36*dx*t1*vw2 + 36*dx*vt1*w1 - 36*dx*vt1*w2))/(140*dx^2),  (A*E*(12*vw1 - 12*vw2 + dx*vt1 + dx*vt2))/(10*dx^2), -(3*A*E*(14*dx*vu2 - 14*dx*vu1 + 72*vw1*w1 - 72*vw1*w2 - 72*vw2*w1 + 72*vw2*w2 + 3*dx^2*t1*vt1 + 3*dx^2*t2*vt2 + 9*dx*t1*vw1 - 9*dx*t1*vw2 + 9*dx*t2*vw1 - 9*dx*t2*vw2 + 9*dx*vt1*w1 - 9*dx*vt1*w2 + 9*dx*vt2*w1 - 9*dx*vt2*w2))/(35*dx^3),                                                     (A*E*(14*dx*vu2 - 14*dx*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 + 3*dx^2*t1*vt1 + 3*dx^2*t1*vt2 + 3*dx^2*t2*vt1 - 3*dx^2*t2*vt2 + 36*dx*t2*vw1 - 36*dx*t2*vw2 + 36*dx*vt2*w1 - 36*dx*vt2*w2))/(140*dx^2);
                -(A*E*(3*vw1 - 3*vw2 + 4*dx*vt1 - dx*vt2))/(30*dx),                   (A*E*(14*dx*vu2 - 14*dx*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 - 3*dx^2*t1*vt1 + 3*dx^2*t1*vt2 + 3*dx^2*t2*vt1 + 3*dx^2*t2*vt2 + 36*dx*t1*vw1 - 36*dx*t1*vw2 + 36*dx*vt1*w1 - 36*dx*vt1*w2))/(140*dx^2), -(A*E*(56*dx*vu1 - 56*dx*vu2 - 108*vw1*w1 + 108*vw1*w2 + 108*vw2*w1 - 108*vw2*w2 - 72*dx^2*t1*vt1 + 9*dx^2*t1*vt2 + 9*dx^2*t2*vt1 - 6*dx^2*t2*vt2 + 9*dx*t1*vw1 - 9*dx*t1*vw2 - 9*dx*t2*vw1 + 9*dx*t2*vw2 + 9*dx*vt1*w1 - 9*dx*vt1*w2 - 9*dx*vt2*w1 + 9*dx*vt2*w2))/(420*dx),    (A*E*(3*vw1 - 3*vw2 + 4*dx*vt1 - dx*vt2))/(30*dx),                  -(A*E*(14*dx*vu2 - 14*dx*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 - 3*dx^2*t1*vt1 + 3*dx^2*t1*vt2 + 3*dx^2*t2*vt1 + 3*dx^2*t2*vt2 + 36*dx*t1*vw1 - 36*dx*t1*vw2 + 36*dx*vt1*w1 - 36*dx*vt1*w2))/(140*dx^2),                                                                                                  (A*E*(14*vu1 - 14*vu2 + 9*t1*vw1 - 9*t1*vw2 + 9*t2*vw1 - 9*t2*vw2 + 9*vt1*w1 - 9*vt1*w2 + 9*vt2*w1 - 9*vt2*w2 - 9*dx*t1*vt1 + 6*dx*t1*vt2 + 6*dx*t2*vt1 - 9*dx*t2*vt2))/420;
                0,                                                                                                                                                                            ((A*E*(36*vw1 - 36*vw2))/30 + (A*E*dx*(3*vt1 + 3*vt2))/30)/dx^2,                                                                                                                                                                                                                       (A*E*(4*vt1 - vt2))/30 + (A*E*(3*vw1 - 3*vw2))/(30*dx),                                                    0,                                                                                                                                                                           -((A*E*(36*vw1 - 36*vw2))/30 + (A*E*dx*(3*vt1 + 3*vt2))/30)/dx^2,                                                                                                                                                                                                                       (A*E*(3*vw1 - 3*vw2))/(30*dx) - (A*E*(vt1 - 4*vt2))/30;
                (A*E*(12*vw1 - 12*vw2 + dx*vt1 + dx*vt2))/(10*dx^2), -(3*A*E*(14*dx*vu2 - 14*dx*vu1 + 72*vw1*w1 - 72*vw1*w2 - 72*vw2*w1 + 72*vw2*w2 + 3*dx^2*t1*vt1 + 3*dx^2*t2*vt2 + 9*dx*t1*vw1 - 9*dx*t1*vw2 + 9*dx*t2*vw1 - 9*dx*t2*vw2 + 9*dx*vt1*w1 - 9*dx*vt1*w2 + 9*dx*vt2*w1 - 9*dx*vt2*w2))/(35*dx^3),                                                    -(A*E*(14*dx*vu2 - 14*dx*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 - 3*dx^2*t1*vt1 + 3*dx^2*t1*vt2 + 3*dx^2*t2*vt1 + 3*dx^2*t2*vt2 + 36*dx*t1*vw1 - 36*dx*t1*vw2 + 36*dx*vt1*w1 - 36*dx*vt1*w2))/(140*dx^2), -(A*E*(12*vw1 - 12*vw2 + dx*vt1 + dx*vt2))/(10*dx^2),  (3*A*E*(14*dx*vu2 - 14*dx*vu1 + 72*vw1*w1 - 72*vw1*w2 - 72*vw2*w1 + 72*vw2*w2 + 3*dx^2*t1*vt1 + 3*dx^2*t2*vt2 + 9*dx*t1*vw1 - 9*dx*t1*vw2 + 9*dx*t2*vw1 - 9*dx*t2*vw2 + 9*dx*vt1*w1 - 9*dx*vt1*w2 + 9*dx*vt2*w1 - 9*dx*vt2*w2))/(35*dx^3),                                                    -(A*E*(14*dx*vu2 - 14*dx*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 + 3*dx^2*t1*vt1 + 3*dx^2*t1*vt2 + 3*dx^2*t2*vt1 - 3*dx^2*t2*vt2 + 36*dx*t2*vw1 - 36*dx*t2*vw2 + 36*dx*vt2*w1 - 36*dx*vt2*w2))/(140*dx^2);
                -(A*E*(3*vw1 - 3*vw2 - dx*vt1 + 4*dx*vt2))/(30*dx),                   (A*E*(14*dx*vu2 - 14*dx*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 + 3*dx^2*t1*vt1 + 3*dx^2*t1*vt2 + 3*dx^2*t2*vt1 - 3*dx^2*t2*vt2 + 36*dx*t2*vw1 - 36*dx*t2*vw2 + 36*dx*vt2*w1 - 36*dx*vt2*w2))/(140*dx^2),                                                                                                  (A*E*(14*vu1 - 14*vu2 + 9*t1*vw1 - 9*t1*vw2 + 9*t2*vw1 - 9*t2*vw2 + 9*vt1*w1 - 9*vt1*w2 + 9*vt2*w1 - 9*vt2*w2 - 9*dx*t1*vt1 + 6*dx*t1*vt2 + 6*dx*t2*vt1 - 9*dx*t2*vt2))/420,    (A*E*(3*vw1 - 3*vw2 - dx*vt1 + 4*dx*vt2))/(30*dx),                  -(A*E*(14*dx*vu2 - 14*dx*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 + 3*dx^2*t1*vt1 + 3*dx^2*t1*vt2 + 3*dx^2*t2*vt1 - 3*dx^2*t2*vt2 + 36*dx*t2*vw1 - 36*dx*t2*vw2 + 36*dx*vt2*w1 - 36*dx*vt2*w2))/(140*dx^2), -(A*E*(56*dx*vu1 - 56*dx*vu2 - 108*vw1*w1 + 108*vw1*w2 + 108*vw2*w1 - 108*vw2*w2 - 6*dx^2*t1*vt1 + 9*dx^2*t1*vt2 + 9*dx^2*t2*vt1 - 72*dx^2*t2*vt2 - 9*dx*t1*vw1 + 9*dx*t1*vw2 + 9*dx*t2*vw1 - 9*dx*t2*vw2 - 9*dx*vt1*w1 + 9*dx*vt1*w2 + 9*dx*vt2*w1 - 9*dx*vt2*w2))/(420*dx)];
        end
        
        function S = S(self,x)
            l = self.dx;
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            
            u1 = x(1);
            w1 = x(2);
            t1 = x(3);
            u2 = x(4);
            w2 = x(5);
            t2 = x(6);
            
            S = [                                                                                                                                                                                                                                                                             -(A*E*(2*l^2*t1^2 - l^2*t1*t2 + 2*l^2*t2^2 + 3*l*t1*w1 - 3*l*t1*w2 + 3*l*t2*w1 - 3*l*t2*w2 + 18*w1^2 - 36*w1*w2 + 18*w2^2))/(30*l^2);
                (A*E*(864*w1*w2^2 - 864*w1^2*w2 + 288*w1^3 - 288*w2^3 - l^3*t1^3 - l^3*t2^3 - 28*l^2*t1*u1 + 28*l^2*t1*u2 - 28*l^2*t2*u1 + 28*l^2*t2*u2 + 108*l*t1*w1^2 + 108*l*t1*w2^2 + 108*l*t2*w1^2 + 108*l*t2*w2^2 + 3*l^3*t1*t2^2 + 3*l^3*t1^2*t2 + 36*l^2*t1^2*w1 - 36*l^2*t1^2*w2 + 36*l^2*t2^2*w1 - 36*l^2*t2^2*w2 - 336*l*u1*w1 + 336*l*u1*w2 + 336*l*u2*w1 - 336*l*u2*w2 - 216*l*t1*w1*w2 - 216*l*t2*w1*w2))/(280*l^3);
                (A*E*(324*w1*w2^2 - 324*w1^2*w2 + 108*w1^3 - 108*w2^3 + 24*l^3*t1^3 - 3*l^3*t2^3 - 112*l^2*t1*u1 + 112*l^2*t1*u2 + 28*l^2*t2*u1 - 28*l^2*t2*u2 + 108*l*t1*w1^2 + 108*l*t1*w2^2 + 6*l^3*t1*t2^2 - 9*l^3*t1^2*t2 - 9*l^2*t1^2*w1 + 9*l^2*t1^2*w2 + 9*l^2*t2^2*w1 - 9*l^2*t2^2*w2 - 84*l*u1*w1 + 84*l*u1*w2 + 84*l*u2*w1 - 84*l*u2*w2 - 216*l*t1*w1*w2 + 18*l^2*t1*t2*w1 - 18*l^2*t1*t2*w2))/(840*l^2);
                (A*E*(2*l^2*t1^2 - l^2*t1*t2 + 2*l^2*t2^2 + 3*l*t1*w1 - 3*l*t1*w2 + 3*l*t2*w1 - 3*l*t2*w2 + 18*w1^2 - 36*w1*w2 + 18*w2^2))/(30*l^2);
                -(A*E*(864*w1*w2^2 - 864*w1^2*w2 + 288*w1^3 - 288*w2^3 - l^3*t1^3 - l^3*t2^3 - 28*l^2*t1*u1 + 28*l^2*t1*u2 - 28*l^2*t2*u1 + 28*l^2*t2*u2 + 108*l*t1*w1^2 + 108*l*t1*w2^2 + 108*l*t2*w1^2 + 108*l*t2*w2^2 + 3*l^3*t1*t2^2 + 3*l^3*t1^2*t2 + 36*l^2*t1^2*w1 - 36*l^2*t1^2*w2 + 36*l^2*t2^2*w1 - 36*l^2*t2^2*w2 - 336*l*u1*w1 + 336*l*u1*w2 + 336*l*u2*w1 - 336*l*u2*w2 - 216*l*t1*w1*w2 - 216*l*t2*w1*w2))/(280*l^3);
                (A*E*(324*w1*w2^2 - 324*w1^2*w2 + 108*w1^3 - 108*w2^3 - 3*l^3*t1^3 + 24*l^3*t2^3 + 28*l^2*t1*u1 - 28*l^2*t1*u2 - 112*l^2*t2*u1 + 112*l^2*t2*u2 + 108*l*t2*w1^2 + 108*l*t2*w2^2 - 9*l^3*t1*t2^2 + 6*l^3*t1^2*t2 + 9*l^2*t1^2*w1 - 9*l^2*t1^2*w2 - 9*l^2*t2^2*w1 + 9*l^2*t2^2*w2 - 84*l*u1*w1 + 84*l*u1*w2 + 84*l*u2*w1 - 84*l*u2*w2 - 216*l*t2*w1*w2 + 18*l^2*t1*t2*w1 - 18*l^2*t1*t2*w2))/(840*l^2)];
            
            
        end
        
        function S_2 = S_2(self,x)
            l = self.dx;
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            
            u1 = x(1);
            w1 = x(2);
            t1 = x(3);
            u2 = x(4);
            w2 = x(5);
            t2 = x(6);
            S_2 =  [- (A*E*(2*t1^2 - t1*t2 + 2*t2^2))/30 - (3*A*E*(w1 - w2)^2)/(5*l^2) - (A*E*(t1 + t2)*(w1 - w2))/(10*l);
                -(A*E*(u1 - u2)*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2);
                - (A*E*(4*t1 - t2)*(u1 - u2))/30 - (A*E*(u1 - u2)*(w1 - w2))/(10*l);
                (A*E*(2*t1^2 - t1*t2 + 2*t2^2))/30 + (3*A*E*(w1 - w2)^2)/(5*l^2) + (A*E*(t1 + t2)*(w1 - w2))/(10*l);
                (A*E*(u1 - u2)*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2);
                (A*E*(t1 - 4*t2)*(u1 - u2))/30 - (A*E*(u1 - u2)*(w1 - w2))/(10*l)];
        end
        
        function T2 = T2_local(self,SIZE)
            dx = self.dx;
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            
            T2 = tenzeros(SIZE);
            
            
            
            T2(:,:,2) = [                 0, -(3*A*E)/(5*dx^2), -(A*E)/(20*dx),                 0,  (3*A*E)/(5*dx^2), -(A*E)/(20*dx);
                -(6*A*E)/(5*dx^2),                 0,              0,  (6*A*E)/(5*dx^2),                 0,              0;
                -(A*E)/(10*dx),                 0,              0,     (A*E)/(10*dx),                 0,              0;
                0,  (3*A*E)/(5*dx^2),  (A*E)/(20*dx),                 0, -(3*A*E)/(5*dx^2),  (A*E)/(20*dx);
                (6*A*E)/(5*dx^2),                 0,              0, -(6*A*E)/(5*dx^2),                 0,              0;
                -(A*E)/(10*dx),                 0,              0,     (A*E)/(10*dx),                 0,              0];
            
            
            T2(:,:,3) = [              0, -(A*E)/(20*dx), -(A*E)/15,              0,  (A*E)/(20*dx),  (A*E)/60;
                -(A*E)/(10*dx),              0,         0,  (A*E)/(10*dx),              0,         0;
                -(2*A*E)/15,              0,         0,     (2*A*E)/15,              0,         0;
                0,  (A*E)/(20*dx),  (A*E)/15,              0, -(A*E)/(20*dx), -(A*E)/60;
                (A*E)/(10*dx),              0,         0, -(A*E)/(10*dx),              0,         0;
                (A*E)/30,              0,         0,      -(A*E)/30,              0,         0];
            
            
            T2(:,:,5) = [                 0,  (3*A*E)/(5*dx^2),  (A*E)/(20*dx),                 0, -(3*A*E)/(5*dx^2),  (A*E)/(20*dx);
                (6*A*E)/(5*dx^2),                 0,              0, -(6*A*E)/(5*dx^2),                 0,              0;
                (A*E)/(10*dx),                 0,              0,    -(A*E)/(10*dx),                 0,              0;
                0, -(3*A*E)/(5*dx^2), -(A*E)/(20*dx),                 0,  (3*A*E)/(5*dx^2), -(A*E)/(20*dx);
                -(6*A*E)/(5*dx^2),                 0,              0,  (6*A*E)/(5*dx^2),                 0,              0;
                (A*E)/(10*dx),                 0,              0,    -(A*E)/(10*dx),                 0,              0];
            
            
            T2(:,:,6) = [              0, -(A*E)/(20*dx),  (A*E)/60,              0,  (A*E)/(20*dx), -(A*E)/15;
                -(A*E)/(10*dx),              0,         0,  (A*E)/(10*dx),              0,         0;
                (A*E)/30,              0,         0,      -(A*E)/30,              0,         0;
                0,  (A*E)/(20*dx), -(A*E)/60,              0, -(A*E)/(20*dx),  (A*E)/15;
                (A*E)/(10*dx),              0,         0, -(A*E)/(10*dx),              0,         0;
                -(2*A*E)/15,              0,         0,     (2*A*E)/15,              0,         0];
            
            
        end
            
        function S_3 = S_3(self,x)
            l = self.dx;
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            
            w1 = x(2);
            t1 = x(3);
            w2 = x(5);
            t2 = x(6);
            
            S_3 =  [                                                                                                                                                                                                                                                                                            0;
                (A*E*(- l^3*t1^3 + 3*l^3*t1^2*t2 + 3*l^3*t1*t2^2 - l^3*t2^3 + 36*l^2*t1^2*w1 - 36*l^2*t1^2*w2 + 36*l^2*t2^2*w1 - 36*l^2*t2^2*w2 + 108*l*t1*w1^2 - 216*l*t1*w1*w2 + 108*l*t1*w2^2 + 108*l*t2*w1^2 - 216*l*t2*w1*w2 + 108*l*t2*w2^2 + 288*w1^3 - 864*w1^2*w2 + 864*w1*w2^2 - 288*w2^3))/(280*l^3);
                (A*E*(8*l^3*t1^3 - 3*l^3*t1^2*t2 + 2*l^3*t1*t2^2 - l^3*t2^3 - 3*l^2*t1^2*w1 + 3*l^2*t1^2*w2 + 6*l^2*t1*t2*w1 - 6*l^2*t1*t2*w2 + 3*l^2*t2^2*w1 - 3*l^2*t2^2*w2 + 36*l*t1*w1^2 - 72*l*t1*w1*w2 + 36*l*t1*w2^2 + 36*w1^3 - 108*w1^2*w2 + 108*w1*w2^2 - 36*w2^3))/(280*l^2);
                0;
                -(A*E*(- l^3*t1^3 + 3*l^3*t1^2*t2 + 3*l^3*t1*t2^2 - l^3*t2^3 + 36*l^2*t1^2*w1 - 36*l^2*t1^2*w2 + 36*l^2*t2^2*w1 - 36*l^2*t2^2*w2 + 108*l*t1*w1^2 - 216*l*t1*w1*w2 + 108*l*t1*w2^2 + 108*l*t2*w1^2 - 216*l*t2*w1*w2 + 108*l*t2*w2^2 + 288*w1^3 - 864*w1^2*w2 + 864*w1*w2^2 - 288*w2^3))/(280*l^3);
                (A*E*(- l^3*t1^3 + 2*l^3*t1^2*t2 - 3*l^3*t1*t2^2 + 8*l^3*t2^3 + 3*l^2*t1^2*w1 - 3*l^2*t1^2*w2 + 6*l^2*t1*t2*w1 - 6*l^2*t1*t2*w2 - 3*l^2*t2^2*w1 + 3*l^2*t2^2*w2 + 36*l*t2*w1^2 - 72*l*t2*w1*w2 + 36*l*t2*w2^2 + 36*w1^3 - 108*w1^2*w2 + 108*w1*w2^2 - 36*w2^3))/(280*l^2)];
        end
        
        function T3 = T3_local(self,SIZE)
            dx = self.dx;
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            
            T3 = tenzeros(SIZE);            
            
            T3(:,:,2,2) = [ 0,                   0,                  0, 0,                   0,                  0;
             0,  (36*A*E)/(35*dx^3),  (9*A*E)/(70*dx^2), 0, -(36*A*E)/(35*dx^3),  (9*A*E)/(70*dx^2);
             0,   (9*A*E)/(70*dx^2),    (3*A*E)/(70*dx), 0,  -(9*A*E)/(70*dx^2),                  0;
             0,                   0,                  0, 0,                   0,                  0;
             0, -(36*A*E)/(35*dx^3), -(9*A*E)/(70*dx^2), 0,  (36*A*E)/(35*dx^3), -(9*A*E)/(70*dx^2);
             0,   (9*A*E)/(70*dx^2),                  0, 0,  -(9*A*E)/(70*dx^2),    (3*A*E)/(70*dx)];
            
            
            T3(:,:,3,2) = [ 0,                  0,                0, 0,                  0,         0;
             0,  (9*A*E)/(70*dx^2),  (3*A*E)/(70*dx), 0, -(9*A*E)/(70*dx^2),         0;
             0,    (3*A*E)/(70*dx),       -(A*E)/280, 0,   -(3*A*E)/(70*dx), (A*E)/280;
             0,                  0,                0, 0,                  0,         0;
             0, -(9*A*E)/(70*dx^2), -(3*A*E)/(70*dx), 0,  (9*A*E)/(70*dx^2),         0;
             0,                  0,        (A*E)/280, 0,                  0, (A*E)/280];
                        
            
            T3(:,:,5,2) = [ 0,                   0,                  0, 0,                   0,                  0;
             0, -(36*A*E)/(35*dx^3), -(9*A*E)/(70*dx^2), 0,  (36*A*E)/(35*dx^3), -(9*A*E)/(70*dx^2);
             0,  -(9*A*E)/(70*dx^2),   -(3*A*E)/(70*dx), 0,   (9*A*E)/(70*dx^2),                  0;
             0,                   0,                  0, 0,                   0,                  0;
             0,  (36*A*E)/(35*dx^3),  (9*A*E)/(70*dx^2), 0, -(36*A*E)/(35*dx^3),  (9*A*E)/(70*dx^2);
             0,  -(9*A*E)/(70*dx^2),                  0, 0,   (9*A*E)/(70*dx^2),   -(3*A*E)/(70*dx)];
            
            
            T3(:,:,6,2) = [ 0,                  0,         0, 0,                  0,                0;
             0,  (9*A*E)/(70*dx^2),         0, 0, -(9*A*E)/(70*dx^2),  (3*A*E)/(70*dx);
             0,                  0, (A*E)/280, 0,                  0,        (A*E)/280;
             0,                  0,         0, 0,                  0,                0;
             0, -(9*A*E)/(70*dx^2),         0, 0,  (9*A*E)/(70*dx^2), -(3*A*E)/(70*dx);
             0,    (3*A*E)/(70*dx), (A*E)/280, 0,   -(3*A*E)/(70*dx),       -(A*E)/280];
            
                        
            
            T3(:,:,2,3) =[ 0,                  0,                0, 0,                  0,         0;
             0,  (9*A*E)/(70*dx^2),  (3*A*E)/(70*dx), 0, -(9*A*E)/(70*dx^2),         0;
             0,    (3*A*E)/(70*dx),       -(A*E)/280, 0,   -(3*A*E)/(70*dx), (A*E)/280;
             0,                  0,                0, 0,                  0,         0;
             0, -(9*A*E)/(70*dx^2), -(3*A*E)/(70*dx), 0,  (9*A*E)/(70*dx^2),         0;
             0,                  0,        (A*E)/280, 0,                  0, (A*E)/280];
            
            
            T3(:,:,3,3) =[ 0,                0,             0, 0,                0,             0;
             0,  (3*A*E)/(70*dx),    -(A*E)/280, 0, -(3*A*E)/(70*dx),     (A*E)/280;
             0,       -(A*E)/280,   (A*E*dx)/35, 0,        (A*E)/280, -(A*E*dx)/280;
             0,                0,             0, 0,                0,             0;
             0, -(3*A*E)/(70*dx),     (A*E)/280, 0,  (3*A*E)/(70*dx),    -(A*E)/280;
             0,        (A*E)/280, -(A*E*dx)/280, 0,       -(A*E)/280,  (A*E*dx)/420];
            
                       
            
            T3(:,:,5,3) =[ 0,                  0,                0, 0,                  0,          0;
             0, -(9*A*E)/(70*dx^2), -(3*A*E)/(70*dx), 0,  (9*A*E)/(70*dx^2),          0;
             0,   -(3*A*E)/(70*dx),        (A*E)/280, 0,    (3*A*E)/(70*dx), -(A*E)/280;
             0,                  0,                0, 0,                  0,          0;
             0,  (9*A*E)/(70*dx^2),  (3*A*E)/(70*dx), 0, -(9*A*E)/(70*dx^2),          0;
             0,                  0,       -(A*E)/280, 0,                  0, -(A*E)/280];
                        
            T3(:,:,6,3) =[ 0,         0,             0, 0,          0,             0;
             0,         0,     (A*E)/280, 0,          0,     (A*E)/280;
             0, (A*E)/280, -(A*E*dx)/280, 0, -(A*E)/280,  (A*E*dx)/420;
             0,         0,             0, 0,          0,             0;
             0,         0,    -(A*E)/280, 0,          0,    -(A*E)/280;
             0, (A*E)/280,  (A*E*dx)/420, 0, -(A*E)/280, -(A*E*dx)/280];
                        
            
            
            T3(:,:,2,5) =[ 0,                   0,                  0, 0,                   0,                  0;
             0, -(36*A*E)/(35*dx^3), -(9*A*E)/(70*dx^2), 0,  (36*A*E)/(35*dx^3), -(9*A*E)/(70*dx^2);
             0,  -(9*A*E)/(70*dx^2),   -(3*A*E)/(70*dx), 0,   (9*A*E)/(70*dx^2),                  0;
             0,                   0,                  0, 0,                   0,                  0;
             0,  (36*A*E)/(35*dx^3),  (9*A*E)/(70*dx^2), 0, -(36*A*E)/(35*dx^3),  (9*A*E)/(70*dx^2);
             0,  -(9*A*E)/(70*dx^2),                  0, 0,   (9*A*E)/(70*dx^2),   -(3*A*E)/(70*dx)];
            
            
            T3(:,:,3,5) =[ 0,                  0,                0, 0,                  0,          0;
             0, -(9*A*E)/(70*dx^2), -(3*A*E)/(70*dx), 0,  (9*A*E)/(70*dx^2),          0;
             0,   -(3*A*E)/(70*dx),        (A*E)/280, 0,    (3*A*E)/(70*dx), -(A*E)/280;
             0,                  0,                0, 0,                  0,          0;
             0,  (9*A*E)/(70*dx^2),  (3*A*E)/(70*dx), 0, -(9*A*E)/(70*dx^2),          0;
             0,                  0,       -(A*E)/280, 0,                  0, -(A*E)/280];
            
            
           
            
            T3(:,:,5,5) =[ 0,                   0,                  0, 0,                   0,                  0;
             0,  (36*A*E)/(35*dx^3),  (9*A*E)/(70*dx^2), 0, -(36*A*E)/(35*dx^3),  (9*A*E)/(70*dx^2);
             0,   (9*A*E)/(70*dx^2),    (3*A*E)/(70*dx), 0,  -(9*A*E)/(70*dx^2),                  0;
             0,                   0,                  0, 0,                   0,                  0;
             0, -(36*A*E)/(35*dx^3), -(9*A*E)/(70*dx^2), 0,  (36*A*E)/(35*dx^3), -(9*A*E)/(70*dx^2);
             0,   (9*A*E)/(70*dx^2),                  0, 0,  -(9*A*E)/(70*dx^2),    (3*A*E)/(70*dx)];
            
            
            T3(:,:,6,5) =[ 0,                  0,          0, 0,                  0,                0;
             0, -(9*A*E)/(70*dx^2),          0, 0,  (9*A*E)/(70*dx^2), -(3*A*E)/(70*dx);
             0,                  0, -(A*E)/280, 0,                  0,       -(A*E)/280;
             0,                  0,          0, 0,                  0,                0;
             0,  (9*A*E)/(70*dx^2),          0, 0, -(9*A*E)/(70*dx^2),  (3*A*E)/(70*dx);
             0,   -(3*A*E)/(70*dx), -(A*E)/280, 0,    (3*A*E)/(70*dx),        (A*E)/280];
            
            
            
            T3(:,:,2,6) =[ 0,                  0,         0, 0,                  0,                0;
             0,  (9*A*E)/(70*dx^2),         0, 0, -(9*A*E)/(70*dx^2),  (3*A*E)/(70*dx);
             0,                  0, (A*E)/280, 0,                  0,        (A*E)/280;
             0,                  0,         0, 0,                  0,                0;
             0, -(9*A*E)/(70*dx^2),         0, 0,  (9*A*E)/(70*dx^2), -(3*A*E)/(70*dx);
             0,    (3*A*E)/(70*dx), (A*E)/280, 0,   -(3*A*E)/(70*dx),       -(A*E)/280];
            
            
            T3(:,:,3,6) =[ 0,         0,             0, 0,          0,             0;
             0,         0,     (A*E)/280, 0,          0,     (A*E)/280;
             0, (A*E)/280, -(A*E*dx)/280, 0, -(A*E)/280,  (A*E*dx)/420;
             0,         0,             0, 0,          0,             0;
             0,         0,    -(A*E)/280, 0,          0,    -(A*E)/280;
             0, (A*E)/280,  (A*E*dx)/420, 0, -(A*E)/280, -(A*E*dx)/280];
            
            
            
            T3(:,:,5,6) =[ 0,                  0,          0, 0,                  0,                0;
             0, -(9*A*E)/(70*dx^2),          0, 0,  (9*A*E)/(70*dx^2), -(3*A*E)/(70*dx);
             0,                  0, -(A*E)/280, 0,                  0,       -(A*E)/280;
             0,                  0,          0, 0,                  0,                0;
             0,  (9*A*E)/(70*dx^2),          0, 0, -(9*A*E)/(70*dx^2),  (3*A*E)/(70*dx);
             0,   -(3*A*E)/(70*dx), -(A*E)/280, 0,    (3*A*E)/(70*dx),        (A*E)/280];
            
            
            T3(:,:,6,6) =[ 0,                0,             0, 0,                0,             0;
             0,  (3*A*E)/(70*dx),     (A*E)/280, 0, -(3*A*E)/(70*dx),    -(A*E)/280;
             0,        (A*E)/280,  (A*E*dx)/420, 0,       -(A*E)/280, -(A*E*dx)/280;
             0,                0,             0, 0,                0,             0;
             0, -(3*A*E)/(70*dx),    -(A*E)/280, 0,  (3*A*E)/(70*dx),     (A*E)/280;
             0,       -(A*E)/280, -(A*E*dx)/280, 0,        (A*E)/280,   (A*E*dx)/35];
            
        end
        
        function dS = dS(self,x)
            l = self.dx;
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            
            u1 = x(1);
            w1 = x(2);
            t1 = x(3);
            u2 = x(4);
            w2 = x(5);
            t2 = x(6);
            dS = [                                                0,                                                                                                                   -(A*E*(36*w1 - 36*w2 + 3*l*t1 + 3*l*t2))/(30*l^2),                                                                                                                       -(A*E*(3*w1 - 3*w2 + 4*l*t1 - l*t2))/(30*l),                                                0,                                                                                                                    (A*E*(36*w1 - 36*w2 + 3*l*t1 + 3*l*t2))/(30*l^2),                                                                                                                       -(A*E*(3*w1 - 3*w2 - l*t1 + 4*l*t2))/(30*l);
                -(A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2),  (A*E*(336*l*u2 - 336*l*u1 - 1728*w1*w2 + 864*w1^2 + 864*w2^2 + 36*l^2*t1^2 + 36*l^2*t2^2 + 216*l*t1*w1 - 216*l*t1*w2 + 216*l*t2*w1 - 216*l*t2*w2))/(280*l^3),                      (A*E*(28*l*u2 - 28*l*u1 - 216*w1*w2 + 108*w1^2 + 108*w2^2 - 3*l^2*t1^2 + 3*l^2*t2^2 + 6*l^2*t1*t2 + 72*l*t1*w1 - 72*l*t1*w2))/(280*l^2),  (A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2), -(A*E*(336*l*u2 - 336*l*u1 - 1728*w1*w2 + 864*w1^2 + 864*w2^2 + 36*l^2*t1^2 + 36*l^2*t2^2 + 216*l*t1*w1 - 216*l*t1*w2 + 216*l*t2*w1 - 216*l*t2*w2))/(280*l^3),                      (A*E*(28*l*u2 - 28*l*u1 - 216*w1*w2 + 108*w1^2 + 108*w2^2 + 3*l^2*t1^2 - 3*l^2*t2^2 + 6*l^2*t1*t2 + 72*l*t2*w1 - 72*l*t2*w2))/(280*l^2);
                -(A*E*(3*w1 - 3*w2 + 4*l*t1 - l*t2))/(30*l),                     (A*E*(84*l*u2 - 84*l*u1 - 648*w1*w2 + 324*w1^2 + 324*w2^2 - 9*l^2*t1^2 + 9*l^2*t2^2 + 18*l^2*t1*t2 + 216*l*t1*w1 - 216*l*t1*w2))/(840*l^2), (A*E*(56*l*u2 - 56*l*u1 - 108*w1*w2 + 54*w1^2 + 54*w2^2 + 36*l^2*t1^2 + 3*l^2*t2^2 - 9*l^2*t1*t2 - 9*l*t1*w1 + 9*l*t1*w2 + 9*l*t2*w1 - 9*l*t2*w2))/(420*l),    (A*E*(3*w1 - 3*w2 + 4*l*t1 - l*t2))/(30*l),                    -(A*E*(84*l*u2 - 84*l*u1 - 648*w1*w2 + 324*w1^2 + 324*w2^2 - 9*l^2*t1^2 + 9*l^2*t2^2 + 18*l^2*t1*t2 + 216*l*t1*w1 - 216*l*t1*w2))/(840*l^2),                                                         -(A*E*(28*u2 - 28*u1 - 18*t1*w1 + 18*t1*w2 - 18*t2*w1 + 18*t2*w2 + 9*l*t1^2 + 9*l*t2^2 - 12*l*t1*t2))/840;
                0,                                                                                                                    (A*E*(36*w1 - 36*w2 + 3*l*t1 + 3*l*t2))/(30*l^2),                                                                                                                        (A*E*(3*w1 - 3*w2 + 4*l*t1 - l*t2))/(30*l),                                                0,                                                                                                                   -(A*E*(36*w1 - 36*w2 + 3*l*t1 + 3*l*t2))/(30*l^2),                                                                                                                        (A*E*(3*w1 - 3*w2 - l*t1 + 4*l*t2))/(30*l);
                (A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2), -(A*E*(336*l*u2 - 336*l*u1 - 1728*w1*w2 + 864*w1^2 + 864*w2^2 + 36*l^2*t1^2 + 36*l^2*t2^2 + 216*l*t1*w1 - 216*l*t1*w2 + 216*l*t2*w1 - 216*l*t2*w2))/(280*l^3),                     -(A*E*(28*l*u2 - 28*l*u1 - 216*w1*w2 + 108*w1^2 + 108*w2^2 - 3*l^2*t1^2 + 3*l^2*t2^2 + 6*l^2*t1*t2 + 72*l*t1*w1 - 72*l*t1*w2))/(280*l^2), -(A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2),  (A*E*(336*l*u2 - 336*l*u1 - 1728*w1*w2 + 864*w1^2 + 864*w2^2 + 36*l^2*t1^2 + 36*l^2*t2^2 + 216*l*t1*w1 - 216*l*t1*w2 + 216*l*t2*w1 - 216*l*t2*w2))/(280*l^3),                     -(A*E*(28*l*u2 - 28*l*u1 - 216*w1*w2 + 108*w1^2 + 108*w2^2 + 3*l^2*t1^2 - 3*l^2*t2^2 + 6*l^2*t1*t2 + 72*l*t2*w1 - 72*l*t2*w2))/(280*l^2);
                -(A*E*(3*w1 - 3*w2 - l*t1 + 4*l*t2))/(30*l),                     (A*E*(84*l*u2 - 84*l*u1 - 648*w1*w2 + 324*w1^2 + 324*w2^2 + 9*l^2*t1^2 - 9*l^2*t2^2 + 18*l^2*t1*t2 + 216*l*t2*w1 - 216*l*t2*w2))/(840*l^2),                                                         -(A*E*(28*u2 - 28*u1 - 18*t1*w1 + 18*t1*w2 - 18*t2*w1 + 18*t2*w2 + 9*l*t1^2 + 9*l*t2^2 - 12*l*t1*t2))/840,    (A*E*(3*w1 - 3*w2 - l*t1 + 4*l*t2))/(30*l),                    -(A*E*(84*l*u2 - 84*l*u1 - 648*w1*w2 + 324*w1^2 + 324*w2^2 + 9*l^2*t1^2 - 9*l^2*t2^2 + 18*l^2*t1*t2 + 216*l*t2*w1 - 216*l*t2*w2))/(840*l^2), (A*E*(56*l*u2 - 56*l*u1 - 108*w1*w2 + 54*w1^2 + 54*w2^2 + 3*l^2*t1^2 + 36*l^2*t2^2 - 9*l^2*t1*t2 + 9*l*t1*w1 - 9*l*t1*w2 - 9*l*t2*w1 + 9*l*t2*w2))/(420*l)];
        end
        function dS_2 = dS_2(self,x)
            l = self.dx;
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            
            u1 = x(1);
            w1 = x(2);
            t1 = x(3);
            u2 = x(4);
            w2 = x(5);
            t2 = x(6);
            dS_2 =[                                                0, - (3*A*E*(2*w1 - 2*w2))/(5*l^2) - (A*E*(t1 + t2))/(10*l), - (A*E*(4*t1 - t2))/30 - (A*E*(w1 - w2))/(10*l),                                                0,   (3*A*E*(2*w1 - 2*w2))/(5*l^2) + (A*E*(t1 + t2))/(10*l), (A*E*(t1 - 4*t2))/30 - (A*E*(w1 - w2))/(10*l);
                -(A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2),                                -(6*A*E*(u1 - u2))/(5*l^2),                         -(A*E*(u1 - u2))/(10*l),  (A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2),                                 (6*A*E*(u1 - u2))/(5*l^2),                       -(A*E*(u1 - u2))/(10*l);
                - (A*E*(4*t1 - t2))/30 - (A*E*(w1 - w2))/(10*l),                                   -(A*E*(u1 - u2))/(10*l),                            -(2*A*E*(u1 - u2))/15,   (A*E*(4*t1 - t2))/30 + (A*E*(w1 - w2))/(10*l),                                    (A*E*(u1 - u2))/(10*l),                             (A*E*(u1 - u2))/30;
                0,   (3*A*E*(2*w1 - 2*w2))/(5*l^2) + (A*E*(t1 + t2))/(10*l),   (A*E*(4*t1 - t2))/30 + (A*E*(w1 - w2))/(10*l),                                                0, - (3*A*E*(2*w1 - 2*w2))/(5*l^2) - (A*E*(t1 + t2))/(10*l), (A*E*(w1 - w2))/(10*l) - (A*E*(t1 - 4*t2))/30;
                (A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2),                                 (6*A*E*(u1 - u2))/(5*l^2),                          (A*E*(u1 - u2))/(10*l), -(A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2),                                -(6*A*E*(u1 - u2))/(5*l^2),                        (A*E*(u1 - u2))/(10*l);
                (A*E*(t1 - 4*t2))/30 - (A*E*(w1 - w2))/(10*l),                                   -(A*E*(u1 - u2))/(10*l),                               (A*E*(u1 - u2))/30,   (A*E*(w1 - w2))/(10*l) - (A*E*(t1 - 4*t2))/30,                                    (A*E*(u1 - u2))/(10*l),                          -(2*A*E*(u1 - u2))/15];
        end
        
        function dS_3 = dS_3(self,x)
            l = self.dx;
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            
            w1 = x(2);
            t1 = x(3);
            w2 = x(5);
            t2 = x(6);
            dS_3 = [ 0,                                                                                                                                              0,                                                                                                                                           0, 0,                                                                                                                                              0,                                                                                                                                           0;
                0,  (A*E*(36*l^2*t1^2 + 36*l^2*t2^2 + 216*l*t1*w1 - 216*l*t1*w2 + 216*l*t2*w1 - 216*l*t2*w2 + 864*w1^2 - 1728*w1*w2 + 864*w2^2))/(280*l^3),                      (3*A*E*(- l^2*t1^2 + 2*l^2*t1*t2 + l^2*t2^2 + 24*l*t1*w1 - 24*l*t1*w2 + 36*w1^2 - 72*w1*w2 + 36*w2^2))/(280*l^2), 0, -(A*E*(36*l^2*t1^2 + 36*l^2*t2^2 + 216*l*t1*w1 - 216*l*t1*w2 + 216*l*t2*w1 - 216*l*t2*w2 + 864*w1^2 - 1728*w1*w2 + 864*w2^2))/(280*l^3),                        (3*A*E*(l^2*t1^2 + 2*l^2*t1*t2 - l^2*t2^2 + 24*l*t2*w1 - 24*l*t2*w2 + 36*w1^2 - 72*w1*w2 + 36*w2^2))/(280*l^2);
                0,                    (A*E*(- 3*l^2*t1^2 + 6*l^2*t1*t2 + 3*l^2*t2^2 + 72*l*t1*w1 - 72*l*t1*w2 + 108*w1^2 - 216*w1*w2 + 108*w2^2))/(280*l^2), (A*E*(12*l^2*t1^2 - 3*l^2*t1*t2 + l^2*t2^2 - 3*l*t1*w1 + 3*l*t1*w2 + 3*l*t2*w1 - 3*l*t2*w2 + 18*w1^2 - 36*w1*w2 + 18*w2^2))/(140*l), 0,                   -(A*E*(- 3*l^2*t1^2 + 6*l^2*t1*t2 + 3*l^2*t2^2 + 72*l*t1*w1 - 72*l*t1*w2 + 108*w1^2 - 216*w1*w2 + 108*w2^2))/(280*l^2),                                                     -(A*E*(6*t1*w2 - 6*t1*w1 - 6*t2*w1 + 6*t2*w2 + 3*l*t1^2 + 3*l*t2^2 - 4*l*t1*t2))/280;
                0,                                                                                                                                              0,                                                                                                                                           0, 0,                                                                                                                                              0,                                                                                                                                           0;
                0, -(A*E*(36*l^2*t1^2 + 36*l^2*t2^2 + 216*l*t1*w1 - 216*l*t1*w2 + 216*l*t2*w1 - 216*l*t2*w2 + 864*w1^2 - 1728*w1*w2 + 864*w2^2))/(280*l^3),                     -(3*A*E*(- l^2*t1^2 + 2*l^2*t1*t2 + l^2*t2^2 + 24*l*t1*w1 - 24*l*t1*w2 + 36*w1^2 - 72*w1*w2 + 36*w2^2))/(280*l^2), 0,  (A*E*(36*l^2*t1^2 + 36*l^2*t2^2 + 216*l*t1*w1 - 216*l*t1*w2 + 216*l*t2*w1 - 216*l*t2*w2 + 864*w1^2 - 1728*w1*w2 + 864*w2^2))/(280*l^3),                       -(3*A*E*(l^2*t1^2 + 2*l^2*t1*t2 - l^2*t2^2 + 24*l*t2*w1 - 24*l*t2*w2 + 36*w1^2 - 72*w1*w2 + 36*w2^2))/(280*l^2);
                0,                      (A*E*(3*l^2*t1^2 + 6*l^2*t1*t2 - 3*l^2*t2^2 + 72*l*t2*w1 - 72*l*t2*w2 + 108*w1^2 - 216*w1*w2 + 108*w2^2))/(280*l^2),                                                     -(A*E*(6*t1*w2 - 6*t1*w1 - 6*t2*w1 + 6*t2*w2 + 3*l*t1^2 + 3*l*t2^2 - 4*l*t1*t2))/280, 0,                     -(A*E*(3*l^2*t1^2 + 6*l^2*t1*t2 - 3*l^2*t2^2 + 72*l*t2*w1 - 72*l*t2*w2 + 108*w1^2 - 216*w1*w2 + 108*w2^2))/(280*l^2), (A*E*(l^2*t1^2 - 3*l^2*t1*t2 + 12*l^2*t2^2 + 3*l*t1*w1 - 3*l*t1*w2 - 3*l*t2*w1 + 3*l*t2*w2 + 18*w1^2 - 36*w1*w2 + 18*w2^2))/(140*l)];
            
        end
        
    end
    
end
