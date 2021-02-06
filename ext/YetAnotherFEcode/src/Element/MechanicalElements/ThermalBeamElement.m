classdef ThermalBeamElement < Element
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
        function self = ThermalBeamElement(b, h, Material)
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
        function [xe,Te] = extract_element_data(self,x,T)
            % x is a vector of full DOFs
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            xe = x(index,:);
            Te = T(self.nodeIDs,:);
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
        
        function [K, F] = tangent_stiffness_and_force(self,x,t)
            % this function computes the element stiffness matrix and
            % internal force vector in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            % x : vector of full DOFs in global coordinates
            
            [x_e, t_e] = self.extract_element_data(x,t);
            T_e = self.transformationMatrix;
            % Displacements in local coordinates
            q = T_e*x_e;
            % Forces and Stiffness in local coordinates
            [K_local, F_local] = self.tangent_stiffness_and_force_local(q,t_e);
            % Forces and Stiffness in global coordinates
            K = T_e.' * K_local * T_e;
            F = T_e.' * F_local;
        end
        
        function [F] = internal_force(self,x,t)
            % this function computes the element stiffness matrix and
            % internal force vector in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            % x : vector of full DOFs in global coordinates
            
            [x_e, t_e] = self.extract_element_data(x,t);
            T_e = self.transformationMatrix;
            % Displacements in local coordinates
            q = T_e*x_e;
            % Forces and Stiffness in local coordinates
            [F_local] = self.internal_force_local(q,t_e);
            % Force in global coordinates
            F = T_e.' * F_local;
        end
        
        function [K] = stiffness_derivative(self,x,t,v)
            % this function computes the element stiffness derivative matrix
            % in the global coordinates when the
            % nodes : matrix containing Nodal coordinates
            % x : element DOF values in global coordinates
            [x_e, t_e] = self.extract_element_data(x,t);
            [v_e, ~] = self.extract_element_data(v,t);
            
            T_e = self.transformationMatrix;
            
            q = T_e * x_e;
            eta = T_e * v_e;
            [K_local] = self.stiffness_derivative_local(q,t_e,eta);
            K = Te.' * K_local * Te;
        end
        
        
        %% Local level methods
        
        function K = tangent_stiffness_matrix_local(self,x,T)
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            l = self.dx;
            I = self.areaMoment;
            alpha = self.Material.ThermalExpansionCoefficient;
            
            u1 = x(1);
            w1 = x(2);
            t1 = x(3);
            u2 = x(4);
            w2 = x(5);
            t2 = x(6);
            T1 = T(1);
            T2 = T(2);
            
            K =  [                                         (A*E)/l,                                                                                                                                                                                    -((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                          - (A*E*(4*t1 - t2))/30 - (A*E*(3*w1 - 3*w2))/(30*l),                                        -(A*E)/l,                                                                                                                                                                                     ((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                            (A*E*(t1 - 4*t2))/30 - (A*E*(3*w1 - 3*w2))/(30*l);
                -(A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2),  (E*(3360*I + 864*A*w1^2 + 864*A*w2^2 + 36*A*l^2*t1^2 + 36*A*l^2*t2^2 - 336*A*l*u1 + 336*A*l*u2 - 1728*A*w1*w2 - 168*A*T1*alpha*l^2 - 168*A*T2*alpha*l^2 + 216*A*l*t1*w1 - 216*A*l*t1*w2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(280*l^3),                                               (E*(1680*I + 108*A*w1^2 + 108*A*w2^2 - 3*A*l^2*t1^2 + 3*A*l^2*t2^2 - 28*A*l*u1 + 28*A*l*u2 - 216*A*w1*w2 - 28*A*T2*alpha*l^2 + 6*A*l^2*t1*t2 + 72*A*l*t1*w1 - 72*A*l*t1*w2))/(280*l^2),  (A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2), -(E*(3360*I + 864*A*w1^2 + 864*A*w2^2 + 36*A*l^2*t1^2 + 36*A*l^2*t2^2 - 336*A*l*u1 + 336*A*l*u2 - 1728*A*w1*w2 - 168*A*T1*alpha*l^2 - 168*A*T2*alpha*l^2 + 216*A*l*t1*w1 - 216*A*l*t1*w2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(280*l^3),                                               (E*(1680*I + 108*A*w1^2 + 108*A*w2^2 + 3*A*l^2*t1^2 - 3*A*l^2*t2^2 - 28*A*l*u1 + 28*A*l*u2 - 216*A*w1*w2 - 28*A*T1*alpha*l^2 + 6*A*l^2*t1*t2 + 72*A*l*t2*w1 - 72*A*l*t2*w2))/(280*l^2);
                -(A*E*(3*w1 - 3*w2 + 4*l*t1 - l*t2))/(30*l),                                              (E*(5040*I + 324*A*w1^2 + 324*A*w2^2 - 9*A*l^2*t1^2 + 9*A*l^2*t2^2 - 84*A*l*u1 + 84*A*l*u2 - 648*A*w1*w2 - 84*A*T2*alpha*l^2 + 18*A*l^2*t1*t2 + 216*A*l*t1*w1 - 216*A*l*t1*w2))/(840*l^2), (E*(1680*I + 54*A*w1^2 + 54*A*w2^2 + 36*A*l^2*t1^2 + 3*A*l^2*t2^2 - 56*A*l*u1 + 56*A*l*u2 - 108*A*w1*w2 - 42*A*T1*alpha*l^2 - 14*A*T2*alpha*l^2 - 9*A*l^2*t1*t2 - 9*A*l*t1*w1 + 9*A*l*t1*w2 + 9*A*l*t2*w1 - 9*A*l*t2*w2))/(420*l),    (A*E*(3*w1 - 3*w2 + 4*l*t1 - l*t2))/(30*l),                                             -(E*(5040*I + 324*A*w1^2 + 324*A*w2^2 - 9*A*l^2*t1^2 + 9*A*l^2*t2^2 - 84*A*l*u1 + 84*A*l*u2 - 648*A*w1*w2 - 84*A*T2*alpha*l^2 + 18*A*l^2*t1*t2 + 216*A*l*t1*w1 - 216*A*l*t1*w2))/(840*l^2),                                   (E*(1680*I - 9*A*l^2*t1^2 - 9*A*l^2*t2^2 + 28*A*l*u1 - 28*A*l*u2 + 14*A*T1*alpha*l^2 + 14*A*T2*alpha*l^2 + 12*A*l^2*t1*t2 + 18*A*l*t1*w1 - 18*A*l*t1*w2 + 18*A*l*t2*w1 - 18*A*l*t2*w2))/(840*l);
                -(A*E)/l,                                                                                                                                                                                     ((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                            (A*E*(4*t1 - t2))/30 + (A*E*(3*w1 - 3*w2))/(30*l),                                         (A*E)/l,                                                                                                                                                                                    -((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                            (A*E*(3*w1 - 3*w2))/(30*l) - (A*E*(t1 - 4*t2))/30;
                (A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2), -(E*(3360*I + 864*A*w1^2 + 864*A*w2^2 + 36*A*l^2*t1^2 + 36*A*l^2*t2^2 - 336*A*l*u1 + 336*A*l*u2 - 1728*A*w1*w2 - 168*A*T1*alpha*l^2 - 168*A*T2*alpha*l^2 + 216*A*l*t1*w1 - 216*A*l*t1*w2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(280*l^3),                                              -(E*(1680*I + 108*A*w1^2 + 108*A*w2^2 - 3*A*l^2*t1^2 + 3*A*l^2*t2^2 - 28*A*l*u1 + 28*A*l*u2 - 216*A*w1*w2 - 28*A*T2*alpha*l^2 + 6*A*l^2*t1*t2 + 72*A*l*t1*w1 - 72*A*l*t1*w2))/(280*l^2), -(A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2),  (E*(3360*I + 864*A*w1^2 + 864*A*w2^2 + 36*A*l^2*t1^2 + 36*A*l^2*t2^2 - 336*A*l*u1 + 336*A*l*u2 - 1728*A*w1*w2 - 168*A*T1*alpha*l^2 - 168*A*T2*alpha*l^2 + 216*A*l*t1*w1 - 216*A*l*t1*w2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(280*l^3),                                              -(E*(1680*I + 108*A*w1^2 + 108*A*w2^2 + 3*A*l^2*t1^2 - 3*A*l^2*t2^2 - 28*A*l*u1 + 28*A*l*u2 - 216*A*w1*w2 - 28*A*T1*alpha*l^2 + 6*A*l^2*t1*t2 + 72*A*l*t2*w1 - 72*A*l*t2*w2))/(280*l^2);
                -(A*E*(3*w1 - 3*w2 - l*t1 + 4*l*t2))/(30*l),                                              (E*(5040*I + 324*A*w1^2 + 324*A*w2^2 + 9*A*l^2*t1^2 - 9*A*l^2*t2^2 - 84*A*l*u1 + 84*A*l*u2 - 648*A*w1*w2 - 84*A*T1*alpha*l^2 + 18*A*l^2*t1*t2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(840*l^2),                                   (E*(1680*I - 9*A*l^2*t1^2 - 9*A*l^2*t2^2 + 28*A*l*u1 - 28*A*l*u2 + 14*A*T1*alpha*l^2 + 14*A*T2*alpha*l^2 + 12*A*l^2*t1*t2 + 18*A*l*t1*w1 - 18*A*l*t1*w2 + 18*A*l*t2*w1 - 18*A*l*t2*w2))/(840*l),    (A*E*(3*w1 - 3*w2 - l*t1 + 4*l*t2))/(30*l),                                             -(E*(5040*I + 324*A*w1^2 + 324*A*w2^2 + 9*A*l^2*t1^2 - 9*A*l^2*t2^2 - 84*A*l*u1 + 84*A*l*u2 - 648*A*w1*w2 - 84*A*T1*alpha*l^2 + 18*A*l^2*t1*t2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(840*l^2), (E*(1680*I + 54*A*w1^2 + 54*A*w2^2 + 3*A*l^2*t1^2 + 36*A*l^2*t2^2 - 56*A*l*u1 + 56*A*l*u2 - 108*A*w1*w2 - 14*A*T1*alpha*l^2 - 42*A*T2*alpha*l^2 - 9*A*l^2*t1*t2 + 9*A*l*t1*w1 - 9*A*l*t1*w2 - 9*A*l*t2*w1 + 9*A*l*t2*w2))/(420*l)];
            
        end
        
        function F = internal_force_local(self,x,T)
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            l = self.dx;
            I = self.areaMoment;
            alpha = self.Material.THERMAL_EXPANSION_COEFFICIENT;
            
            u1 = x(1);
            w1 = x(2);
            t1 = x(3);
            u2 = x(4);
            w2 = x(5);
            t2 = x(6);
            T1 = T(1);
            T2 = T(2);
            
            F =   [(A*E*(- 2*t1^2 + t1*t2 - 2*t2^2 + 15*T1*alpha + 15*T2*alpha))/30 - ((A*E*(18*w1^2 - 36*w1*w2 + 18*w2^2))/30 - (A*E*l*(30*u1 - 30*u2 - 3*t1*w1 + 3*t1*w2 - 3*t2*w1 + 3*t2*w2))/30)/l^2;
                (E*(3360*I*w1 - 3360*I*w2 + 288*A*w1^3 - 288*A*w2^3 - A*l^3*t1^3 - A*l^3*t2^3 + 1680*I*l*t1 + 1680*I*l*t2 + 864*A*w1*w2^2 - 864*A*w1^2*w2 - 28*A*l^2*t1*u1 + 28*A*l^2*t1*u2 - 28*A*l^2*t2*u1 + 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 + 3*A*l^3*t1*t2^2 + 3*A*l^3*t1^2*t2 + 36*A*l^2*t1^2*w1 - 36*A*l^2*t1^2*w2 + 36*A*l^2*t2^2*w1 - 36*A*l^2*t2^2*w2 - 336*A*l*u1*w1 + 336*A*l*u1*w2 + 336*A*l*u2*w1 - 336*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 216*A*l*t2*w1*w2 - 28*A*T1*alpha*l^3*t2 - 28*A*T2*alpha*l^3*t1 - 168*A*T1*alpha*l^2*w1 + 168*A*T1*alpha*l^2*w2 - 168*A*T2*alpha*l^2*w1 + 168*A*T2*alpha*l^2*w2))/(280*l^3);
                (E*(5040*I*w1 - 5040*I*w2 + 108*A*w1^3 - 108*A*w2^3 + 24*A*l^3*t1^3 - 3*A*l^3*t2^3 + 3360*I*l*t1 + 1680*I*l*t2 + 324*A*w1*w2^2 - 324*A*w1^2*w2 - 112*A*l^2*t1*u1 + 112*A*l^2*t1*u2 + 28*A*l^2*t2*u1 - 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 6*A*l^3*t1*t2^2 - 9*A*l^3*t1^2*t2 - 9*A*l^2*t1^2*w1 + 9*A*l^2*t1^2*w2 + 9*A*l^2*t2^2*w1 - 9*A*l^2*t2^2*w2 - 84*A*l*u1*w1 + 84*A*l*u1*w2 + 84*A*l*u2*w1 - 84*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 84*A*T1*alpha*l^3*t1 + 14*A*T1*alpha*l^3*t2 - 28*A*T2*alpha*l^3*t1 + 14*A*T2*alpha*l^3*t2 - 84*A*T2*alpha*l^2*w1 + 84*A*T2*alpha*l^2*w2 + 18*A*l^2*t1*t2*w1 - 18*A*l^2*t1*t2*w2))/(840*l^2);
                ((A*E*(18*w1^2 - 36*w1*w2 + 18*w2^2))/30 - (A*E*l*(30*u1 - 30*u2 - 3*t1*w1 + 3*t1*w2 - 3*t2*w1 + 3*t2*w2))/30)/l^2 - (A*E*(- 2*t1^2 + t1*t2 - 2*t2^2 + 15*T1*alpha + 15*T2*alpha))/30;
                -(E*(3360*I*w1 - 3360*I*w2 + 288*A*w1^3 - 288*A*w2^3 - A*l^3*t1^3 - A*l^3*t2^3 + 1680*I*l*t1 + 1680*I*l*t2 + 864*A*w1*w2^2 - 864*A*w1^2*w2 - 28*A*l^2*t1*u1 + 28*A*l^2*t1*u2 - 28*A*l^2*t2*u1 + 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 + 3*A*l^3*t1*t2^2 + 3*A*l^3*t1^2*t2 + 36*A*l^2*t1^2*w1 - 36*A*l^2*t1^2*w2 + 36*A*l^2*t2^2*w1 - 36*A*l^2*t2^2*w2 - 336*A*l*u1*w1 + 336*A*l*u1*w2 + 336*A*l*u2*w1 - 336*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 216*A*l*t2*w1*w2 - 28*A*T1*alpha*l^3*t2 - 28*A*T2*alpha*l^3*t1 - 168*A*T1*alpha*l^2*w1 + 168*A*T1*alpha*l^2*w2 - 168*A*T2*alpha*l^2*w1 + 168*A*T2*alpha*l^2*w2))/(280*l^3);
                (E*(5040*I*w1 - 5040*I*w2 + 108*A*w1^3 - 108*A*w2^3 - 3*A*l^3*t1^3 + 24*A*l^3*t2^3 + 1680*I*l*t1 + 3360*I*l*t2 + 324*A*w1*w2^2 - 324*A*w1^2*w2 + 28*A*l^2*t1*u1 - 28*A*l^2*t1*u2 - 112*A*l^2*t2*u1 + 112*A*l^2*t2*u2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 - 9*A*l^3*t1*t2^2 + 6*A*l^3*t1^2*t2 + 9*A*l^2*t1^2*w1 - 9*A*l^2*t1^2*w2 - 9*A*l^2*t2^2*w1 + 9*A*l^2*t2^2*w2 - 84*A*l*u1*w1 + 84*A*l*u1*w2 + 84*A*l*u2*w1 - 84*A*l*u2*w2 - 216*A*l*t2*w1*w2 + 14*A*T1*alpha*l^3*t1 - 28*A*T1*alpha*l^3*t2 + 14*A*T2*alpha*l^3*t1 - 84*A*T2*alpha*l^3*t2 - 84*A*T1*alpha*l^2*w1 + 84*A*T1*alpha*l^2*w2 + 18*A*l^2*t1*t2*w1 - 18*A*l^2*t1*t2*w2))/(840*l^2)];
            
        end
        
        function [K,F] = tangent_stiffness_and_force_local(self,x,T)
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            l = self.dx;
            I = self.areaMoment;
            alpha = self.Material.THERMAL_EXPANSION_COEFFICIENT;
            
            u1 = x(1);
            w1 = x(2);
            t1 = x(3);
            u2 = x(4);
            w2 = x(5);
            t2 = x(6);
            T1 = T(1);
            T2 = T(2);
            
            
            K =  [                                         (A*E)/l,                                                                                                                                                                                    -((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                          - (A*E*(4*t1 - t2))/30 - (A*E*(3*w1 - 3*w2))/(30*l),                                        -(A*E)/l,                                                                                                                                                                                     ((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                            (A*E*(t1 - 4*t2))/30 - (A*E*(3*w1 - 3*w2))/(30*l);
                -(A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2),  (E*(3360*I + 864*A*w1^2 + 864*A*w2^2 + 36*A*l^2*t1^2 + 36*A*l^2*t2^2 - 336*A*l*u1 + 336*A*l*u2 - 1728*A*w1*w2 - 168*A*T1*alpha*l^2 - 168*A*T2*alpha*l^2 + 216*A*l*t1*w1 - 216*A*l*t1*w2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(280*l^3),                                               (E*(1680*I + 108*A*w1^2 + 108*A*w2^2 - 3*A*l^2*t1^2 + 3*A*l^2*t2^2 - 28*A*l*u1 + 28*A*l*u2 - 216*A*w1*w2 - 28*A*T2*alpha*l^2 + 6*A*l^2*t1*t2 + 72*A*l*t1*w1 - 72*A*l*t1*w2))/(280*l^2),  (A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2), -(E*(3360*I + 864*A*w1^2 + 864*A*w2^2 + 36*A*l^2*t1^2 + 36*A*l^2*t2^2 - 336*A*l*u1 + 336*A*l*u2 - 1728*A*w1*w2 - 168*A*T1*alpha*l^2 - 168*A*T2*alpha*l^2 + 216*A*l*t1*w1 - 216*A*l*t1*w2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(280*l^3),                                               (E*(1680*I + 108*A*w1^2 + 108*A*w2^2 + 3*A*l^2*t1^2 - 3*A*l^2*t2^2 - 28*A*l*u1 + 28*A*l*u2 - 216*A*w1*w2 - 28*A*T1*alpha*l^2 + 6*A*l^2*t1*t2 + 72*A*l*t2*w1 - 72*A*l*t2*w2))/(280*l^2);
                -(A*E*(3*w1 - 3*w2 + 4*l*t1 - l*t2))/(30*l),                                              (E*(5040*I + 324*A*w1^2 + 324*A*w2^2 - 9*A*l^2*t1^2 + 9*A*l^2*t2^2 - 84*A*l*u1 + 84*A*l*u2 - 648*A*w1*w2 - 84*A*T2*alpha*l^2 + 18*A*l^2*t1*t2 + 216*A*l*t1*w1 - 216*A*l*t1*w2))/(840*l^2), (E*(1680*I + 54*A*w1^2 + 54*A*w2^2 + 36*A*l^2*t1^2 + 3*A*l^2*t2^2 - 56*A*l*u1 + 56*A*l*u2 - 108*A*w1*w2 - 42*A*T1*alpha*l^2 - 14*A*T2*alpha*l^2 - 9*A*l^2*t1*t2 - 9*A*l*t1*w1 + 9*A*l*t1*w2 + 9*A*l*t2*w1 - 9*A*l*t2*w2))/(420*l),    (A*E*(3*w1 - 3*w2 + 4*l*t1 - l*t2))/(30*l),                                             -(E*(5040*I + 324*A*w1^2 + 324*A*w2^2 - 9*A*l^2*t1^2 + 9*A*l^2*t2^2 - 84*A*l*u1 + 84*A*l*u2 - 648*A*w1*w2 - 84*A*T2*alpha*l^2 + 18*A*l^2*t1*t2 + 216*A*l*t1*w1 - 216*A*l*t1*w2))/(840*l^2),                                   (E*(1680*I - 9*A*l^2*t1^2 - 9*A*l^2*t2^2 + 28*A*l*u1 - 28*A*l*u2 + 14*A*T1*alpha*l^2 + 14*A*T2*alpha*l^2 + 12*A*l^2*t1*t2 + 18*A*l*t1*w1 - 18*A*l*t1*w2 + 18*A*l*t2*w1 - 18*A*l*t2*w2))/(840*l);
                -(A*E)/l,                                                                                                                                                                                     ((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                            (A*E*(4*t1 - t2))/30 + (A*E*(3*w1 - 3*w2))/(30*l),                                         (A*E)/l,                                                                                                                                                                                    -((A*E*(36*w1 - 36*w2))/30 + (A*E*l*(3*t1 + 3*t2))/30)/l^2,                                                                                                                                                                                            (A*E*(3*w1 - 3*w2))/(30*l) - (A*E*(t1 - 4*t2))/30;
                (A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2), -(E*(3360*I + 864*A*w1^2 + 864*A*w2^2 + 36*A*l^2*t1^2 + 36*A*l^2*t2^2 - 336*A*l*u1 + 336*A*l*u2 - 1728*A*w1*w2 - 168*A*T1*alpha*l^2 - 168*A*T2*alpha*l^2 + 216*A*l*t1*w1 - 216*A*l*t1*w2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(280*l^3),                                              -(E*(1680*I + 108*A*w1^2 + 108*A*w2^2 - 3*A*l^2*t1^2 + 3*A*l^2*t2^2 - 28*A*l*u1 + 28*A*l*u2 - 216*A*w1*w2 - 28*A*T2*alpha*l^2 + 6*A*l^2*t1*t2 + 72*A*l*t1*w1 - 72*A*l*t1*w2))/(280*l^2), -(A*E*(12*w1 - 12*w2 + l*t1 + l*t2))/(10*l^2),  (E*(3360*I + 864*A*w1^2 + 864*A*w2^2 + 36*A*l^2*t1^2 + 36*A*l^2*t2^2 - 336*A*l*u1 + 336*A*l*u2 - 1728*A*w1*w2 - 168*A*T1*alpha*l^2 - 168*A*T2*alpha*l^2 + 216*A*l*t1*w1 - 216*A*l*t1*w2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(280*l^3),                                              -(E*(1680*I + 108*A*w1^2 + 108*A*w2^2 + 3*A*l^2*t1^2 - 3*A*l^2*t2^2 - 28*A*l*u1 + 28*A*l*u2 - 216*A*w1*w2 - 28*A*T1*alpha*l^2 + 6*A*l^2*t1*t2 + 72*A*l*t2*w1 - 72*A*l*t2*w2))/(280*l^2);
                -(A*E*(3*w1 - 3*w2 - l*t1 + 4*l*t2))/(30*l),                                              (E*(5040*I + 324*A*w1^2 + 324*A*w2^2 + 9*A*l^2*t1^2 - 9*A*l^2*t2^2 - 84*A*l*u1 + 84*A*l*u2 - 648*A*w1*w2 - 84*A*T1*alpha*l^2 + 18*A*l^2*t1*t2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(840*l^2),                                   (E*(1680*I - 9*A*l^2*t1^2 - 9*A*l^2*t2^2 + 28*A*l*u1 - 28*A*l*u2 + 14*A*T1*alpha*l^2 + 14*A*T2*alpha*l^2 + 12*A*l^2*t1*t2 + 18*A*l*t1*w1 - 18*A*l*t1*w2 + 18*A*l*t2*w1 - 18*A*l*t2*w2))/(840*l),    (A*E*(3*w1 - 3*w2 - l*t1 + 4*l*t2))/(30*l),                                             -(E*(5040*I + 324*A*w1^2 + 324*A*w2^2 + 9*A*l^2*t1^2 - 9*A*l^2*t2^2 - 84*A*l*u1 + 84*A*l*u2 - 648*A*w1*w2 - 84*A*T1*alpha*l^2 + 18*A*l^2*t1*t2 + 216*A*l*t2*w1 - 216*A*l*t2*w2))/(840*l^2), (E*(1680*I + 54*A*w1^2 + 54*A*w2^2 + 3*A*l^2*t1^2 + 36*A*l^2*t2^2 - 56*A*l*u1 + 56*A*l*u2 - 108*A*w1*w2 - 14*A*T1*alpha*l^2 - 42*A*T2*alpha*l^2 - 9*A*l^2*t1*t2 + 9*A*l*t1*w1 - 9*A*l*t1*w2 - 9*A*l*t2*w1 + 9*A*l*t2*w2))/(420*l)];
            
            F =   [(A*E*(- 2*t1^2 + t1*t2 - 2*t2^2 + 15*T1*alpha + 15*T2*alpha))/30 - ((A*E*(18*w1^2 - 36*w1*w2 + 18*w2^2))/30 - (A*E*l*(30*u1 - 30*u2 - 3*t1*w1 + 3*t1*w2 - 3*t2*w1 + 3*t2*w2))/30)/l^2;
                (E*(3360*I*w1 - 3360*I*w2 + 288*A*w1^3 - 288*A*w2^3 - A*l^3*t1^3 - A*l^3*t2^3 + 1680*I*l*t1 + 1680*I*l*t2 + 864*A*w1*w2^2 - 864*A*w1^2*w2 - 28*A*l^2*t1*u1 + 28*A*l^2*t1*u2 - 28*A*l^2*t2*u1 + 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 + 3*A*l^3*t1*t2^2 + 3*A*l^3*t1^2*t2 + 36*A*l^2*t1^2*w1 - 36*A*l^2*t1^2*w2 + 36*A*l^2*t2^2*w1 - 36*A*l^2*t2^2*w2 - 336*A*l*u1*w1 + 336*A*l*u1*w2 + 336*A*l*u2*w1 - 336*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 216*A*l*t2*w1*w2 - 28*A*T1*alpha*l^3*t2 - 28*A*T2*alpha*l^3*t1 - 168*A*T1*alpha*l^2*w1 + 168*A*T1*alpha*l^2*w2 - 168*A*T2*alpha*l^2*w1 + 168*A*T2*alpha*l^2*w2))/(280*l^3);
                (E*(5040*I*w1 - 5040*I*w2 + 108*A*w1^3 - 108*A*w2^3 + 24*A*l^3*t1^3 - 3*A*l^3*t2^3 + 3360*I*l*t1 + 1680*I*l*t2 + 324*A*w1*w2^2 - 324*A*w1^2*w2 - 112*A*l^2*t1*u1 + 112*A*l^2*t1*u2 + 28*A*l^2*t2*u1 - 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 6*A*l^3*t1*t2^2 - 9*A*l^3*t1^2*t2 - 9*A*l^2*t1^2*w1 + 9*A*l^2*t1^2*w2 + 9*A*l^2*t2^2*w1 - 9*A*l^2*t2^2*w2 - 84*A*l*u1*w1 + 84*A*l*u1*w2 + 84*A*l*u2*w1 - 84*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 84*A*T1*alpha*l^3*t1 + 14*A*T1*alpha*l^3*t2 - 28*A*T2*alpha*l^3*t1 + 14*A*T2*alpha*l^3*t2 - 84*A*T2*alpha*l^2*w1 + 84*A*T2*alpha*l^2*w2 + 18*A*l^2*t1*t2*w1 - 18*A*l^2*t1*t2*w2))/(840*l^2);
                ((A*E*(18*w1^2 - 36*w1*w2 + 18*w2^2))/30 - (A*E*l*(30*u1 - 30*u2 - 3*t1*w1 + 3*t1*w2 - 3*t2*w1 + 3*t2*w2))/30)/l^2 - (A*E*(- 2*t1^2 + t1*t2 - 2*t2^2 + 15*T1*alpha + 15*T2*alpha))/30;
                -(E*(3360*I*w1 - 3360*I*w2 + 288*A*w1^3 - 288*A*w2^3 - A*l^3*t1^3 - A*l^3*t2^3 + 1680*I*l*t1 + 1680*I*l*t2 + 864*A*w1*w2^2 - 864*A*w1^2*w2 - 28*A*l^2*t1*u1 + 28*A*l^2*t1*u2 - 28*A*l^2*t2*u1 + 28*A*l^2*t2*u2 + 108*A*l*t1*w1^2 + 108*A*l*t1*w2^2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 + 3*A*l^3*t1*t2^2 + 3*A*l^3*t1^2*t2 + 36*A*l^2*t1^2*w1 - 36*A*l^2*t1^2*w2 + 36*A*l^2*t2^2*w1 - 36*A*l^2*t2^2*w2 - 336*A*l*u1*w1 + 336*A*l*u1*w2 + 336*A*l*u2*w1 - 336*A*l*u2*w2 - 216*A*l*t1*w1*w2 - 216*A*l*t2*w1*w2 - 28*A*T1*alpha*l^3*t2 - 28*A*T2*alpha*l^3*t1 - 168*A*T1*alpha*l^2*w1 + 168*A*T1*alpha*l^2*w2 - 168*A*T2*alpha*l^2*w1 + 168*A*T2*alpha*l^2*w2))/(280*l^3);
                (E*(5040*I*w1 - 5040*I*w2 + 108*A*w1^3 - 108*A*w2^3 - 3*A*l^3*t1^3 + 24*A*l^3*t2^3 + 1680*I*l*t1 + 3360*I*l*t2 + 324*A*w1*w2^2 - 324*A*w1^2*w2 + 28*A*l^2*t1*u1 - 28*A*l^2*t1*u2 - 112*A*l^2*t2*u1 + 112*A*l^2*t2*u2 + 108*A*l*t2*w1^2 + 108*A*l*t2*w2^2 - 9*A*l^3*t1*t2^2 + 6*A*l^3*t1^2*t2 + 9*A*l^2*t1^2*w1 - 9*A*l^2*t1^2*w2 - 9*A*l^2*t2^2*w1 + 9*A*l^2*t2^2*w2 - 84*A*l*u1*w1 + 84*A*l*u1*w2 + 84*A*l*u2*w1 - 84*A*l*u2*w2 - 216*A*l*t2*w1*w2 + 14*A*T1*alpha*l^3*t1 - 28*A*T1*alpha*l^3*t2 + 14*A*T2*alpha*l^3*t1 - 84*A*T2*alpha*l^3*t2 - 84*A*T1*alpha*l^2*w1 + 84*A*T1*alpha*l^2*w2 + 18*A*l^2*t1*t2*w1 - 18*A*l^2*t1*t2*w2))/(840*l^2)];
            
        end
        
        function [K] = stiffness_derivative_local(self,x,~,v)
            E = self.Material.YOUNGS_MODULUS;
            A = self.area;
            l = self.dx;
            
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
            
            K =  [                                                    0,                                                                                                                                                                           -((A*E*(36*vw1 - 36*vw2))/30 + (A*E*l*(3*vt1 + 3*vt2))/30)/l^2,                                                                                                                                                                                                                     - (A*E*(4*vt1 - vt2))/30 - (A*E*(3*vw1 - 3*vw2))/(30*l),                                                    0,                                                                                                                                                                            ((A*E*(36*vw1 - 36*vw2))/30 + (A*E*l*(3*vt1 + 3*vt2))/30)/l^2,                                                                                                                                                                                                                       (A*E*(vt1 - 4*vt2))/30 - (A*E*(3*vw1 - 3*vw2))/(30*l);
                -(A*E*(12*vw1 - 12*vw2 + l*vt1 + l*vt2))/(10*l^2),  (3*A*E*(14*l*vu2 - 14*l*vu1 + 72*vw1*w1 - 72*vw1*w2 - 72*vw2*w1 + 72*vw2*w2 + 3*l^2*t1*vt1 + 3*l^2*t2*vt2 + 9*l*t1*vw1 - 9*l*t1*vw2 + 9*l*t2*vw1 - 9*l*t2*vw2 + 9*l*vt1*w1 - 9*l*vt1*w2 + 9*l*vt2*w1 - 9*l*vt2*w2))/(35*l^3),                                                     (A*E*(14*l*vu2 - 14*l*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 - 3*l^2*t1*vt1 + 3*l^2*t1*vt2 + 3*l^2*t2*vt1 + 3*l^2*t2*vt2 + 36*l*t1*vw1 - 36*l*t1*vw2 + 36*l*vt1*w1 - 36*l*vt1*w2))/(140*l^2),  (A*E*(12*vw1 - 12*vw2 + l*vt1 + l*vt2))/(10*l^2), -(3*A*E*(14*l*vu2 - 14*l*vu1 + 72*vw1*w1 - 72*vw1*w2 - 72*vw2*w1 + 72*vw2*w2 + 3*l^2*t1*vt1 + 3*l^2*t2*vt2 + 9*l*t1*vw1 - 9*l*t1*vw2 + 9*l*t2*vw1 - 9*l*t2*vw2 + 9*l*vt1*w1 - 9*l*vt1*w2 + 9*l*vt2*w1 - 9*l*vt2*w2))/(35*l^3),                                                     (A*E*(14*l*vu2 - 14*l*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 + 3*l^2*t1*vt1 + 3*l^2*t1*vt2 + 3*l^2*t2*vt1 - 3*l^2*t2*vt2 + 36*l*t2*vw1 - 36*l*t2*vw2 + 36*l*vt2*w1 - 36*l*vt2*w2))/(140*l^2);
                -(A*E*(3*vw1 - 3*vw2 + 4*l*vt1 - l*vt2))/(30*l),                   (A*E*(14*l*vu2 - 14*l*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 - 3*l^2*t1*vt1 + 3*l^2*t1*vt2 + 3*l^2*t2*vt1 + 3*l^2*t2*vt2 + 36*l*t1*vw1 - 36*l*t1*vw2 + 36*l*vt1*w1 - 36*l*vt1*w2))/(140*l^2), -(A*E*(56*l*vu1 - 56*l*vu2 - 108*vw1*w1 + 108*vw1*w2 + 108*vw2*w1 - 108*vw2*w2 - 72*l^2*t1*vt1 + 9*l^2*t1*vt2 + 9*l^2*t2*vt1 - 6*l^2*t2*vt2 + 9*l*t1*vw1 - 9*l*t1*vw2 - 9*l*t2*vw1 + 9*l*t2*vw2 + 9*l*vt1*w1 - 9*l*vt1*w2 - 9*l*vt2*w1 + 9*l*vt2*w2))/(420*l),    (A*E*(3*vw1 - 3*vw2 + 4*l*vt1 - l*vt2))/(30*l),                  -(A*E*(14*l*vu2 - 14*l*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 - 3*l^2*t1*vt1 + 3*l^2*t1*vt2 + 3*l^2*t2*vt1 + 3*l^2*t2*vt2 + 36*l*t1*vw1 - 36*l*t1*vw2 + 36*l*vt1*w1 - 36*l*vt1*w2))/(140*l^2),                                                                                                  (A*E*(14*vu1 - 14*vu2 + 9*t1*vw1 - 9*t1*vw2 + 9*t2*vw1 - 9*t2*vw2 + 9*vt1*w1 - 9*vt1*w2 + 9*vt2*w1 - 9*vt2*w2 - 9*l*t1*vt1 + 6*l*t1*vt2 + 6*l*t2*vt1 - 9*l*t2*vt2))/420;
                0,                                                                                                                                                                            ((A*E*(36*vw1 - 36*vw2))/30 + (A*E*l*(3*vt1 + 3*vt2))/30)/l^2,                                                                                                                                                                                                                       (A*E*(4*vt1 - vt2))/30 + (A*E*(3*vw1 - 3*vw2))/(30*l),                                                    0,                                                                                                                                                                           -((A*E*(36*vw1 - 36*vw2))/30 + (A*E*l*(3*vt1 + 3*vt2))/30)/l^2,                                                                                                                                                                                                                       (A*E*(3*vw1 - 3*vw2))/(30*l) - (A*E*(vt1 - 4*vt2))/30;
                (A*E*(12*vw1 - 12*vw2 + l*vt1 + l*vt2))/(10*l^2), -(3*A*E*(14*l*vu2 - 14*l*vu1 + 72*vw1*w1 - 72*vw1*w2 - 72*vw2*w1 + 72*vw2*w2 + 3*l^2*t1*vt1 + 3*l^2*t2*vt2 + 9*l*t1*vw1 - 9*l*t1*vw2 + 9*l*t2*vw1 - 9*l*t2*vw2 + 9*l*vt1*w1 - 9*l*vt1*w2 + 9*l*vt2*w1 - 9*l*vt2*w2))/(35*l^3),                                                    -(A*E*(14*l*vu2 - 14*l*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 - 3*l^2*t1*vt1 + 3*l^2*t1*vt2 + 3*l^2*t2*vt1 + 3*l^2*t2*vt2 + 36*l*t1*vw1 - 36*l*t1*vw2 + 36*l*vt1*w1 - 36*l*vt1*w2))/(140*l^2), -(A*E*(12*vw1 - 12*vw2 + l*vt1 + l*vt2))/(10*l^2),  (3*A*E*(14*l*vu2 - 14*l*vu1 + 72*vw1*w1 - 72*vw1*w2 - 72*vw2*w1 + 72*vw2*w2 + 3*l^2*t1*vt1 + 3*l^2*t2*vt2 + 9*l*t1*vw1 - 9*l*t1*vw2 + 9*l*t2*vw1 - 9*l*t2*vw2 + 9*l*vt1*w1 - 9*l*vt1*w2 + 9*l*vt2*w1 - 9*l*vt2*w2))/(35*l^3),                                                    -(A*E*(14*l*vu2 - 14*l*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 + 3*l^2*t1*vt1 + 3*l^2*t1*vt2 + 3*l^2*t2*vt1 - 3*l^2*t2*vt2 + 36*l*t2*vw1 - 36*l*t2*vw2 + 36*l*vt2*w1 - 36*l*vt2*w2))/(140*l^2);
                -(A*E*(3*vw1 - 3*vw2 - l*vt1 + 4*l*vt2))/(30*l),                   (A*E*(14*l*vu2 - 14*l*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 + 3*l^2*t1*vt1 + 3*l^2*t1*vt2 + 3*l^2*t2*vt1 - 3*l^2*t2*vt2 + 36*l*t2*vw1 - 36*l*t2*vw2 + 36*l*vt2*w1 - 36*l*vt2*w2))/(140*l^2),                                                                                                  (A*E*(14*vu1 - 14*vu2 + 9*t1*vw1 - 9*t1*vw2 + 9*t2*vw1 - 9*t2*vw2 + 9*vt1*w1 - 9*vt1*w2 + 9*vt2*w1 - 9*vt2*w2 - 9*l*t1*vt1 + 6*l*t1*vt2 + 6*l*t2*vt1 - 9*l*t2*vt2))/420,    (A*E*(3*vw1 - 3*vw2 - l*vt1 + 4*l*vt2))/(30*l),                  -(A*E*(14*l*vu2 - 14*l*vu1 + 108*vw1*w1 - 108*vw1*w2 - 108*vw2*w1 + 108*vw2*w2 + 3*l^2*t1*vt1 + 3*l^2*t1*vt2 + 3*l^2*t2*vt1 - 3*l^2*t2*vt2 + 36*l*t2*vw1 - 36*l*t2*vw2 + 36*l*vt2*w1 - 36*l*vt2*w2))/(140*l^2), -(A*E*(56*l*vu1 - 56*l*vu2 - 108*vw1*w1 + 108*vw1*w2 + 108*vw2*w1 - 108*vw2*w2 - 6*l^2*t1*vt1 + 9*l^2*t1*vt2 + 9*l^2*t2*vt1 - 72*l^2*t2*vt2 - 9*l*t1*vw1 + 9*l*t1*vw2 + 9*l*t2*vw1 - 9*l*t2*vw2 - 9*l*vt1*w1 + 9*l*vt1*w2 + 9*l*vt2*w1 - 9*l*vt2*w2))/(420*l)];
        end
        
        
        
    end
    
end
