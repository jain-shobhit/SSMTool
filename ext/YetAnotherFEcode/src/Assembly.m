classdef Assembly < handle
    properties
        parallelized % whether to compute element contributions in parallel        
        Mesh    % Object of Mesh class (handle to it)
        DATA    % Miscellaneous data structure with arbitrary, user-defined 
                % fields, which may be useful to precompute and store in  
                % the Assembly Object e.g. DATA.K, DATA.M, DATA.C etc.
    end

    methods
        function self = Assembly(Mesh,varargin)
            self.Mesh = Mesh;   
            if nargin == 2
                self.parallelized = varargin{1};
            else
                self.parallelized = false;
            end
        end
        
        function set.parallelized(self,val)
            self.parallelized = val;
            adjust_parallel_settings(self);
        end
        
        function adjust_parallel_settings(self)
            ps = parallel.Settings;
            if self.parallelized
                ps.Pool.AutoCreate = true;
            else
                % delete exisiting parallel pool
                poolobj = gcp('nocreate');
                delete(poolobj);
                % prevent autocreation
                ps.Pool.AutoCreate = false;               
            end                
        end
        
        function [K, f] = tangent_stiffness_and_force(self, varargin)
            [K, f] = self.matrix_and_vector('tangent_stiffness_and_force',...
                        varargin{:});
        end

        function [f] = internal_force(self, varargin)
            f = self.vector('internal_force',varargin{:});
        end

        function [F] = uniform_body_force(self,varargin)
            n_e = self.Mesh.nElements;
            F = cell(n_e,1); % values
            index = cell(n_e,1);
            
            [elementWeights,~] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % extracting domain elements from the supplied set
            domainElements = ~[self.Mesh.Elements(elementSet).isBoundary];
            
            Elements = self.Mesh.Elements;
            parfor j = elementSet(domainElements)
                thisElement = Elements(j).Object;
                index{j} = thisElement.iDOFs;            
                F{j} = thisElement.uniformBodyForce;
            end
            index = vertcat(index{:});
            F = vertcat(F{:});
            F = sparse(index, ones(length(index),1), F, self.Mesh.nDOFs, 1);
        end

        function M = mass_matrix(self, varargin)
            M = self.matrix('mass_matrix',varargin{:});
        end

        function C = damping_matrix(self, varargin)
            C = self.matrix('damping_matrix',varargin{:});
        end

        function K = stiffness_matrix(self, varargin)
            K = self.matrix('stiffness_matrix',varargin{:});
        end

        function [Kd] = stiffness_derivative(self, varargin)
            Kd = self.matrix('stiffness_derivative',varargin{:});
        end

        function [K] = matrix(self,elementMethodName,varargin)
            % This function assembles a generic finite element matrix from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level matrix Ke.
            % For this to work, a method named elementMethodName which
            % returns the appropriate matrix must be defined for all
            % element types in the FE Mesh.

            n_e = self.Mesh.nElements;

            I = cell(n_e,1); % row indices
            J = cell(n_e,1); % column indices
            K = cell(n_e,1); % values
            Elements = self.Mesh.Elements;
            
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions            
            parfor j = elementSet
                thisElement = Elements(j).Object;
                
                index = thisElement.iDOFs;
                d = length(index);
                I{j} = kron(true(d,1), index);
                J{j} = kron(index, true(d,1));

                [Ke] = elementWeights(j) * thisElement.(elementMethodName)(inputs{:});
                K{j} = Ke(:);
            end

            I = vertcat(I{:});
            J = vertcat(J{:});
            K = vertcat(K{:});

            K = sparse(I, J, K, self.Mesh.nDOFs, self.Mesh.nDOFs);
        end

        function [F] = vector(self,elementMethodName,varargin)
            % This function assembles a generic finite element vector from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level vector Fe.
            % For this to work, a method named elementMethodName which
            % returns the appropriate vector must be defined for all
            % element types in the FE Mesh.

            n_e = self.Mesh.nElements;            
            index = cell(n_e,1); % indices
            F = cell(n_e,1); % values
            Elements = self.Mesh.Elements; % Elements array
            
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            parfor j = elementSet
                thisElement = Elements(j).Object;
                
                index{j} = thisElement.iDOFs;
                F{j} = elementWeights(j) * thisElement.(elementMethodName)(inputs{:});
            end
            
            % Assembling
            index = vertcat(index{:});
            F = vertcat(F{:});
            F = sparse(index, ones(length(index),1), F, self.Mesh.nDOFs, 1);

        end
        
        function [S] = scalar(self,elementMethodName,varargin)
            % This function assembles a generic finite element scalar from
            % its element level. 
            % elementMethodName is a string input containing the name of
            % the method that returns the element level scalar.
            % For this to work, a method named elementMethodName which
            % returns the appropriate scalar must be defined for all
            % element types in the FE Mesh.
            
            n_e = self.Mesh.nElements;
            
            S = zeros(n_e,1); % values
            Elements = self.Mesh.Elements; % Elements array            
            
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            parfor j = elementSet
                thisElement = Elements(j).Object;
                
                S(j) = elementWeights(j) * thisElement.(elementMethodName)(inputs{:});
            end            
        end
        
        function [K, f] = matrix_and_vector(self,elementMethodName, varargin)
            % This function assembles a generic finite element matrix and 
            % vector from its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level matrix Ke.
            % For this to work, a method named elementMethodName which
            % returns the appropriate matrix must be defined for all
            % element types in the FE Mesh.
            n_e = self.Mesh.nElements;
            
            index = cell(n_e,1);
            I = cell(n_e,1); % row indices
            J = cell(n_e,1); % column indices
            K = cell(n_e,1); % values
            f = cell(n_e,1); % values
            Elements = self.Mesh.Elements;
            
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            parfor j = elementSet
                thisElement = Elements(j).Object;
                
                index{j} = thisElement.iDOFs;
                d = length(index{j});
                I{j} = kron(true(d,1), index{j});
                J{j} = kron(index{j}, true(d,1));
                
                [Ke, fe] =  thisElement.(elementMethodName)(inputs{:});
                K{j} = elementWeights(j) * Ke(:);
                f{j} = elementWeights(j) * fe;
            end
            
            index = vertcat(index{:});
            I = vertcat(I{:});
            J = vertcat(J{:});
            K = vertcat(K{:});
            f = vertcat(f{:});
            
            K = sparse(I, J, K, self.Mesh.nDOFs, self.Mesh.nDOFs);
            f = sparse(index, ones(length(index),1), f, self.Mesh.nDOFs, 1);
        end
        
        function [T] = tensor(self,elementMethodName,SIZE,sumDIMS,varargin)
            % This function assembles a generic finite element vector from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level vector Fe.
            % For this to work, a method named elementMethodName which
            % returns the appropriate vector must be defined for all
            % element types in the FE Mesh.

            n_e = self.Mesh.nElements;

            subs = cell(n_e,1);
            T = cell(n_e,1); % values
            
            Elements = self.Mesh.Elements;
            
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            parfor j = elementSet
                thisElement = Elements(j).Object;
                
                [Te, globalSUBS] = thisElement.(elementMethodName)(inputs{:});

                [subs{j}, T{j}] = sparsify(elementWeights(j) * Te, globalSUBS, sumDIMS);
            end

            subs = vertcat(subs{:});
            if ~isempty(sumDIMS)
                subs(:,sumDIMS) = sort(subs(:,sumDIMS),2);
            end
            T = vertcat(T{:});
            T = sptensor(subs, T, SIZE);

        end


        function Kc = constrain_matrix(self,K)
            if ~isempty(self.Mesh.EBC)
                Kc = K(self.Mesh.EBC.unconstrainedDOFs,self.Mesh.EBC.unconstrainedDOFs);
            else
                Kc = K;
            end
        end

        function vc = constrain_vector(self,v)
            if ~isempty(self.Mesh.EBC)
                vc = v(self.Mesh.EBC.unconstrainedDOFs,:);
            else
                vc = v;
            end
        end
        
        function Tc = constrain_tensor(self,T,varargin)            
            if ~isempty(self.Mesh.EBC)
                narginchk(2,3)
                n = ndims(T);
                SIZE = size(T);
                if nargin>2
                    constrainedDIMS = varargin{1}; 
                else
                    constrainedDIMS = 1:n;
                end                
                unchangedDIMS = setdiff(1:n,constrainedDIMS);                
                dofs = self.Mesh.EBC.unconstrainedDOFs;
                subs = cell(1,n);
                subs(constrainedDIMS) = {dofs}; 
                for j = unchangedDIMS
                    subs{j} = 1:SIZE(j);
                end
                Tc = T(subs{:});                
            else
                Tc = T;
            end
        end

        function v = unconstrain_vector(self,vc)
            if ~isempty(self.Mesh.EBC)
                v = self.Mesh.EBC.B * vc;
                dofs = self.Mesh.EBC.constrainedDOFs(:,1);
                vals = self.Mesh.EBC.constrainedDOFs(:,2);
                v(dofs,:) = repmat(vals,[1,size(vc,2)]);
            else
                v = vc;
            end
        end
        
        function ind_cons = free2constrained_index(self, ind_free)
            % this function changes a DOF index, referred to the complete
            % DOF indexes (free system), into the corresponding index in
            % the constrained system. ind_free can also be a vector.
            % Example:  free system DOFs = [1 2 3 4]
            %           assume DOF #1 is constrained, we'll have:
            %           constrained system DOFs = [2 3 4].
            %           Then DOF#3 corresponds to   ind_free=3
            %                                       ind_cons=2
            ind_cons = zeros(size(ind_free));
            for ii = 1 : length(ind_free)
                ind_cons(ii) = find( self.Mesh.EBC.unconstrainedDOFs == ind_free(ii));
            end
        end

        function u = solve_system(self,K,f)
            if ~isempty(self.Mesh.EBC)
                [K_bc, f_bc] = self.Mesh.EBC.apply(K,f);
                u = K_bc\f_bc;
            else
                u = K\f;
            end
        end
        
        function [weights, input] = parse_inputs(obj,varargin)
            % this function parses the optional argument 'weights' which
            % defines the element weights to be used in FE Assembly. This
            % should be passed as a name value pair. The remaining arguments
            % after parsing weights are also returned.
            
            n_e = obj.Mesh.nElements;
            
            % locate position of 'weights' string
            iWeights = find(strcmpi(varargin,'weights'));            
            

            if ~isempty(iWeights)
                % extract user-supplied value of weights
                weights = varargin{iWeights + 1};
                
                % check input for weights
                assert( isvector(weights) && numel(weights) == n_e && ...
                    (isnumeric(weights) || islogical(weights)), ...
                    ['Incorrect input for weights: it should be a vector with ' ...
                    num2str(n_e) ' elements']);
                
                if iscolumn(weights)
                    weights = transpose(weights);
                end

                % extract the remaining input
                input = varargin([1:iWeights-1,iWeights+2:end]);
            else
                input = varargin;
                weights = true(1,n_e);
            end
        end
    end
end
