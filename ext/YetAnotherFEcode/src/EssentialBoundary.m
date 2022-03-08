classdef EssentialBoundary < handle
    % class for defining Dirichlet boundary conditions
    
    properties
        nDOFs               % Total number of DOFs in the Mesh
        nDOFPerNode         % Number of DOFs per node in the Mesh
        constrainedDOFs     % 2 column matrix: 1st column contains the
                                % indices of Dirichlet DOFs, 2nd column
                                % contains the solution values at the
                                % corresponding DOFs
        unconstrainedDOFs   %  
        B                   % Boolean matrix to constrain Matrix to free
                                % DOFs after removing rows and columns with
                                % indices corresponding to the Dirichlet DOFs.
    end
    
    methods
        function self = EssentialBoundary(nDOFs,nDOFPerNode)
            self.nDOFs = nDOFs;
            self.nDOFPerNode = nDOFPerNode;
        end
        
        function apply_Dirichlet_BC(self,constrainedNodes,constrainedDOF,value)
            % DirichletNodes : vector containing node IDs
            % ConstrainedDOF : vector (or scalar) containing the index of the
            %                    nodal DOFs whose value is prescribed
            % value : scalar which gives the Dirichlet BC value
            
            if ~isscalar(value)
                error('Please specify a scalar value for Dirichlet boundary condition')
            end
            
            n = self.nDOFPerNode;
            nc = length(constrainedNodes);
            nd = length(constrainedDOF);
            
            if max(constrainedDOF) > n
                error('Wrong input for ConstrainedDOF')
            end
            
            % collect the DOFs for Dirichlet BC
            cDOFs = zeros(nc,nd);
            for j = 1:nd
                cDOFs(:,j) = n * (constrainedNodes - 1) + constrainedDOF(j); % overloading the + operator (adding scalar to vector)
            end
            
            % check previously exisiting Dirichlet DOFs
            if ~isempty(self.constrainedDOFs)
                [bool, loc] = ismember(cDOFs(:), self.constrainedDOFs(:,1));
                if any(bool)
                    warning('Some DOFs have previously been assigned values: These would be overwritten by the new value')
                    repeat = setdiff(loc,0);
                    self.constrainedDOFs(repeat,:) = []; % remove the conflicting Dirichlet DOFs from the list
                end
            end
            
            % add new Dirichlet DOFs to the Dirichlet array
            newDirichletDOFs  = unique(cDOFs(:));
            n = length(newDirichletDOFs);
            val = value*ones(n,1);
            newconstrainedDOFs = [newDirichletDOFs val];
            self.constrainedDOFs = [self.constrainedDOFs; newconstrainedDOFs];
            
            % update FreeDOFs
            self.unconstrainedDOFs = setdiff(1 : self.nDOFs, self.constrainedDOFs(:,1));
            
            % update the boolean matrix
            nf = length(self.unconstrainedDOFs);
            self.B = sparse(self.unconstrainedDOFs,1:nf,true,self.nDOFs,nf);
        end
        
        function [K, F] = apply(self,K,F)
            n = length(K);
            dofs = self.constrainedDOFs(:,1);
            vals = self.constrainedDOFs(:,2);
            % correct RHS
            F = F - K(:,dofs) * vals;
            % set to zero the rows and colums corresponding to Dirichlet DOFs
            K(dofs,:) = 0;
            K(:,dofs) = 0;
            % set to one the diagonal entries corresponding to Dirichlet DOFs
            diagonal_indices = (n+1)*(dofs-1)+1;
            K(diagonal_indices) = 1;
            % set to Dirichlet value the entries corresponding to Dirichlet
            % DOFs in F
            F(dofs) = vals;
        end       
        

    end
    
end

