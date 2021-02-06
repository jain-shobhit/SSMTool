classdef (Abstract) Element < handle
    properties (Abstract)
        % every subclass must have the following properties
        nodes       % global coordinates of element nodes
        nodeIDs     % the index location of element nodes    
        nDOFPerNode
        nNodes      
        nDim        % dimensionality of the element 
                    % (may be different from spatial dimensionality of the mesh)
    end
    
    properties (Dependent)
        iDOFs       % indices of the global DOFs associated to the element 
    end
    
    methods (Abstract)
        % every subclass must implement the below methods on its own (in
        % addition to its own methods        
        [K,F] = tangent_stiffness_and_force(self,varargin)
        varargout    = extract_element_data(self,varargin)    % get the element unknowns from global unkowns
    end    
    
    methods
        function idx = get.iDOFs(self)
            idx = get_index(self.nodeIDs, self.nDOFPerNode);
        end
    end
end
