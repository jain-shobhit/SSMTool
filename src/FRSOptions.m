classdef FRSOptions < matlab.mixin.SetGet
    %ManifoldOptions This class sets the different options available for
    %the Manifold class

    properties        
        rhoMax   = 1                   % number of discrete rho values in the range [0, rhomax] (relevant for method == 'level set') 
        meshDens = 100                 % number of discrete parameter (Omega/epsilon) values in the parameter range         
        method   = 'level set'         % 'level set', 'continuation' 
        calFRS   = false               % whether compute FRS
    end

    methods                
        function set.method(obj,method)
            switch lower(method)
                case 'level set'
                    obj.method = 'level set';
                case 'continuation'
                    obj.method = 'continuation';
                otherwise
                    error('Unknown "method": set "level set" / "continuation" types')
            end
        end        
    end
end