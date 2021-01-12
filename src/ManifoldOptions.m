classdef ManifoldOptions < matlab.mixin.SetGet
    %ManifoldOptions This class sets the different options available for
    %the Manifold class
    
    properties
        notation  = 'tensor' % 'multiindex'
        paramStyle = 'normalform' % 'graph'
        reltol = 0.5
        IRtol = 0.05        
        
        contribNonAuto = true      % true: non-autonomous contributions are computed in manifold parametrization
                                   % false: non-autonomous contibutions not computed in the manifold parametrization 
                                   % and only leading-order reduced
                                   % dynamics parametrization is computed.
                                   % (relevant in FRC computation)
    end
    
    methods
        function set.notation(obj,notation)
            switch lower(notation)
                case 'tensor'
                    obj.notation = 'tensor';
                case 'multiindex'
                    obj.notation = 'multiindex';
                otherwise
                    error('Unknown notation type: set tensor or multiindex notation types')
            end
        end
        
        
        function set.contribNonAuto(obj,contribNonAuto)
            validateattributes(contribNonAuto,{'logical'},{'nonempty'})
            obj.contribNonAuto = contribNonAuto;
        end 
        
    end
end

