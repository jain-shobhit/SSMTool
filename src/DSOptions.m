classdef DSOptions < matlab.mixin.SetGet
    %DSOptions Options for DynamicalSystems class
    properties
        notation  = 'tensor' % 'multiindex'
        Nmax = 100; % maximum dimensionality up to which all eigenvalues would be computed for first order systems
        Emax = 10;   % If all eigenvalues are not computed, then we use only first E_max eigenvalues for checking outer resonance.
        outDOF = []; % output degree of freedom
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
    end
end

