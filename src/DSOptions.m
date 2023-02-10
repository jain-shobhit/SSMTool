classdef DSOptions < matlab.mixin.SetGet
    %DSOptions Options for DynamicalSystems class
    properties
        notation  = 'tensor' % 'multiindex'
        Nmax = 100; % maximum dimensionality up to which all eigenvalues would be computed for first order systems
        Emax = 10;   % If all eigenvalues are not computed, then we use only first E_max eigenvalues for checking outer resonance.
        outDOF = []; % output degree of freedom
        RayleighDamping  = true; % damping matrix of second-order system
        HarmonicForce = true; % external forcing
        lambdaThreshold = 1e16; % Threshold for stiff eigenmodes (will be removed)
        BaseExcitation = false; % harmonic forcing in the form \epsilon\Omega^2 f^{ext}(\Omega t)
        sigma  = 0; % used in eigenvalue computation in eigs (compute Emax eigenvalues around sigma, which is zero by default - smallestabs)
        RemoveZeros = true; % remove zero eigenvalues in linear spectral analysis (false if parameter-dependent SSM/LSM is computed)
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

