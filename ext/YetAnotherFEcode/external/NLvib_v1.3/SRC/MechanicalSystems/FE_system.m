%========================================================================
% DESCRIPTION: 
% Mechanical system governed by the equations of motion
% 
%   M*\ddot u + D*\dot u + K*u + fnl(u) = fex
% 
%       where fnl(u) is a polynomial in u (stiffness nonlinearity)
% 
%========================================================================
% This file is an extension of YetAnotherFEcode for NLvib.
% class: FE_system < MechanicalSystem
% constructor: obj = FE_system(myAssembly, Fext, nonlinearity_type)
%   where
%       myAssembly: is the Assembly from YetAnotherFEcode
%       Fext:       is the external forcing
%       nonlinearity_type: can be either "FE" or "custom". In the "FE" case
%                   only elastic internal forces are computed. With the
%                   "CUSTOM" version, one can define an ad hoc function to
%                   compute the nonlinear terms in fnl. See the examble
%                   below.
%
% *** EXAMPLE OF THE CUSTOM FUNCTION
% function [Kt, fnl, varargout] = fnl_CUSTOM(myAssembly, varargin)
%     % This function returns the nonlinear vector and its Jacobian. The
%     % default input is the Assembly of the system (see YetAnotherFEcode),
%     % but additional inputs can be considered. 
%     % NOTE 1: if additional inputs (outputs) are required (e.g. the 
%     % velocity vector), the "custom" caseS in HB_residual.m 
%     % (HB_nonlinear_forces_AFT function, global nonlinearity) and 
%     % shooting_residual.m (nonlinear_forces) MUST BE MODIFIED manually
%     % adding the additional input to fnl_CUSTOM. By default, fnl_CUSTOM
%     % is implemented in the aforementioned functions with the Assembly
%     % input and the displacement vector X. 
%     % (see lines ~400-410 in HB_residual.m)
%     % (see lines ~277-284 in shooting_residual.m)
%     % NOTE 2: fnl and Kt must contain ONLY the nonlinear terms! (not K0*X
%     % for intance).
%     
%     X = varargin{1};
%     ...
%     
%     % use custom class methods
%     Kt =  myAssembly.custom_Kt_method( X );
%     fnl = myAssembly.custom_fnl_method( X );
%     
%     % or define custom functions
%     Kt =  custom_Kt_function( X );
%     fnl = custom_fnl_function( X );
%     
%     varargout{1} = X; % (optional)
% end
%
% Created: 21 April 2021
% Jacopo Marconi, Politecnico di Milano
%========================================================================

classdef FE_system < MechanicalSystem
    methods
        
        % Constructor
        function obj = FE_system(myAssembly, Fext, nonlinearity_type)
            
            ntot = myAssembly.Mesh.nDOFs;
            n = length(myAssembly.Mesh.EBC.unconstrainedDOFs);
            
            % Handle incomplete input
            if nargin < 2
                Fext = zeros(n, 1);
                nonlinearity_type = 'fe';
            elseif nargin < 3
                nonlinearity_type = 'fe';
            elseif nargin == 3
                if ~(strcmp(nonlinearity_type, 'fe') || strcmp(nonlinearity_type, 'custom'))
                    error(['Wrong nonlinearity type. Available options '...
                        'are "fe" or "custom". For the "custom" version'...
                        ',specific functions must be defined outside of'...
                        ' this class.'])
                end
            end
            
            % K: stiffness matrix, M: mass matrix, D: damping matrix
            K = myAssembly.DATA.K;
            M = myAssembly.DATA.M;
            D = myAssembly.DATA.D;
            
            % check sizes (constrained matrices are used)        
            if size(K,1) == ntot
                K = myAssembly.constrain_matrix(K);
            end
            if size(M,1) == ntot
                M = myAssembly.constrain_matrix(M);
            end
            if size(D,1) == ntot
                D = myAssembly.constrain_matrix(D);
            end
            if length(Fext) == ntot
                Fext = myAssembly.constrain_vector(Fext);
            end
            
            % Define nonlinearity
            nonlinearity = struct('type',nonlinearity_type,'K',K,'islocal',0,...
                'assembly', myAssembly);
            nonlinearity.force_direction = Fext*0+1; % for shooting method
            if strcmp(nonlinearity_type, 'custom')
                nonlinearity.custom_function = myAssembly.DATA.fnl_CUSTOM;
            end
            
            nonlinearity = { nonlinearity };
            
            % Call parent class constructor
            obj = obj@MechanicalSystem(M,D,K,nonlinearity,Fext);
        end
    end
end



