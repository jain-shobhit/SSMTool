%========================================================================
% DESCRIPTION: 
% Mechanical system comprised of a horizontal chain of point masses 
% connected via springs, dampers and nonlinear elements.
% 
% The constructor determines the system matrices M,D,K.
% 
% INPUT SPECIFICATION
%       mi      n-vector of mass values, from left to right
%       ki      n+1-vector of spring values, including the left and right
%               spring to the ground; if there is no spring to the ground,
%               set the value to zero
%       di      damper values analogous to ki
%========================================================================
% This file is part of NLvib.
% 
% If you use NLvib, please refer to the book:
%   M. Krack, J. Gross: Harmonic Balance for Nonlinear Vibration
%   Problems. Springer, 2019. https://doi.org/10.1007/978-3-030-14023-6.
% 
% COPYRIGHT AND LICENSING: 
% NLvib Version 1.1 Copyright (C) 2019  Malte Krack  
%										(malte.krack@ila.uni-stuttgart.de) 
%                     					Johann Gross 
%										(johann.gross@ila.uni-stuttgart.de)
%                     					University of Stuttgart
% This program comes with ABSOLUTELY NO WARRANTY. 
% NLvib is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%========================================================================
classdef ChainOfOscillators < MechanicalSystem
    methods
        function obj = ChainOfOscillators(mi,di,ki,...
                nonlinear_elements,Fex1)
            %% Constructor
            
            % Setup system matrices
            M = diag(mi);
            K = diag(ki(1:end-1)+ki(2:end)) - ...
                diag(ki(2:end-1),-1) - diag(ki(2:end-1),1);
            D = diag(di(1:end-1)+di(2:end)) - ...
                diag(di(2:end-1),-1) - diag(di(2:end-1),1);
            
            % Handle incomplete input
            if nargin<=4
                Fex1 = zeros(length(mi),1);
            end
            
            % Call parent class constructor
            obj = obj@MechanicalSystem(M,D,K,...
                nonlinear_elements,Fex1);
        end
    end
end