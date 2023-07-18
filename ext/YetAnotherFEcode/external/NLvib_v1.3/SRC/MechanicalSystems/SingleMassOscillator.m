%========================================================================
% DESCRIPTION: 
% Single mass oscillator with linear spring and damper, possible nonlinear
% elements and excitation.
%========================================================================
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
classdef SingleMassOscillator < ChainOfOscillators
    methods
        function obj = SingleMassOscillator(m,d,k,...
                nonlinear_elements,varargin)
            % Make sure force directions are given
            for i=1:length(nonlinear_elements)
                nonlinear_elements{i}.force_direction = 1;
            end
            %% A single mass oscillator is just a special case of a 
            % chain of oscillators (add zero damping and stiffness to the
            % right)
            obj = obj@ChainOfOscillators(m,[d 0],[k 0],...
                nonlinear_elements,varargin{:});
        end
    end
end