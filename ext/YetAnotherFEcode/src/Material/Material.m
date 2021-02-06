classdef Material < matlab.mixin.SetGet
    properties
        % SI units follows (fictitious material properties)
        YOUNGS_MODULUS = 70e9   % [Pa]
        POISSONS_RATIO = 0.3    % [-]
        THERMAL_EXPANSION_COEFFICIENT = 11.7e-6  % [K^(-1)]
        DENSITY = 2700         % [Kg/m^3]
        DAMPING_MODULUS = 1e8   % [Pa-s]
        THERMAL_CONDUCTIVITY = 205  % [W/m K]
        HEAT_CAPACITY = 900         % [J/kg K]
        CONVECTION_COEFFICIENT = 50 % [W/m^2 K]
        EMISSIVITY = 0             % [-] 
        % Feel free to add more material properties here with their default
        % values, allow for validation of user's input by defining set 
        % methods below
        
    end
    
    properties (Constant)
       STEFAN_BOLTZMANN = 5.670373e-8 
    end
    
    methods        
        function set.YOUNGS_MODULUS(self,val)
            validateattributes(val,{'numeric'},{'nonempty','positive'})
            self.YOUNGS_MODULUS = val;
        end
        
        function set.POISSONS_RATIO(self,val)
            validateattributes(val,{'numeric'},{'nonempty','positive'})
            self.POISSONS_RATIO = val;
        end
    end
end
