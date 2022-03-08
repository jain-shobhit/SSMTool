classdef KirchoffMaterial < Material
    properties 
        % Along with standard properties for the Material class, we have:
        
        PLANE_STRESS = true % whether the plane stess assumption is applied
                            % for 2-D continuum
        % Lami's constants 
        lambda                 
        mu
    end
    methods
        function self = KirchoffMaterial(varargin)
            % call Material Class constructor
            self@Material(varargin{:})
        end
        
        function lambda = get.lambda(self)
            E = self.YOUNGS_MODULUS;
            nu = self.POISSONS_RATIO;
            lambda = nu*E / ((1 + nu) * (1 - 2*nu));
        end
        
        function mu = get.mu(self)
            E = self.YOUNGS_MODULUS;
            nu = self.POISSONS_RATIO;
            mu  = E / (2*(1 + nu));
        end
        
        function D = get_stress_strain_matrix_2D(self)
            E = self.YOUNGS_MODULUS;
            nu = self.POISSONS_RATIO;

            if self.PLANE_STRESS
                D = E/(1-nu^2)*[1   nu  0;
                                nu  1   0;
                                0   0   (1-nu)/2];
            else
                D = [self.lambda + 2*self.mu    self.lambda,               0;
                    self.lambda                 self.lambda + 2*self.mu    0;
                    0                           0                          self.mu];
            end

        end
        
        function D = get_stress_strain_matrix_3D(self)
            D =[self.lambda + 2*self.mu, self.lambda, self.lambda, 0, 0, 0;
                self.lambda, self.lambda + 2*self.mu, self.lambda, 0, 0, 0;
                self.lambda, self.lambda, self.lambda + 2*self.mu, 0, 0, 0;
                0, 0, 0, self.mu, 0, 0;
                0, 0, 0, 0, self.mu, 0;
                0, 0, 0, 0, 0, self.mu];
        
        end
    end
end