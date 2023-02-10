function y = compute_auto_invariance_error(obj,W,R,rhosamp,orders,ntheta,varargin)
% COMPUTE_AUTO_INVARIANCE_ERROR This function calculates the invariance
% error measure for obtained SSM (and reduced dynamics on it). It is the
% average of residuals of the invariance equation evaluated at a collection
% of points on the SSM. This implementation only support 2D and 4D SSMs.
% For 2D SSM, the sampling points are given by
%           p = [\rho*exp(i*theta_j) \rho*exp(-i*theta_j)], 1<=j<ntheta
% where {\theta_j} follows a regular mesh.
% For 4D SSM, the samping points are given by
%  p = [\rho_1*exp(i*theta^1_j) \rho_1*exp(-i*theta^1_j)
%       \rho_2*exp(i*theta^2_k) \rho_2*exp(-i*theta^2_k)]
% where \rho_1=\rho\cos\alpha_i, \rho_2=\rho\sin\alpha_i, 1<=i<=nalpha
%
% The errors will be computed at each radius in rhosamp and each expansion
% order in orders. The maximum element in orders should not be larger than
% the order of expansion in W and R.
%
% Y = compute_auto_invariance_error(obj,W,R,rhosamp,orders,ntheta,varargin)
%
% W: expansion of autonomous SSM
% R: expansion of reduced dynamics on SSM
% rhosamp: sampling radius for parameterization coordinates, at which error
% measures are computed
% orders: ordera at which invariance errors are computed
% ntheta: number of sample points along angle directions
% varargin: nalpha if the SSM is four-dimensional

dimSSM = obj.dimManifold;
norder = numel(orders);
nrho   = numel(rhosamp);
y      = zeros(norder,nrho);
theta  = linspace(0,2*pi,ntheta+1);
theta  = theta(1:end-1);
assert(max(orders)<=numel(R), 'Some requested expansion orders are higher than the approximation order of SSM');
assert(dimSSM==2||dimSSM==4, 'Only support 2D and 4D SSMs');
if dimSSM==4
    nalpha = varargin{1}; 
    alphas = linspace(0,pi/2,nalpha); 
end
% loop over orders
for ko=1:norder
    orderk = orders(ko);
    Wk = W(1:orderk);
    Rk = R(1:orderk);
    % loop over rhosamp
    for krho=1:nrho
        fprintf('calculate error at order %d and rho=%d\n',orderk,rhosamp(krho));
        if dimSSM==2
            % construct sampling points
            pj  = rhosamp(krho)*exp(1i*theta);
            % evaluate residuals
            res = obj.compuate_invariance_residual(Wk,Rk,pj,'auto');
        else
            res = zeros(nalpha,ntheta,ntheta);
            % construct sampling points
            for i=1:nalpha
                rho1 = rhosamp(krho)*cos(alphas(i));
                rho2 = rhosamp(krho)*sin(alphas(i));
                for j=1:ntheta
                    p1 = rho1*exp(1i*theta(j))*ones(1,ntheta);
                    p3 = rho2*exp(1i*theta);
                    pij = [p1;conj(p1);p3;conj(p3)];
                    res(i,j,:) = obj.compuate_invariance_residual(Wk,Rk,pij,'auto');
                end                
            end
            res = res(:);
        end
        % record avearged residual
        y(ko,krho) = mean(res);
    end
end

end
        
        
        