function [A,B,F] = build_model(varargin)

% parameter set
gamma = 1;
beta  = 1.25;
b     = 0.5;
C     = 0.1;

% build Charney-DeVore model
m = 1:2;
alpham    = 8*sqrt(2)/pi*m.^2./(4*m.^2-1).*(b^2+m.^2-1)./(b^2+m.^2);
betam     = beta*b^2./(b^2+m.^2);
deltam    = 64*sqrt(2)/(15*pi)*(b^2-m.^2+1)./(b^2+m.^2);
gammaastm = gamma*4*m*sqrt(2)*b./(4*m.^2-1)/pi;
gammam    = gamma*4*m.^3./(4*m.^2-1)*sqrt(2)*b./pi./(b^2+m.^2);
epsilon   = 16*sqrt(2)/(5*pi);

% linear part
A = [-C 0 gammaastm(1) 0 0 0
    0 -C betam(1) 0 0 0
    -gammam(1) -betam(1) -C 0 0 0
    0 0 0 -C 0 gammaastm(2)
    0 0 0 0 -C betam(2)
    0 0 0 -gammam(2) -betam(2) -C];
B = eye(6);

% quadratic part
F2 = sptensor([6 6 6]);
F2(2,1,3) = -alpham(1);
F2(2,4,6) = -deltam(1);
F2(3,1,2) = alpham(1);
F2(3,4,5) = deltam(1);
F2(4,2,6) = epsilon;
F2(4,3,5) = -epsilon;
F2(5,1,6) = -alpham(2);
F2(5,3,4) = -deltam(2);
F2(6,1,5) = alpham(2);
F2(6,2,4) = deltam(2);
F = {F2};

% shift of origin
if numel(varargin)>0
   sol = varargin{1};
   x = sol.x;
   A = A+spmatrix(ttv(F2,x,2)+ttv(F2,x,3));
end

end