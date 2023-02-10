function y = lambda_constraints(t,x,p,type)

q1 = x(1,:);
q2 = x(2,:);
q3 = x(3,:);
v1 = x(4,:);
v2 = x(5,:);
v3 = x(6,:);
epf = p(1,:);
om  = p(2,:);

om1 = 2;
om2 = 3;
om3 = 5;
zeta1 = 0.01;
zeta2 = 0.05;
zeta3 = 0.10;

alpha = 1;
beta  = 1;

omav = 0.5*(om1^2+om2^2+om3^2);

F1 = 2*zeta1*om1*v1+om1^2*q1+0.5*om1^2*(3*q1.^2+q2.^2+q3.^2)+om2^2*q1.*q2+...
    +om3^2*q1.*q3+omav*q1.*(q1.^2+q2.^2+q3.^2);
F2 = 2*zeta2*om2*v2+om2^2*q2+0.5*om2^2*(3*q2.^2+q1.^2+q3.^2)+om1^2*q1.*q2+...
    +om3^2*q2.*q3+omav*q2.*(q1.^2+q2.^2+q3.^2);
F3 = 2*zeta3*om3*v3+om3^2*q3+0.5*om3^2*(3*q3.^2+q1.^2+q2.^2)+om1^2*q1.*q3+...
    +om2^2*q2.*q3+omav*q3.*(q1.^2+q2.^2+q3.^2);

F1 = epf.*cos(om.*t)-F1;
F2 = -F2;
F3 = -F3;

switch type
    case 'cubic'
        g  = q3-q1.^3-q2.^3;
        dg = v3-3*q1.^2.*v1-3*q2.^2.*v2;
        c  = 6*q1.*v1.^2+6*q2.*v2.^2-alpha*g-beta*dg; % (alpha,beta) for stabilization
        y   = -(c+3*q1.^2.*F1+3*q2.^2.*F2-F3)./(9*q1.^4+9*q2.^4+1);
        
    case 'sphere'
        g  = q1.^2+q2.^2+(q3-1).^2-1;
        dg = 2*q1.*v1+2*q2.*v2+2*(q3-1).*v3;
        c  = -2*(v1.^2+v2.^2+v3.^2)-alpha*g-beta*dg; % (alpha,beta) for stabilization
        y   = -(c-2*q1.*F1-2*q2.*F2-2*(q3-1).*F3)/4;    
        
    otherwise
        error('type should be either cubic or sphere');       
end

end
