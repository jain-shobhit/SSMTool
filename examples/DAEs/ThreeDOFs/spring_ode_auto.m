function y = spring_ode_auto(x,p,type)

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
v1 = x(4,:);
v2 = x(5,:);
v3 = x(6,:);
zeta = p(1,:);

om1 = 2;
om2 = 3;
om3 = 5;

alpha = 1;
beta  = 1;

omav = 0.5*(om1^2+om2^2+om3^2);

F1 = 2*zeta.*om1.*v1+om1^2*x1+0.5*om1^2*(3*x1.^2+x2.^2+x3.^2)+om2^2*x1.*x2+...
    +om3^2*x1.*x3+omav*x1.*(x1.^2+x2.^2+x3.^2);
F2 = 2*zeta.*om2.*v2+om2^2*x2+0.5*om2^2*(3*x2.^2+x1.^2+x3.^2)+om1^2*x1.*x2+...
    +om3^2*x2.*x3+omav*x2.*(x1.^2+x2.^2+x3.^2);
F3 = 2*zeta.*om3.*v3+om3^2*x3+0.5*om3^2*(3*x3.^2+x1.^2+x2.^2)+om1^2*x1.*x3+...
    +om2^2*x2.*x3+omav*x3.*(x1.^2+x2.^2+x3.^2);

F1 = -F1;
F2 = -F2;
F3 = -F3;

switch type
    case 'cubic'
        g  = x3-x1.^3-x2.^3;
        dg = v3-3*x1.^2.*v1-3*x2.^2.*v2;
        c  = 6*x1.*v1.^2+6*x2.*v2.^2-alpha*g-beta*dg; % (alpha,beta) for stabilization
        temp   = (c+3*x1.^2.*F1+3*x2.^2.*F2-F3)./(9*x1.^4+9*x2.^4+1);
        y(1,:) = v1;
        y(2,:) = v2;
        y(3,:) = v3;
        y(4,:) = F1-3*x1.^2.*temp;
        y(5,:) = F2-3*x2.^2.*temp;
        y(6,:) = F3+temp;
        
    case 'sphere'
        g  = x1.^2+x2.^2+(x3-1).^2-1;
        dg = 2*x1.*v1+2*x2.*v2+2*(x3-1).*v3;
        c  = -2*(v1.^2+v2.^2+v3.^2)-alpha*g-beta*dg; % (alpha,beta) for stabilization
        temp   = (c-2*x1.*F1-2*x2.*F2-2*(x3-1).*F3)/4;
        y(1,:) = v1;
        y(2,:) = v2;
        y(3,:) = v3;
        y(4,:) = F1+2*x1.*temp;
        y(5,:) = F2+2*x2.*temp;
        y(6,:) = F3+2*(x3-1).*temp;    
        
    case 'none'
        y(1,:) = v1;
        y(2,:) = v2;
        y(3,:) = v3;
        y(4,:) = F1;
        y(5,:) = F2;
        y(6,:) = F3;          
        
    otherwise
        error('type should be cubic/sphere/none');       
end

end
