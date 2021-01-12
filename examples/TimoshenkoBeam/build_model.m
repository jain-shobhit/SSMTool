function [M,C,K,fnl,fext, outdof] = build_model(nElements)

% Geometrically Nonlinear Timoshenko beam 
n = 5*nElements + 1;
forcing_dof = n - 1;
outdof = forcing_dof;
%% 
% construct forcing amplitude vector

P = 120 * 1e7; % forcing amplitude vector
fext = sparse(n,1);
fext(forcing_dof) = P;

[~,K,f,x,xd,C,M,~]=FEM_Timoshenko(nElements,0);

M = sparse(double(M)); %#ok<*NODEF>
C = sparse(double(C));
K = sparse(double(K));
ndof = length(x);

%% get polynomial stiffness coefficients 

f = simplify(subs(f,xd,sym(zeros(ndof,1))));
[powers,coefficients] = get_coefficients(f,x);

% Obtaining tensor from multi-index coefficients 
degree = sum(powers,2); % get degree associated to each multi-index
%% 
% extract multi-indices and tensors corresponding to degree 2

idx_2 = find(degree==2);
If2 = powers(idx_2,:);
Cf2 = coefficients(idx_2,:);
f2 = multi_index_to_tensor(Cf2.',If2);

%% 
% extract multi-indices and tensors corresponding to degree 3
idx_3 = find(degree==3);
If3 = powers(idx_3,:);
Cf3 = coefficients(idx_3,:);
f3 = multi_index_to_tensor(Cf3.',If3);

fnl = {f2,f3};
end

function [E,P] = get_coefficients(f,x)
% this function returns the nonlinear coefficients required in the format
% for Harmonic Balance using NLvib.
ndof = length(x);
E = [];
P = [];

for i=1:ndof
    [c,t] = coeffs(f(i),x.');
    for j = 1:length(c)
        exponent = get_exponent(t(j),x);
        coefficient = double(c(j));
        
        if isempty(E) % taking care of the trivial case
            is = false;
        else
            [is, loc] = ismember(exponent,E,'rows');
        end
        
        if is
            P(loc,i) = coefficient;
        else
            p = zeros(1,ndof);
            p(i) = coefficient;
            P = [P; p];
            E = [E;exponent];
        end        
    end    
end
E(1,:) = [];
P(1,:) = [];


    function e = get_exponent(monomial,x)          
        factors = factor(monomial);
        e = sum(double(jacobian(factors.',x)));        
    end
end

