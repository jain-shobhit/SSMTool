function [M,C,K,fnl,fext,outdof] = build_model(nElements,isViscoelastic)

if isViscoelastic
    name = append('Vis_Damp_Timo_model_nE=',num2str(nElements),'.mat');
else
    name = append('Lin_Damp_Timo_model_nE=',num2str(nElements),'.mat');
end
try
    load(name,'M','C','K','fnl','fext','outdof');
    fnl = {sptensor(double(fnl{1})),sptensor(double(fnl{2}))};
catch
    % Geometrically Nonlinear Timoshenko beam
    n = 5*nElements + 1;
    forcing_dof = n - 1;
    outdof = forcing_dof;
    %%
    % construct forcing amplitude vector
    P = 120 * 1e7; % moment amplitude
    f_0 = sparse(n,1);
    f_0(forcing_dof) = P;
    fext = f_0;

    [~,K,f,x,xd,C,M,~]=FEM_Timoshenko(nElements,0);

    M = sparse(double(M)); %#ok<*NODEF>
    C = sparse(double(C));
    K = sparse(double(K));
    ndof = length(x);
    %% Additional mass at 0.25L
    mass = 0;
    moment_of_inertia = 0;
    % determine DOF
    node_idx = round(nElements/4)+1;
    mDOFs = 5*(node_idx-1)-3+[1, 4]; % translational DOFs
    rDOF = 5*(node_idx-1)-3+3;       % rotational DOF
    M(mDOFs,mDOFs) = M(mDOFs,mDOFs) + mass*eye(2,2); % adding mass to translational DOF
    M(rDOF,rDOF)   = M(rDOF,rDOF)+moment_of_inertia;

    %% get polynomial stiffness coefficients
    if isViscoelastic
        [powers,coefficients] = get_coefficients(f,[x;xd]);
    else
        f = simplify(subs(f,xd,sym(zeros(ndof,1))));
        [powers,coefficients] = get_coefficients(f,x);
    end

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

    save(name,'M','C','K','fnl','fext','outdof','f');
end
end

function [E,P] = get_coefficients(f,x)
% this function returns the nonlinear coefficients required in the format
% for Harmonic Balance using NLvib.
n = length(f);
E = [];
P = [];

for i=1:n
    [c,t] = coeffs(f(i),x.');
    c = double(c);
    for j = 1:length(c)
        exponent = get_exponent(t(j),x);
        coefficient = c(j);

        if isempty(E) % taking care of the trivial case
            is = false;
        else
            [is, loc] = ismember(exponent,E,'rows');
        end

        if is
            P(loc,i) = coefficient;
        else
            p = zeros(1,n);
            p(i) = coefficient;
            P = [P; p];
            E = [E;exponent];
        end
    end
end


    function e = get_exponent(monomial,x)
        factors = factor(monomial);
        e = sum(double(jacobian(factors.',x)));
    end
end

